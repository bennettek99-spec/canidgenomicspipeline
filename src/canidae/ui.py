# ruff: noqa: E501
"""Local, dependency-free browser interface for laptop-safe CANIS runs.

The UI intentionally binds only to loopback and keeps all state in the local project:
there is no account, cloud service, telemetry, or hosted server.  It is a thin controller
over the same validated config and pipeline used by ``canidae run``.
"""

from __future__ import annotations

import json
import threading
import webbrowser
from dataclasses import dataclass, field
from datetime import datetime
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

from canidae.core.config import GlobalConfig
from canidae.core.datastore import DataStore
from canidae.core.resources import ResourceManager
from canidae.pipeline import run_pipeline


class UiRequestError(ValueError):
    """A friendly validation error returned to the local browser UI."""


PANEL_ESTIMATES = {
    "2k": {"download": "under 1 GB", "disk": "2-3 GB", "runtime": "5-15 minutes"},
    "10k": {"download": "about 1-3 GB", "disk": "4-7 GB", "runtime": "15-45 minutes"},
    "25k": {"download": "about 2-8 GB", "disk": "8-16 GB", "runtime": "30-90 minutes"},
}


def estimate_laptop_run(preset: str) -> dict[str, str]:
    """Return a transparent preliminary estimate before remote byte-range preflight."""
    if preset not in PANEL_ESTIMATES:
        raise UiRequestError("Choose a reduced panel size of 2k, 10k, or 25k.")
    estimate = dict(PANEL_ESTIMATES[preset])
    estimate["note"] = (
        "The reduced-panel stage calculates the exact indexed-range estimate before downloading "
        "source VCF ranges and enforces its 9 GB ceiling."
    )
    return estimate


def build_ui_config(payload: dict[str, Any], project_root: Path) -> GlobalConfig:
    """Translate simple UI controls into the normal validated pipeline configuration."""
    mode = str(payload.get("dataset_mode", "local"))
    if mode not in {"local", "redwolf_public"}:
        raise UiRequestError("Choose either a local VCF or the public Red Wolf/jackal project.")
    analyses = {str(value) for value in payload.get("analyses", [])}
    preset = str(payload.get("panel_preset", "10k"))
    if mode == "redwolf_public":
        estimate_laptop_run(preset)

    workers = _bounded_int(payload.get("max_workers", 2), default=2, minimum=1, maximum=8)
    memory_mb = _bounded_int(payload.get("memory_mb", 8192), default=8192, minimum=1024, maximum=65536)
    pipeline = _analysis_pipeline(
        analyses, entry="reduced_panel" if mode == "redwolf_public" else "ingest"
    )
    overrides: dict[str, Any] = {
        "project_name": str(payload.get("project_name") or "canidae-ui-run"),
        "paths.root": str(project_root),
        "executor.max_workers": workers,
        "resource_manager.max_workers": workers,
        "resource_manager.max_memory_mb": memory_mb,
        "resource_manager.temperature_friendly": True,
        "resource_manager.max_threads_per_stage": min(workers, 4),
        "pipeline": pipeline,
        "stages.qc.min_sample_call_rate": float(payload.get("min_sample_call_rate", 0.7)),
        "stages.qc.min_site_call_rate": float(payload.get("min_site_call_rate", 0.8)),
        "stages.report.title": str(payload.get("report_title") or "CANIS laptop analysis"),
    }
    if mode == "local":
        vcf = str(payload.get("vcf_path", "")).strip()
        sheet = str(payload.get("sample_sheet", "")).strip()
        if not vcf or not sheet:
            raise UiRequestError("Local VCF mode needs both a VCF/BCF path and a sample-sheet CSV path.")
        overrides.update({
            "stages.ingest.callset": vcf,
            "stages.ingest.sample_sheet": sheet,
            "stages.ingest.dataset_id": str(payload.get("dataset_id") or "local_ui_dataset"),
            "stages.ingest.reference_build": str(payload.get("reference_build") or "unknown"),
            "stages.ingest.panel_relative": bool(payload.get("panel_relative", False)),
        })
    else:
        sheet = str(payload.get("sample_sheet", "")).strip()
        if not sheet:
            sheet = "configs/examples/redwolf_jackal_aadr_samples.csv"
        selected = _sample_ids(str(payload.get("selected_samples", "")))
        overrides.update({
            "stages.reduced_panel.sample_sheet": sheet,
            "stages.reduced_panel.preset": preset,
            "stages.reduced_panel.selected_samples": selected,
            "stages.reduced_panel.confirm_large_transfer": bool(payload.get("confirm_large_transfer", False)),
            "stages.reduced_panel.dataset_id": "PRJNA448733_reduced_panel",
            "stages.reduced_panel.reference_id": "CanFam3.1",
        })
    return GlobalConfig.load(overrides=overrides)


def _analysis_pipeline(analyses: set[str], *, entry: str) -> list[str]:
    """Report-compatible core metrics plus optional analyses selected in the browser."""
    # The report's core sections require PCA/FST/diversity, so those remain on to avoid a
    # broken package even when optional analyses are unchecked.
    pipeline = [entry, "qc", "load_genotypes"]
    if "tree" in analyses or "distance" in analyses:
        pipeline.append("distance")
    pipeline.extend(["pca", "fst", "diversity"])
    if "tree" in analyses:
        pipeline.append("nj_tree")
    if "admixture" in analyses:
        pipeline.append("admixture")
    if "introgression" in analyses:
        pipeline.append("dstats")
    if "local_ancestry" in analyses:
        pipeline.append("local_ancestry")
    pipeline.append("report")
    return pipeline


def _sample_ids(text: str) -> list[str]:
    return [value.strip() for value in text.split(",") if value.strip()]


def _bounded_int(value: Any, *, default: int, minimum: int, maximum: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    return min(max(parsed, minimum), maximum)


@dataclass(slots=True)
class RunController:
    """Own one local run and expose cooperative pause/resume/cancel state to the UI."""

    project_root: Path
    state: str = "idle"
    message: str = "Choose a dataset and start a laptop-safe run."
    run_id: str | None = None
    config: GlobalConfig | None = None
    error: str | None = None
    outputs: list[dict[str, str]] = field(default_factory=list)
    _thread: threading.Thread | None = field(default=None, init=False, repr=False)
    _lock: threading.RLock = field(default_factory=threading.RLock, init=False, repr=False)

    def start(self, payload: dict[str, Any]) -> dict[str, Any]:
        with self._lock:
            if self._thread is not None and self._thread.is_alive():
                raise UiRequestError("A run is already active. Pause, stop, or wait for it to finish.")
            config = build_ui_config(payload, self.project_root)
            self.config = config
            self.run_id = f"ui-{datetime.now().strftime('%Y%m%d-%H%M%S')}"
            self.state = "running"
            self.message = "Pipeline started. Stages are promoted only after validation."
            self.error = None
            self.outputs = []
            self._thread = threading.Thread(target=self._run, name="canidae-ui-run", daemon=True)
            self._thread.start()
            return self.status()

    def pause(self) -> dict[str, Any]:
        manager = self._manager()
        manager.request_pause()
        with self._lock:
            self.state = "pause_requested"
            self.message = "Pause requested. The pipeline stops safely after the current stage batch."
        return self.status()

    def stop(self) -> dict[str, Any]:
        manager = self._manager()
        manager.request_cancel()
        with self._lock:
            self.state = "stop_requested"
            self.message = "Stop requested. Completed validated stages remain cached for recovery."
        return self.status()

    def resume(self) -> dict[str, Any]:
        with self._lock:
            if self.config is None or self.run_id is None:
                raise UiRequestError("There is no prior UI run to resume.")
            if self._thread is not None and self._thread.is_alive():
                raise UiRequestError("The current run is still active.")
            manager = self._manager()
            manager.resume()
            manager.clear_cancel()
            self.state = "running"
            self.message = "Resuming with the same run ID; unchanged stages will use the safe cache."
            self.error = None
            self._thread = threading.Thread(target=self._run, name="canidae-ui-resume", daemon=True)
            self._thread.start()
            return self.status()

    def status(self) -> dict[str, Any]:
        with self._lock:
            return {
                "state": self.state,
                "message": self.message,
                "run_id": self.run_id,
                "error": self.error,
                "outputs": list(self.outputs),
                "active": bool(self._thread is not None and self._thread.is_alive()),
            }

    def _manager(self) -> ResourceManager:
        if self.config is None or self.run_id is None:
            raise UiRequestError("Start a run before using run controls.")
        return ResourceManager(self.config.resource_manager, run_dir=self.config.paths.run_root / self.run_id)

    def _run(self) -> None:
        assert self.config is not None and self.run_id is not None
        try:
            report = run_pipeline(self.config, run_id=self.run_id)
            outputs = _outputs_for(self.config)
            with self._lock:
                self.outputs = outputs
                if report.cancelled:
                    self.state = "cancelled"
                    self.message = "Stopped safely. Remove the stop marker or use Resume to continue."
                elif report.paused:
                    self.state = "paused"
                    self.message = "Paused safely. Use Resume when ready."
                elif report.failed:
                    self.state = "failed"
                    self.error = "; ".join(f"{stage}: {error}" for stage, error in report.failed.items())
                    self.message = "The manifest contains the error and a recovery instruction."
                else:
                    self.state = "completed"
                    self.message = "Completed. Open the self-contained report or any companion output below."
        except Exception as exc:
            with self._lock:
                self.state = "failed"
                self.error = f"{type(exc).__name__}: {exc}"
                self.message = "The run stopped before completion. Check the manifest and correct the plain-language error."


def _outputs_for(config: GlobalConfig) -> list[dict[str, str]]:
    store = DataStore(config.paths.data_root / "store")
    return [
        {"label": f"{artifact.role}: {artifact.path.name}", "path": str(artifact.path)}
        for artifact in sorted(store.all(), key=lambda artifact: (artifact.kind.value, artifact.role))
        if artifact.path.exists()
    ]


class _UiHandler(BaseHTTPRequestHandler):
    controller: RunController

    def do_GET(self) -> None:
        if self.path == "/":
            self._send_html(_HTML)
        elif self.path == "/api/status":
            self._send_json(self.controller.status())
        elif self.path.startswith("/api/estimate"):
            preset = self.path.split("preset=", 1)[-1] if "preset=" in self.path else "10k"
            self._send_json(estimate_laptop_run(preset))
        else:
            self._send_json({"error": "Not found"}, status=HTTPStatus.NOT_FOUND)

    def do_POST(self) -> None:
        try:
            payload = self._body()
            if self.path == "/api/start":
                response = self.controller.start(payload)
            elif self.path == "/api/pause":
                response = self.controller.pause()
            elif self.path == "/api/resume":
                response = self.controller.resume()
            elif self.path == "/api/stop":
                response = self.controller.stop()
            else:
                self._send_json({"error": "Not found"}, status=HTTPStatus.NOT_FOUND)
                return
            self._send_json(response)
        except UiRequestError as exc:
            self._send_json({"error": str(exc)}, status=HTTPStatus.BAD_REQUEST)
        except Exception as exc:  # defensive boundary: never show a raw traceback in the UI
            self._send_json({"error": f"Unexpected local UI error: {type(exc).__name__}: {exc}"},
                            status=HTTPStatus.INTERNAL_SERVER_ERROR)

    def log_message(self, _format: str, *_args: Any) -> None:
        return

    def _body(self) -> dict[str, Any]:
        length = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(length) if length else b"{}"
        value = json.loads(raw.decode("utf-8"))
        if not isinstance(value, dict):
            raise UiRequestError("The UI request must be an object.")
        return value

    def _send_json(self, payload: dict[str, Any], *, status: HTTPStatus = HTTPStatus.OK) -> None:
        body = json.dumps(payload, default=str).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_html(self, content: str) -> None:
        body = content.encode("utf-8")
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


def serve_ui(project_root: Path, *, host: str = "127.0.0.1", port: int = 8765,
             open_browser: bool = True) -> None:
    """Launch the browser UI and serve until the user interrupts it with Ctrl+C."""
    controller = RunController(Path(project_root).resolve())

    class Handler(_UiHandler):
        pass

    Handler.controller = controller
    server = ThreadingHTTPServer((host, port), Handler)
    url = f"http://{host}:{server.server_port}"
    if open_browser:
        webbrowser.open(url)
    print(f"CANIS UI is running at {url} (Ctrl+C to stop the local UI server).")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nCANIS UI server stopped.")
    finally:
        server.server_close()


_HTML = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CANIS laptop genomics</title><style>
:root{font-family:system-ui,-apple-system,Segoe UI,sans-serif;color:#172033;background:#f6f8fb}body{max-width:1050px;margin:24px auto;padding:0 18px}h1{margin-bottom:4px}.sub{color:#536071}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(270px,1fr));gap:14px}.card{background:white;border:1px solid #dde3eb;border-radius:12px;padding:16px;box-shadow:0 1px 2px #00000008}label{display:block;font-size:.87rem;font-weight:650;margin:10px 0 4px}input,select,textarea{box-sizing:border-box;width:100%;padding:8px;border:1px solid #b8c5d6;border-radius:6px;font:inherit}textarea{min-height:60px}button{margin:8px 7px 0 0;border:0;border-radius:7px;padding:9px 13px;background:#006da8;color:white;font-weight:650;cursor:pointer}button.alt{background:#4b5563}button.warn{background:#b45309}.badge{display:inline-block;padding:3px 8px;border-radius:99px;background:#e0f2fe;color:#075985;font-size:.8rem;font-weight:700}.estimate{background:#f8fafc;border-left:4px solid #0ea5e9;padding:10px;margin-top:10px}.status{white-space:pre-wrap}.outputs a{display:block;margin:.35rem 0;color:#075985}.check{display:flex;align-items:center;gap:7px;font-weight:500}.check input{width:auto}.note{font-size:.84rem;color:#536071}</style></head>
<body><h1>CANIS laptop genomics</h1><p class="sub">Local-only control panel. No account, hosted service, or full-source VCF download.</p>
<div class="grid"><section class="card"><h2>Dataset</h2><label>Source</label><select id="mode" onchange="updateMode()"><option value="local">My local pre-called VCF</option><option value="redwolf_public">Public Red Wolf / jackal indexed VCF</option></select><label>VCF / BCF path (local mode)</label><input id="vcf" placeholder="C:\\data\\canids.vcf.gz"><label>Sample-sheet CSV</label><input id="sheet" placeholder="configs/examples/redwolf_jackal_aadr_samples.csv"><label>Reference build</label><input id="reference" value="CanFam3.1"><label>Dataset label</label><input id="dataset" value="laptop_dataset"><label class="check"><input id="panelRelative" type="checkbox"> This local VCF is a fixed SNP panel (no callable-genome denominator)</label><label>Selected samples (public mode, comma-separated)</label><textarea id="samples" placeholder="Wolf25, Wolf26, GoldenJackal01"></textarea></section>
<section class="card"><h2>Laptop-safe plan</h2><label>Reduced panel (public mode)</label><select id="preset" onchange="estimate()"><option value="2k">2,000 SNPs</option><option value="10k" selected>10,000 SNPs</option><option value="25k">25,000 SNPs</option></select><div id="estimate" class="estimate">Loading estimate…</div><label class="check"><input id="confirm" type="checkbox"> I confirm a transfer above the displayed confirmation threshold, if the exact preflight asks.</label><label>Maximum workers</label><input id="workers" type="number" min="1" max="8" value="2"><label>Memory ceiling (MiB)</label><input id="memory" type="number" min="1024" max="65536" value="8192"><p class="note">The exact indexed-range estimate runs before source VCF ranges are fetched. The hard ceiling is below 10 GB.</p></section>
<section class="card"><h2>Analysis package</h2><p class="note">PCA, F<sub>ST</sub>, and panel-aware diversity are included so the final report remains complete.</p><label class="check"><input type="checkbox" value="distance" checked> Pairwise distance</label><label class="check"><input type="checkbox" value="tree" checked> Neighbor-joining tree</label><label class="check"><input type="checkbox" value="admixture"> Admixture</label><label class="check"><input type="checkbox" value="introgression"> Chromosome-aware D-statistics</label><label class="check"><input type="checkbox" value="local_ancestry"> Chromosome-reset local ancestry</label><label>Report title</label><input id="title" value="CANIS laptop analysis"></section></div>
<section class="card"><h2>Run controls</h2><button onclick="start()">Start safe run</button><button class="alt" onclick="pause()">Pause</button><button class="alt" onclick="resume()">Resume</button><button class="warn" onclick="stopRun()">Stop safely</button><div id="status" class="status estimate">Idle.</div><div id="outputs" class="outputs"></div></section>
<script>const $=id=>document.getElementById(id);async function api(url,method='GET',data={}){const r=await fetch(url,{method,headers:{'Content-Type':'application/json'},body:method==='GET'?undefined:JSON.stringify(data)});const v=await r.json();if(!r.ok)throw Error(v.error||'Request failed');return v}function updateMode(){const pub=$('mode').value==='redwolf_public';$('vcf').disabled=pub;$('preset').disabled=!pub;$('samples').disabled=!pub;$('panelRelative').disabled=pub}async function estimate(){try{const v=await api('/api/estimate?preset='+$('preset').value);$('estimate').textContent='Preliminary transfer: '+v.download+' • disk: '+v.disk+' • runtime: '+v.runtime+'\n'+v.note}catch(e){$('estimate').textContent=e.message}}function payload(){return{dataset_mode:$('mode').value,vcf_path:$('vcf').value,sample_sheet:$('sheet').value,reference_build:$('reference').value,dataset_id:$('dataset').value,panel_relative:$('panelRelative').checked,selected_samples:$('samples').value,panel_preset:$('preset').value,confirm_large_transfer:$('confirm').checked,max_workers:$('workers').value,memory_mb:$('memory').value,report_title:$('title').value,analyses:[...document.querySelectorAll('input[type=checkbox][value]:checked')].map(x=>x.value)}}function render(v){$('status').textContent=(v.state||'')+'\n'+(v.message||'')+(v.error?'\nError: '+v.error:'')+ (v.run_id?'\nRun ID: '+v.run_id:'');$('outputs').innerHTML=(v.outputs||[]).map(o=>'<a href="file:///'+o.path.replaceAll('\\\\','/')+'" target="_blank">'+o.label+'</a>').join('')}async function start(){try{render(await api('/api/start','POST',payload()))}catch(e){$('status').textContent=e.message}}async function pause(){try{render(await api('/api/pause','POST'))}catch(e){$('status').textContent=e.message}}async function resume(){try{render(await api('/api/resume','POST'))}catch(e){$('status').textContent=e.message}}async function stopRun(){try{render(await api('/api/stop','POST'))}catch(e){$('status').textContent=e.message}}async function poll(){try{render(await api('/api/status'))}catch(_){}setTimeout(poll,1800)}updateMode();estimate();poll();</script></body></html>"""


__all__ = ["RunController", "UiRequestError", "build_ui_config", "estimate_laptop_run", "serve_ui"]
