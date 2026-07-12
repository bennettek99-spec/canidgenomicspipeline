"""Command-line interface for CANIS.

Thin, declarative Typer app. Every command resolves a :class:`GlobalConfig` from layered
YAML plus dotted ``--set key=value`` overrides, then delegates to the library. The CLI adds
no analysis logic of its own — it is one of several front-ends (Python API, notebooks) onto
the same core.
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from canidae.core.config import GlobalConfig
from canidae.core.registry import STAGES, load_entry_point_plugins
from canidae.version import __version__

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="CANIS — canid evolutionary genomics pipeline.",
)
console = Console()


def _load_config(
    config_files: list[Path] | None,
    overrides: list[str] | None,
) -> GlobalConfig:
    parsed: dict[str, str] = {}
    for item in overrides or []:
        if "=" not in item:
            raise typer.BadParameter(f"override must be key=value, got '{item}'")
        key, value = item.split("=", 1)
        parsed[key.strip()] = value
    return GlobalConfig.load(*(config_files or []), overrides=parsed)


@app.command()
def version() -> None:
    """Print the CANIS version."""
    console.print(f"CANIS (canidae) [bold]{__version__}[/bold]")


@app.command("config")
def show_config(
    config: list[Path] | None = typer.Option(None, "--config", "-c",
                                                 help="YAML config file(s), merged in order."),
    set_: list[str] | None = typer.Option(None, "--set", "-s",
                                             help="Dotted override, e.g. executor.max_workers=8."),
) -> None:
    """Resolve and print the effective configuration."""
    cfg = _load_config(config, set_)
    console.print(f"[dim]config digest:[/dim] {cfg.digest()[:16]}")
    console.print(cfg.to_yaml())


@app.command()
def stages() -> None:
    """List all registered pipeline stages."""
    from canidae.stages import load_builtin_stages

    load_builtin_stages()
    load_entry_point_plugins()
    table = Table(title="Registered stages")
    table.add_column("name", style="bold")
    table.add_column("config model")
    for name in STAGES.names():
        cls = STAGES.get(name)
        table.add_row(name, getattr(cls.config_model, "__name__", "?"))
    if not len(STAGES):
        console.print("[yellow]No stages registered yet (foundation-only build).[/yellow]")
    else:
        console.print(table)


@app.command()
def run(
    config: list[Path] | None = typer.Option(None, "--config", "-c",
                                                 help="YAML config file(s), merged in order."),
    set_: list[str] | None = typer.Option(None, "--set", "-s",
                                             help="Dotted override, e.g. seed=7."),
    dry_run: bool = typer.Option(False, "--dry-run",
                                 help="Assemble and print the DAG without executing."),
    run_id: str | None = typer.Option(None, "--run-id", help="Explicit run directory name."),
) -> None:
    """Run the configured pipeline (or preview its DAG with --dry-run)."""
    from canidae.pipeline import run_pipeline  # local import keeps CLI startup fast

    cfg = _load_config(config, set_)
    if not cfg.pipeline:
        console.print("[yellow]config.pipeline is empty — nothing to run.[/yellow]")
        raise typer.Exit(code=0)

    report = run_pipeline(cfg, dry_run=dry_run, run_id=run_id)
    verb = "planned" if dry_run else "executed"
    console.print(f"[green]{verb}:[/green] {', '.join(report.executed) or '(none)'}")
    if report.skipped:
        console.print(f"[dim]skipped:[/dim] {', '.join(report.skipped)}")
    if report.paused:
        console.print(f"[yellow]paused safely:[/yellow] {', '.join(report.paused)}")
    if report.cancelled:
        console.print(f"[yellow]stopped safely:[/yellow] {', '.join(report.cancelled)}")
    if report.failed:
        for name, err in report.failed.items():
            console.print(f"[red]failed[/red] {name}: {err}")
        raise typer.Exit(code=1)


@app.command()
def ui(
    host: str = typer.Option("127.0.0.1", help="Loopback address for the local-only UI."),
    port: int = typer.Option(8765, min=0, max=65535, help="Local browser UI port (0 picks one)."),
    no_browser: bool = typer.Option(False, help="Start the UI without opening a browser tab."),
) -> None:
    """Open the local browser interface for laptop-safe pipeline runs."""
    from canidae.ui import serve_ui

    if host not in {"127.0.0.1", "localhost", "::1"}:
        raise typer.BadParameter("the local UI only accepts a loopback host")
    console.print("[green]Starting local-only CANIS UI…[/green]")
    serve_ui(Path.cwd(), host=host, port=port, open_browser=not no_browser)


if __name__ == "__main__":  # pragma: no cover
    app()
