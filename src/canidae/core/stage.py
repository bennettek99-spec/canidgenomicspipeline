"""The Stage contract — the single interface every pipeline module implements.

A stage declares the artifact *kinds/roles* it needs and produces (never file paths),
validates its inputs, and does its work given a :class:`RunContext`. Because stages depend
only on this contract and on the domain model, they are fully decoupled: adding a module is
writing a ``Stage`` subclass plus a config model and registering it.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, ClassVar

from pydantic import BaseModel

from canidae.core.config import GlobalConfig
from canidae.core.datastore import DataStore
from canidae.core.errors import StageInputError
from canidae.core.model import Artifact, ArtifactKind, Cohort
from canidae.core.provenance import ProvenanceRecord, ProvenanceWriter
from canidae.core.runtime import ResourceSpec, ToolRunner


class StageConfig(BaseModel):
    """Base class for per-stage configuration models. Stages subclass and add fields."""

    enabled: bool = True


@dataclass(frozen=True, slots=True)
class ArtifactSpec:
    """A declaration of an artifact a stage consumes or produces, by kind + role.

    ``role`` is a semantic label (e.g. ``"harmonized_callset"``). ``optional`` inputs do
    not block execution when absent.
    """

    kind: ArtifactKind
    role: str
    optional: bool = False

    @property
    def key(self) -> tuple[ArtifactKind, str]:
        return (self.kind, self.role)


@dataclass(slots=True)
class RunContext:
    """Everything a stage's :meth:`Stage.run` is given. The only channel into a stage —
    no globals, no cross-module imports."""

    config: GlobalConfig
    datastore: DataStore
    runner: ToolRunner
    provenance: ProvenanceWriter
    cohort: Cohort | None = None
    run_dir: Any = None  # pathlib.Path; kept loose to avoid import churn
    resource_manager: Any = None
    scratch: dict[str, Any] = field(default_factory=dict)

    def resources_for(self, resources: ResourceSpec | None) -> ResourceSpec:
        return resources or ResourceSpec.from_defaults(self.config.resources)


@dataclass(slots=True)
class StageResult:
    """What a stage returns: produced artifacts, its provenance record, and metrics."""

    artifacts: list[Artifact] = field(default_factory=list)
    provenance: ProvenanceRecord | None = None
    metrics: dict[str, Any] = field(default_factory=dict)
    skipped: bool = False


class Stage(ABC):
    """Abstract base for all pipeline stages."""

    #: Unique registry name, e.g. ``"qc"``. Set on each subclass.
    name: ClassVar[str] = ""
    #: The Pydantic model validating this stage's config block.
    config_model: ClassVar[type[StageConfig]] = StageConfig

    def __init__(self, config: StageConfig | None = None) -> None:
        self.config = config or self.config_model()

    # -- declarative interface (used for DAG assembly) ---------------------------------

    @abstractmethod
    def required_inputs(self) -> list[ArtifactSpec]:
        """Artifact kinds/roles this stage consumes."""

    @abstractmethod
    def produced_outputs(self) -> list[ArtifactSpec]:
        """Artifact kinds/roles this stage produces."""

    def resources(self) -> ResourceSpec | None:
        """Optional resource request; None means 'use the config default'."""
        return None

    # -- execution ---------------------------------------------------------------------

    def validate_inputs(self, ctx: RunContext) -> None:
        """Fail early if required (non-optional) inputs are missing before doing work."""
        missing = [
            spec
            for spec in self.required_inputs()
            if not spec.optional and not ctx.datastore.has(spec.kind, spec.role)
        ]
        if missing:
            pretty = ", ".join(f"{s.kind.value}:{s.role}" for s in missing)
            raise StageInputError(f"stage '{self.name}' missing inputs: {pretty}")

    def gather_inputs(self, ctx: RunContext) -> list[Artifact]:
        """Resolve declared inputs to concrete artifacts (skipping absent optionals)."""
        resolved: list[Artifact] = []
        for spec in self.required_inputs():
            if ctx.datastore.has(spec.kind, spec.role):
                resolved.append(ctx.datastore.get(spec.kind, spec.role))
            elif not spec.optional:
                raise StageInputError(
                    f"stage '{self.name}' input {spec.kind.value}:{spec.role} unavailable"
                )
        return resolved

    def outputs_present(self, ctx: RunContext) -> bool:
        """True when every declared output is already registered — enables resume/skip."""
        specs = self.produced_outputs()
        return bool(specs) and all(ctx.datastore.has(s.kind, s.role) for s in specs)

    @abstractmethod
    def run(self, ctx: RunContext) -> StageResult:
        """Perform the work and return produced artifacts + provenance."""

    def __repr__(self) -> str:  # pragma: no cover
        return f"<Stage {self.name!r}>"
