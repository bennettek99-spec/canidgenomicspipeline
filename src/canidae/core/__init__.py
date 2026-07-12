"""Core foundation layer for CANIS.

Everything in :mod:`canidae.core` is pure-Python, depends on no pipeline stage, and is the
substrate every stage builds on: config, logging, provenance, the domain model, the
datastore, tool runners, the plugin registry, the Stage contract, and the executor.
"""

from __future__ import annotations

from canidae.core.config import GlobalConfig
from canidae.core.datastore import DataStore
from canidae.core.executor import ExecutionReport, NativeExecutor, build_context
from canidae.core.model import (
    AnalysisResult,
    Artifact,
    ArtifactKind,
    Callset,
    CanidTaxon,
    Cohort,
    Dataset,
    FileFormat,
    GenomicInterval,
    Individual,
    ReferenceGenome,
    Sample,
    Taxon,
)
from canidae.core.provenance import ProvenanceRecord, ProvenanceWriter
from canidae.core.registry import ANALYSES, DATASETS, REFERENCES, STAGES
from canidae.core.runtime import LocalRunner, ResourceSpec, ToolSpec, make_runner
from canidae.core.stage import ArtifactSpec, RunContext, Stage, StageConfig, StageResult

__all__ = [
    "ANALYSES",
    "DATASETS",
    "REFERENCES",
    "STAGES",
    "AnalysisResult",
    "Artifact",
    "ArtifactKind",
    "ArtifactSpec",
    "Callset",
    "CanidTaxon",
    "Cohort",
    "DataStore",
    "Dataset",
    "ExecutionReport",
    "FileFormat",
    "GenomicInterval",
    "GlobalConfig",
    "Individual",
    "LocalRunner",
    "NativeExecutor",
    "ProvenanceRecord",
    "ProvenanceWriter",
    "ReferenceGenome",
    "ResourceSpec",
    "RunContext",
    "Sample",
    "Stage",
    "StageConfig",
    "StageResult",
    "Taxon",
    "ToolSpec",
    "build_context",
    "make_runner",
]
