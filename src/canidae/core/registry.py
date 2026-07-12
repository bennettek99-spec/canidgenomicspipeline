"""Registries for the pluggable pieces of CANIS.

Stages, dataset fetchers, reference genomes, and analyses all self-register here. New
modules become available simply by importing their package (which triggers the decorator)
or via Python entry points — no central list to edit. This is the mechanism that lets the
platform grow without restructuring.
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import Generic, TypeVar

from canidae.core.errors import RegistryError
from canidae.core.logging import get_logger

_log = get_logger("registry")

T = TypeVar("T")


class Registry(Generic[T]):
    """A named registry with decorator-based registration."""

    def __init__(self, kind: str) -> None:
        self.kind = kind
        self._items: dict[str, T] = {}

    def register(self, name: str, *, override: bool = False):
        """Decorator: ``@registry.register("qc")`` binds a class/factory under ``name``."""

        def _decorator(obj: T) -> T:
            if name in self._items and not override:
                raise RegistryError(
                    f"{self.kind} '{name}' already registered; pass override=True to replace"
                )
            self._items[name] = obj
            _log.debug("registered %s '%s' -> %r", self.kind, name, obj)
            return obj

        return _decorator

    def add(self, name: str, obj: T, *, override: bool = False) -> None:
        self.register(name, override=override)(obj)

    def get(self, name: str) -> T:
        try:
            return self._items[name]
        except KeyError:
            raise RegistryError(
                f"unknown {self.kind} '{name}'. Registered: {sorted(self._items)}"
            ) from None

    def names(self) -> list[str]:
        return sorted(self._items)

    def __contains__(self, name: str) -> bool:
        return name in self._items

    def __iter__(self) -> Iterator[str]:
        return iter(self.names())

    def __len__(self) -> int:
        return len(self._items)


# The canonical registries. Import and decorate against these.
STAGES: Registry = Registry("stage")
DATASETS: Registry = Registry("dataset_fetcher")
REFERENCES: Registry = Registry("reference")
ANALYSES: Registry = Registry("analysis")


def load_entry_point_plugins(group: str = "canidae.stages") -> int:
    """Discover and import third-party plugins advertised via entry points.

    Returns the number of entry points loaded. Safe to call multiple times.
    """
    try:
        from importlib.metadata import entry_points
    except ImportError:  # pragma: no cover
        return 0
    count = 0
    for ep in entry_points(group=group):
        try:
            ep.load()  # importing the module runs its registration decorators
            count += 1
        except Exception as exc:  # pragma: no cover - defensive
            _log.warning("failed to load plugin %s: %s", ep.name, exc)
    return count
