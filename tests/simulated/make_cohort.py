"""Simulate a small canid cohort with *known* structure for validation.

The demography encodes a ground truth we can assert against:

    * DOG splits from WOLF recently  -> low F_ST(wolf, dog)
    * COYOTE diverges deeply from the wolf/dog ancestor -> high F_ST(coyote, *)

So any correct pipeline must recover  F_ST(coyote,wolf) > F_ST(wolf,dog)  and separate
coyote from wolf+dog on PC1. Used by the Phase-1 integration test; also runnable as a
script to drop a VCF + sample sheet somewhere for manual inspection.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

import msprime

# taxon label for each simulated population
_TAXON = {"WOLF": "gray_wolf", "DOG": "domestic_dog", "COYOTE": "coyote"}
# approximate localities (for the future geographic module)
_LOCALITY = {
    "WOLF": (61.0, 105.0, "eurasia"),
    "DOG": (52.0, 13.0, "europe"),
    "COYOTE": (39.0, -98.0, "north_america"),
}


def simulate_cohort(
    out_dir: Path,
    *,
    n_per_pop: int = 6,
    seed: int = 42,
    sequence_length: float = 2_000_000,
) -> tuple[Path, Path]:
    """Write ``cohort.vcf`` and ``samples.csv`` into ``out_dir``; return their paths."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dem = msprime.Demography()
    dem.add_population(name="WOLF", initial_size=10_000)
    dem.add_population(name="DOG", initial_size=8_000)
    dem.add_population(name="COYOTE", initial_size=12_000)
    dem.add_population(name="WOLFDOG", initial_size=10_000)
    dem.add_population(name="ANC", initial_size=15_000)
    dem.add_population_split(time=2_000, derived=["WOLF", "DOG"], ancestral="WOLFDOG")
    dem.add_population_split(time=15_000, derived=["WOLFDOG", "COYOTE"], ancestral="ANC")

    ts = msprime.sim_ancestry(
        samples={"WOLF": n_per_pop, "DOG": n_per_pop, "COYOTE": n_per_pop},
        demography=dem,
        sequence_length=sequence_length,
        recombination_rate=1e-8,
        random_seed=seed,
    )
    mts = msprime.sim_mutations(ts, rate=1.5e-8, random_seed=seed)

    pop_name = {p.id: p.metadata.get("name", f"pop{p.id}") for p in mts.populations()}
    names: list[str] = []
    counts: dict[str, int] = {}
    for ind in mts.individuals():
        pname = pop_name[mts.node(ind.nodes[0]).population]
        counts[pname] = counts.get(pname, 0) + 1
        names.append(f"{pname}{counts[pname]:02d}")

    vcf_path = out_dir / "cohort.vcf"
    with vcf_path.open("w", encoding="utf-8") as fh:
        mts.write_vcf(fh, individual_names=names, contig_id="1")

    jitter = random.Random(seed)
    sheet_path = out_dir / "samples.csv"
    with sheet_path.open("w", encoding="utf-8") as fh:
        fh.write("sample_id,taxon,population,region,latitude,longitude\n")
        for name in names:
            pop = name[:-2]  # strip the two-digit suffix
            lat, lon, region = _LOCALITY[pop]
            # small per-sample jitter so localities vary within a population
            jlat = round(lat + jitter.uniform(-3.0, 3.0), 4)
            jlon = round(lon + jitter.uniform(-3.0, 3.0), 4)
            fh.write(f"{name},{_TAXON[pop]},{pop.lower()},{region},{jlat},{jlon}\n")

    return vcf_path, sheet_path


# Four-population setup for introgression validation: topology ((WOLF,DOG),COYOTE) with
# JACKAL as outgroup, and (optionally) a COYOTE->DOG gene-flow pulse.
_INTRO_TAXON = {
    "WOLF": "gray_wolf",
    "DOG": "domestic_dog",
    "COYOTE": "coyote",
    "JACKAL": "golden_jackal",
}
_INTRO_LOCALITY = {
    "WOLF": (61.0, 105.0, "eurasia"),
    "DOG": (52.0, 13.0, "europe"),
    "COYOTE": (39.0, -98.0, "north_america"),
    "JACKAL": (22.0, 78.0, "south_asia"),
}


def simulate_introgression_cohort(
    out_dir: Path,
    *,
    n_per_pop: int = 6,
    seed: int = 11,
    admixture_proportion: float = 0.18,
    sequence_length: float = 5_000_000,
) -> tuple[Path, Path]:
    """Write a 4-population cohort with a known COYOTE->DOG gene-flow pulse (set
    ``admixture_proportion=0`` for a no-flow control). Returns (vcf, sample_sheet)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    dem = msprime.Demography()
    for name in ("WOLF", "DOG", "COYOTE", "JACKAL", "WOLFDOG", "INGROUP", "ROOT"):
        dem.add_population(name=name, initial_size=10_000)
    dem.add_population_split(time=1_500, derived=["WOLF", "DOG"], ancestral="WOLFDOG")
    dem.add_population_split(time=8_000, derived=["WOLFDOG", "COYOTE"], ancestral="INGROUP")
    dem.add_population_split(time=25_000, derived=["INGROUP", "JACKAL"], ancestral="ROOT")
    if admixture_proportion > 0:
        # backward-in-time mass migration DOG->COYOTE == forward gene flow COYOTE->DOG
        dem.add_mass_migration(
            time=800, source="DOG", dest="COYOTE", proportion=admixture_proportion
        )
    dem.sort_events()

    ts = msprime.sim_ancestry(
        samples={"WOLF": n_per_pop, "DOG": n_per_pop, "COYOTE": n_per_pop, "JACKAL": n_per_pop},
        demography=dem,
        sequence_length=sequence_length,
        recombination_rate=1e-8,
        random_seed=seed,
    )
    mts = msprime.sim_mutations(ts, rate=1.5e-8, random_seed=seed)

    pop_name = {p.id: p.metadata.get("name", f"pop{p.id}") for p in mts.populations()}
    names: list[str] = []
    counts: dict[str, int] = {}
    for ind in mts.individuals():
        pname = pop_name[mts.node(ind.nodes[0]).population]
        counts[pname] = counts.get(pname, 0) + 1
        names.append(f"{pname}{counts[pname]:02d}")

    vcf_path = out_dir / "cohort.vcf"
    with vcf_path.open("w", encoding="utf-8") as fh:
        mts.write_vcf(fh, individual_names=names, contig_id="1")

    jitter = random.Random(seed)
    sheet_path = out_dir / "samples.csv"
    with sheet_path.open("w", encoding="utf-8") as fh:
        fh.write("sample_id,taxon,population,region,latitude,longitude\n")
        for name in names:
            pop = name[:-2]
            lat, lon, region = _INTRO_LOCALITY[pop]
            jlat = round(lat + jitter.uniform(-3.0, 3.0), 4)
            jlon = round(lon + jitter.uniform(-3.0, 3.0), 4)
            fh.write(f"{name},{_INTRO_TAXON[pop]},{pop.lower()},{region},{jlat},{jlon}\n")

    return vcf_path, sheet_path


if __name__ == "__main__":  # pragma: no cover
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("simulated_cohort")
    vcf, sheet = simulate_cohort(target)
    print(f"wrote {vcf} and {sheet}")
