from __future__ import annotations

import pytest

from canidae.core.model import (
    CanidTaxon,
    Cohort,
    GenomicInterval,
    GeoLocality,
    Individual,
    Sex,
    Taxon,
)


def _ind(id_: str, taxon: CanidTaxon, pop: str) -> Individual:
    return Individual(id=id_, taxon=Taxon(taxon), population=pop)


def test_geolocality_validates_coordinates() -> None:
    with pytest.raises(ValueError):
        GeoLocality(latitude=100.0, longitude=0.0)
    loc = GeoLocality(latitude=45.0, longitude=-90.0)
    assert loc.has_coordinates


def test_interval_region_string() -> None:
    assert GenomicInterval("chr1", 0, 1_000_000).to_region_string() == "chr1:1-1000000"
    assert GenomicInterval("chr2").to_region_string() == "chr2"
    with pytest.raises(ValueError):
        GenomicInterval("chr1", 100, 50)


def test_cohort_grouping() -> None:
    cohort = Cohort(
        id="c1",
        individuals=[
            _ind("w1", CanidTaxon.GRAY_WOLF, "eurasia"),
            _ind("w2", CanidTaxon.GRAY_WOLF, "eurasia"),
            _ind("c1", CanidTaxon.COYOTE, "north_america"),
        ],
    )
    assert cohort.size == 3
    assert set(cohort.populations) == {"eurasia", "north_america"}
    assert cohort.taxa[CanidTaxon.GRAY_WOLF] == 2


def test_cohort_subset() -> None:
    cohort = Cohort(id="c1", individuals=[
        _ind("a", CanidTaxon.DINGO, "au"),
        _ind("b", CanidTaxon.DINGO, "au"),
    ])
    sub = cohort.subset({"a"})
    assert sub.size == 1 and sub.individuals[0].id == "a"


def test_taxon_ingroup_flag() -> None:
    assert CanidTaxon.GRAY_WOLF.is_ingroup_canis
    assert not CanidTaxon.AFRICAN_WILD_DOG.is_ingroup_canis


def test_sex_default_unknown() -> None:
    ind = _ind("x", CanidTaxon.DOMESTIC_DOG, "village")
    assert ind.sex is Sex.UNKNOWN
