"""NHGRI 722-genome sample-name taxonomy (breed / wild / village / mixed)."""

from __future__ import annotations

# Groups whose names contain these substrings are treated as wild/feral canids
# rather than domestic dogs (outgroups + coyote allele frequencies).
WILD_KEYWORDS = (
    "coyote",
    "wolf",
    "jackal",
    "fox",
    "dhole",
    "lycaon",
    "wilddog",
    "dingo",
    "ngsd",
    "newguineasinging",
)
# Free-roaming / unregistered dog populations. PER and BAN are large
# site-coded groups inferred (from sampling publications) to be Peruvian and
# Bangladeshi free-roaming dogs; the inference is flagged in outputs.
VILLAGE_PREFIXES = ("VillDog_",)
VILLAGE_KEYWORDS = ("IndigenousDog",)
INFERRED_VILLAGE_CODES = {"PER": "Peru", "BAN": "Bangladesh"}
# CFA. (24 samples) has no documented meaning; excluded from named candidates.
UNRESOLVED_CODES = {"CFA."}
# Domestic breeds whose names contain wild-canid substrings.
DOMESTIC_EXCEPTIONS = frozenset({"IrishWolfhound"})


def breed_of(sample: str) -> str:
    """Strip trailing digits from a panel sample id to recover the group name."""
    return sample.rstrip("0123456789")


def classify_group(group: str) -> str:
    """Categorize a parsed name-group: wild, village, mixed, ambiguous, or breed."""
    lowered = group.lower()
    if group in DOMESTIC_EXCEPTIONS:
        return "breed"
    if any(keyword in lowered for keyword in WILD_KEYWORDS):
        return "wild"
    if any(group.startswith(prefix) for prefix in VILLAGE_PREFIXES) or any(
        keyword.lower() in lowered for keyword in VILLAGE_KEYWORDS
    ):
        return "village"
    if group in INFERRED_VILLAGE_CODES or group in UNRESOLVED_CODES:
        return "village" if group in INFERRED_VILLAGE_CODES else "ambiguous"
    if "mix" in lowered or "unknown" in lowered:
        return "mixed"
    # Short initialisms (BC, GR, Helsinki_BC, run IDs, sample codes with digits)
    # cannot be mapped to breeds safely; keep them only in the pooled panel.
    if any(character.isdigit() for character in group) or len(group) <= 3:
        return "ambiguous"
    if group == "Helsinki_BC":
        return "ambiguous"
    return "breed"
