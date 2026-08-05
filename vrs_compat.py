#!/usr/bin/env python3
"""Minimal GA4GH VRS v2 computed identifiers used by IsoTag's XV tag.

This module intentionally handles only literal alleles that can be normalized
without fetching reference context. SNVs and substitutions are safe. Pure
insertions/deletions require reference-aware VRS normalization and are rejected.

Replaces the pre-VRS-v2 module (IsoTag v2.3 and earlier), which embedded the
full location object inline in the Allele digest and used the VRS v1
SequenceState type — both produce non-compliant identifiers that will not
match ClinVar/ClinGen VRS-indexed data. See VRS2_MIGRATION.md for details.
"""

import base64
import hashlib
import json
import re
from dataclasses import dataclass
from typing import Any, Dict, Tuple


REFGET_RE = re.compile(r"^(?:ga4gh:)?SQ\.[0-9A-Za-z_-]{32}$")


def sha512t24u(blob: bytes) -> str:
    """Return the GA4GH 24-byte truncated SHA-512 base64url digest."""
    return base64.urlsafe_b64encode(hashlib.sha512(blob).digest()[:24]).decode(
        "ascii"
    )


def canonical_json(obj: Dict[str, Any]) -> bytes:
    """Return canonical bytes for the ASCII-only VRS objects emitted here."""
    return json.dumps(
        obj, ensure_ascii=True, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def ga4gh_digest(obj: Dict[str, Any]) -> str:
    return sha512t24u(canonical_json(obj))


def normalize_substitution(start: int, ref: str, alt: str) -> Tuple[int, str, str]:
    """Trim shared literal flanks.

    If trimming creates an insertion or deletion, reference context is required
    for fully justified VRS normalization and this function raises ValueError.
    """
    ref = ref.upper()
    alt = alt.upper()
    while ref and alt and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        start += 1
    if not ref or not alt:
        raise ValueError(
            "Reference-aware VRS normalization is required for insertions/deletions"
        )
    return start, ref, alt


@dataclass(frozen=True)
class VRSAlleleV2:
    """A normalized literal substitution represented as a VRS v2 Allele."""

    refget_accession: str
    start: int
    ref: str
    alt: str

    def __post_init__(self) -> None:
        if not REFGET_RE.match(self.refget_accession):
            raise ValueError("A complete SQ.<sha512t24u> RefGet accession is required")
        if self.start < 0:
            raise ValueError("VRS coordinates must be non-negative")
        normalized = normalize_substitution(self.start, self.ref, self.alt)
        object.__setattr__(self, "start", normalized[0])
        object.__setattr__(self, "ref", normalized[1])
        object.__setattr__(self, "alt", normalized[2])

    @property
    def bare_refget_accession(self) -> str:
        prefix = "ga4gh:"
        if self.refget_accession.startswith(prefix):
            return self.refget_accession[len(prefix) :]
        return self.refget_accession

    def location_object(self) -> Dict[str, Any]:
        return {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference",
                "refgetAccession": self.bare_refget_accession,
            },
            "start": self.start,
            "end": self.start + len(self.ref),
        }

    def location_digest(self) -> str:
        return ga4gh_digest(self.location_object())

    def allele_digest_object(self) -> Dict[str, Any]:
        return {
            "type": "Allele",
            "location": self.location_digest(),
            "state": {
                "type": "LiteralSequenceExpression",
                "sequence": self.alt,
            },
        }

    def identifier(self) -> str:
        return "ga4gh:VA." + ga4gh_digest(self.allele_digest_object())


def test_vrs_compatibility() -> None:
    """Check official VRS v2 examples and the ClinVar BRCA1 fixture."""
    assert sha512t24u(b"ACGT") == "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"

    official_example = VRSAlleleV2(
        "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl", 44908821, "C", "T"
    )
    assert (
        official_example.location_digest()
        == "wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz"
    )
    assert (
        official_example.identifier()
        == "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt"
    )

    brca1 = VRSAlleleV2(
        "SQ.dLZ15tNO1Ur0IcGjwc3Sdi_0A6Yf4zm7", 43106486, "A", "G"
    )
    assert brca1.identifier() == "ga4gh:VA.gbvJw0s4OeAvloCeAM6BNNvrjFC_Dhc8"


if __name__ == "__main__":
    test_vrs_compatibility()
    print("VRS v2 compatibility checks passed")
