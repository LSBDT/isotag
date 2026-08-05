#!/usr/bin/env python3
"""
IsoTag shared utilities — functions used by multiple IsoTag scripts.

Importing from here ensures a single definition for each shared function;
a bug fix in one place propagates to all scripts automatically.
"""

import hashlib
import base64


def mask_ambiguous_bases(sequence: str, keep_ambiguous: bool = False) -> str:
    """
    Mask ambiguous IUPAC codes with 'N' for consistent RefGet hashing.

    This ensures UCSC hg38 and NCBI GRCh38 produce identical chromosome hashes
    despite differences in ambiguous base representation (R, Y, S, W, K, M, etc.)

    Args:
        sequence: DNA sequence string
        keep_ambiguous: If True, keep R/Y/etc as-is (may cause cross-genome incompatibility)

    Returns:
        Normalized sequence for hashing
    """
    if keep_ambiguous:
        return sequence.upper()

    # Convert ambiguous IUPAC bases to N
    canonical_bases = set('ACGTN')
    normalized = ''.join(base if base in canonical_bases else 'N'
                         for base in sequence.upper())
    return normalized


def sha512t24u_fallback(blob: bytes) -> str:
    """
    Generate base64url-encoded, truncated SHA-512 digest (GA4GH sha512t24u algorithm).

    This is the fallback implementation used when vrs_compat is not installed.
    Prefer importing sha512t24u from vrs_compat when available.
    """
    digest_size = 24
    digest = hashlib.sha512(blob).digest()
    tdigest_b64us = base64.urlsafe_b64encode(digest[:digest_size])
    return tdigest_b64us.decode("ascii")
