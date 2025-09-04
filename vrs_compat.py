#!/usr/bin/env python3
"""
VRS-Compatible Module for IsoTag v1.0

Contains official VRS code copied from GA4GH VRS-Python implementation
for strict VRS compliance in VARSID calculation.

Source: https://github.com/ga4gh/vrs-python
License: Apache License 2.0
"""

import base64
import hashlib
import json
from typing import Dict, Any

try:
    from canonicaljson import encode_canonical_json
except ImportError:
    # Fallback for basic JSON serialization
    def encode_canonical_json(obj):
        return json.dumps(obj, sort_keys=True, separators=(',', ':')).encode('utf-8')


def sha512t24u(blob: bytes) -> str:
    """Generate a base64url-encode, truncated SHA-512 digest for given
    binary data (OFFICIAL VRS IMPLEMENTATION)

    The sha512t24u digest is a convention for constructing and
    formatting digests for use as object identifiers. Specifically::

        * generate a SHA512 digest on binary data
        * truncate at 24 bytes
        * encode using base64url encoding

    Examples:
    >>> sha512t24u(b"")
    'z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc'

    >>> sha512t24u(b"ACGT")
    'aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2'

    Source: ga4gh.core.digests.py
    """
    digest_size = 24
    digest = hashlib.sha512(blob).digest()
    tdigest_b64us = base64.urlsafe_b64encode(digest[:digest_size])
    return tdigest_b64us.decode("ascii")


def vrs_serialize(obj: Dict[str, Any]) -> bytes:
    """Serialize an object for VRS digest computation using canonical JSON
    
    Args:
        obj: Dictionary representation of VRS object
        
    Returns:
        Canonical JSON bytes for hashing
        
    Source: Adapted from ga4gh.core.identifiers.py
    """
    return encode_canonical_json(obj)


def calculate_vrs_digest(vrs_object: Dict[str, Any]) -> str:
    """Calculate VRS digest for any VRS object
    
    Args:
        vrs_object: Dictionary representation of VRS object
        
    Returns:
        32-character VRS digest (without prefix)
        
    This is the core VRS digest calculation used across all VRS objects.
    """
    serialized = vrs_serialize(vrs_object)
    return sha512t24u(serialized)


class VRSVariant:
    """VRS-compliant variant representation for VARSID calculation
    
    This class creates VRS-compliant variant objects that can be used
    to generate official VRS digests for variants.
    """
    
    def __init__(self, refget_id: str, position: int, ref: str, alt: str):
        """Initialize VRS variant
        
        Args:
            refget_id: RefGet sequence identifier (e.g., SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul)
            position: 0-based position on sequence
            ref: Reference allele sequence
            alt: Alternate allele sequence
        """
        self.refget_id = refget_id
        self.position = position
        self.ref = ref
        self.alt = alt
    
    def to_vrs_object(self) -> Dict[str, Any]:
        """Convert to VRS Allele object format for digest calculation
        
        Returns:
            Dictionary representing VRS Allele object
        """
        # Create VRS SequenceLocation
        location = {
            "type": "SequenceLocation",
            "sequenceReference": {
                "type": "SequenceReference", 
                "refgetAccession": self.refget_id
            },
            "start": self.position,
            "end": self.position + len(self.ref)
        }
        
        # Create VRS Allele
        allele = {
            "type": "Allele",
            "location": location,
            "state": {
                "type": "SequenceState",
                "sequence": self.alt
            }
        }
        
        return allele
    
    def calculate_vrs_id(self) -> str:
        """Calculate official VRS ID for this variant
        
        Returns:
            Full VRS identifier (e.g., ga4gh:VA.abc123...)
        """
        vrs_object = self.to_vrs_object()
        digest = calculate_vrs_digest(vrs_object)
        return f"ga4gh:VA.{digest}"
    
    def calculate_digest_only(self) -> str:
        """Calculate just the digest portion (32 characters)
        
        Returns:
            32-character digest for use in IsoTag XV tags
        """
        vrs_object = self.to_vrs_object()
        return calculate_vrs_digest(vrs_object)


def test_vrs_compatibility():
    """Test VRS compatibility with known examples"""
    # Test sha512t24u with known values
    assert sha512t24u(b"ACGT") == "aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2"
    assert sha512t24u(b"") == "z4PhNX7vuL3xVChQ1m2AB9Yg5AULVxXc"
    print("✅ VRS sha512t24u compatibility verified")


if __name__ == "__main__":
    test_vrs_compatibility()
