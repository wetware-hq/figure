"""Headerless field bytes and their content address.

A field is a list of records sorted by Morton code ascending. Each record is
the little-endian code padded to whole bytes, then the raw sample word.
Nothing precedes the first record. The same bytes in another order are a
different file and are rejected by the strict reader; the writer sorts.

Lexicographic order of the little-endian bytes is not code order. Code 1
encodes as 01 00 and code 256 as 00 01, so a byte sort would reverse them.
"""

from __future__ import annotations

import hashlib
from typing import Iterable

from .errors import DuplicateCode, FigError, NonCanonical, OutOfLattice

_B32 = "abcdefghijklmnopqrstuvwxyz234567"


def _uvarint(n: int) -> bytes:
    if n < 0:
        raise FigError("negative varint")
    out = bytearray()
    while True:
        byte = n & 0x7F
        n >>= 7
        if n:
            out.append(byte | 0x80)
        else:
            out.append(byte)
            return bytes(out)


def _base32(data: bytes) -> str:
    acc = 0
    bits = 0
    chars: list[str] = []
    for byte in data:
        acc = (acc << 8) | byte
        bits += 8
        while bits >= 5:
            bits -= 5
            chars.append(_B32[(acc >> bits) & 31])
    if bits:
        chars.append(_B32[(acc << (5 - bits)) & 31])
    return "".join(chars)


def cid_raw(data: bytes) -> str:
    """CIDv1, multicodec raw (0x55), multihash sha2-256, multibase base32."""
    digest = hashlib.sha256(data).digest()
    binary = _uvarint(1) + _uvarint(0x55) + _uvarint(0x12) + _uvarint(len(digest)) + digest
    return "b" + _base32(binary)


def require_cid(body: bytes, declared: str) -> None:
    """A locator is HTTP. The bytes it returns must be the declared field, nothing else."""
    if cid_raw(body) != declared:
        raise FigError("locator bytes do not match the declared cid")


def canonical_bytes(samples: Iterable[tuple[int, bytes]], code_bytes: int) -> bytes:
    if code_bytes < 1 or code_bytes > 32:
        raise FigError("code_bytes must lie in 1..32")
    pairs = list(samples)
    seen: set[int] = set()
    limit = 1 << (8 * code_bytes)
    for code, value in pairs:
        if not isinstance(code, int) or isinstance(code, bool) or code < 0 or code >= limit:
            raise OutOfLattice(f"code {code} does not fit in {code_bytes} bytes")
        if not value:
            raise FigError("empty sample word")
        if code in seen:
            raise DuplicateCode(f"duplicate code {code}")
        seen.add(code)
    pairs.sort(key=lambda item: item[0])
    width = len(pairs[0][1]) if pairs else None
    out = bytearray()
    for code, value in pairs:
        if len(value) != width:
            raise FigError("sample words differ in width")
        out += code.to_bytes(code_bytes, "little")
        out += value
    return bytes(out)


def parse(blob: bytes, code_bytes: int, value_bytes: int) -> list[tuple[int, bytes]]:
    """Strict reader. Rejects a prefix, a tail, a duplicate, and unsorted input."""
    if code_bytes < 1 or value_bytes < 1:
        raise FigError("widths must be positive")
    record = code_bytes + value_bytes
    if len(blob) % record != 0:
        raise NonCanonical("field length is not a whole number of records")
    found: list[tuple[int, bytes]] = []
    previous = -1
    for offset in range(0, len(blob), record):
        code = int.from_bytes(blob[offset : offset + code_bytes], "little")
        if code <= previous:
            raise NonCanonical(f"code {code} is not strictly after {previous}")
        previous = code
        found.append((code, blob[offset + code_bytes : offset + record]))
    return found


def lookup(blob: bytes, code: int, code_bytes: int, value_bytes: int) -> bytes | None:
    """Binary search. A missing code has no byte offset; absence is not a zero word."""
    found = parse(blob, code_bytes, value_bytes)
    lo = 0
    hi = len(found)
    while lo < hi:
        mid = (lo + hi) // 2
        got = found[mid][0]
        if got == code:
            return found[mid][1]
        if got < code:
            lo = mid + 1
        else:
            hi = mid
    return None


def merkle_root(chunks: list[bytes]) -> bytes:
    """RFC 6962 tree. Leaf is SHA-256(0x00 || chunk); parent is SHA-256(0x01 || left || right).

    The field CID is the hash of the whole canonical blob. This root is an
    extra commitment for chunked transport. Chunking does not change the CID.
    """
    if not chunks:
        raise FigError("merkle tree needs a chunk")

    def node(items: list[bytes]) -> bytes:
        if len(items) == 1:
            return hashlib.sha256(b"\x00" + items[0]).digest()
        split = 1 << (len(items).bit_length() - 1)
        if split == len(items):
            split >>= 1
        return hashlib.sha256(b"\x01" + node(items[:split]) + node(items[split:])).digest()

    return node(chunks)


def split_chunks(blob: bytes, chunk_records: int, record_size: int) -> list[bytes]:
    if chunk_records < 1 or record_size < 1:
        raise FigError("chunk size must be positive")
    if len(blob) % record_size != 0:
        raise NonCanonical("field length is not a whole number of records")
    step = chunk_records * record_size
    if not blob:
        return [b""]
    return [blob[i : i + step] for i in range(0, len(blob), step)]
