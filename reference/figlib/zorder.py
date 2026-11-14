"""Morton Z-order on a finite lattice.

Each axis is an element of Z/2^b Z. Interleaving their bits is a permutation
of the bit vector, so it is a linear automorphism of the vector space
(GF(2))^(dims*bits). The code itself lives in the ring Z/2^(dims*bits) Z,
which is a field only when the width is 1. Callers that need a field
operation use the sample alphabet, not the address.
"""

from __future__ import annotations

from .errors import FigError, OutOfLattice

MAX_CODE_BITS = 256


def check_shape(dims: int, bits: int) -> None:
    if dims < 1 or bits < 1:
        raise FigError("dims and bits must be positive")
    if dims * bits > MAX_CODE_BITS:
        raise FigError(f"code wider than {MAX_CODE_BITS} bits")


def code_width(dims: int, bits: int) -> int:
    """Whole bytes used to store one code. There is no bit-packed header."""
    check_shape(dims, bits)
    return (dims * bits + 7) // 8


def encode(coords: tuple[int, ...], bits: int) -> int:
    """Interleave axis bits. Axis 0 occupies the least significant bit of each group."""
    dims = len(coords)
    check_shape(dims, bits)
    limit = 1 << bits
    for x in coords:
        if not isinstance(x, int) or isinstance(x, bool) or x < 0 or x >= limit:
            raise OutOfLattice(f"{x} is outside [0, {limit})")
    code = 0
    for k in range(bits):
        for i, x in enumerate(coords):
            code |= ((x >> k) & 1) << (k * dims + i)
    return code


def decode(code: int, dims: int, bits: int) -> tuple[int, ...]:
    check_shape(dims, bits)
    width = dims * bits
    if not isinstance(code, int) or isinstance(code, bool) or code < 0 or code >= (1 << width):
        raise OutOfLattice(f"code {code} does not fit in {width} bits")
    coords = [0] * dims
    for k in range(bits):
        for i in range(dims):
            if (code >> (k * dims + i)) & 1:
                coords[i] |= 1 << k
    return tuple(coords)


def coarsen(code: int, dims: int, bits: int, levels: int) -> int:
    """Drop `levels` finest bit-planes. The result is the code at quantum * 2^levels.

    A coarse code is a prefix of the fine code: the low `levels * dims` bits
    are discarded and the remaining bits stay in place. `levels == bits`
    collapses the lattice to the single code 0.
    """
    check_shape(dims, bits)
    if not isinstance(levels, int) or isinstance(levels, bool) or levels < 0 or levels > bits:
        raise FigError("coarsen levels must lie in 0..bits")
    width = dims * bits
    if code < 0 or code >= (1 << width):
        raise OutOfLattice(f"code {code} does not fit in {width} bits")
    return code >> (levels * dims)
