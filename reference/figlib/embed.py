"""Embeddings that replace a DICOM pixel coordinate and a CIF atom site.

Both land on the same integer lattice. The Morton code does not contain the
unit. The quantum and the origin are properties of the lattice node.

DICOM Image Position (Patient) is the centre of voxel (column 0, row 0).
The first three image-orientation cosines are the direction of increasing
column. The next three are the direction of increasing row. Pixel Spacing is
(row spacing, column spacing). This module takes those quantities already
separated, as exact rationals of millimetres, and does not renormalise a
bad direction vector.

CIF fractional coordinates are wrapped into [0, 1) before the orthogonalization.
An orthogonal cell (all angles 90°) is exact. Any other cell is evaluated in
IEEE-754 binary64 and then quantized. v1 does not apply space-group symmetry.
"""

from __future__ import annotations

import math
from fractions import Fraction

from .errors import CoordinateCollision, FigError, OutOfLattice

Number = Fraction | int


def F(value: Number) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise FigError(f"{value!r} is not an exact rational")
    return Fraction(value)


def round_half_away(numer: int, denom: int) -> int:
    """Nearest integer to numer/denom. A tie moves away from zero."""
    if denom == 0:
        raise FigError("zero denominator")
    if denom < 0:
        numer = -numer
        denom = -denom
    sign = 1 if numer >= 0 else -1
    numer = abs(numer)
    quotient, remainder = divmod(numer, denom)
    if remainder * 2 >= denom:
        quotient += 1
    return sign * quotient


def round_half_away_float(value: float) -> int:
    if math.isnan(value) or math.isinf(value):
        raise FigError("non-finite coordinate")
    if value >= 0:
        return math.floor(value + 0.5)
    return -math.floor(-value + 0.5)


def quantize(value: Fraction, origin: Fraction, quantum: Fraction) -> int:
    if quantum <= 0:
        raise FigError("quantum must be positive")
    ratio = (value - origin) / quantum
    return round_half_away(ratio.numerator, ratio.denominator)


def _require(coord: int, bits: int) -> int:
    if coord < 0 or coord >= (1 << bits):
        raise OutOfLattice(f"{coord} is outside a {bits}-bit axis")
    return coord


def word_int(data: bytes) -> int:
    return int.from_bytes(data, "little")


def rescale(data: bytes, slope: Number, intercept: Number) -> Fraction:
    """DICOM rescale affine on the integer representative. Not field arithmetic."""
    return F(slope) * word_int(data) + F(intercept)


def dicom_voxel(
    column: int,
    row: int,
    slice_index: int,
    *,
    ipp_mm: tuple[Number, Number, Number],
    column_direction: tuple[Number, Number, Number],
    row_direction: tuple[Number, Number, Number],
    slice_direction: tuple[Number, Number, Number],
    column_spacing_mm: Number,
    row_spacing_mm: Number,
    slice_spacing_mm: Number,
    origin_m: tuple[Number, Number, Number],
    quantum_m: Number,
    bits: int,
) -> tuple[int, int, int]:
    ipp = tuple(F(v) for v in ipp_mm)
    dirs = (
        tuple(F(v) for v in column_direction),
        tuple(F(v) for v in row_direction),
        tuple(F(v) for v in slice_direction),
    )
    steps = (
        F(column) * F(column_spacing_mm),
        F(row) * F(row_spacing_mm),
        F(slice_index) * F(slice_spacing_mm),
    )
    patient_mm = [ipp[comp] for comp in range(3)]
    for axis in range(3):
        for comp in range(3):
            patient_mm[comp] += steps[axis] * dirs[axis][comp]
    quantum = F(quantum_m)
    origin = tuple(F(v) for v in origin_m)
    coords = []
    for comp in range(3):
        metres = patient_mm[comp] / 1000
        coords.append(_require(quantize(metres, origin[comp], quantum), bits))
    return tuple(coords)


def require_injective(coords: list[tuple[int, ...]]) -> None:
    seen: set[tuple[int, ...]] = set()
    for coord in coords:
        if coord in seen:
            raise CoordinateCollision(repr(coord))
        seen.add(coord)


def wrap_frac(value: Number) -> Fraction:
    return F(value) % 1


def cif_orthogonal(
    frac: tuple[Number, Number, Number],
    length_angstrom: tuple[Number, Number, Number],
    *,
    origin_m: tuple[Number, Number, Number] = (0, 0, 0),
    quantum_m: Number,
    bits: int,
) -> tuple[int, int, int]:
    angstrom = Fraction(1, 10**10)
    quantum = F(quantum_m)
    origin = tuple(F(v) for v in origin_m)
    coords = []
    for axis in range(3):
        metres = wrap_frac(frac[axis]) * F(length_angstrom[axis]) * angstrom
        coords.append(_require(quantize(metres, origin[axis], quantum), bits))
    return tuple(coords)


def cif_cell(
    frac: tuple[Number, Number, Number],
    length_angstrom: tuple[Number, Number, Number],
    angles_deg: tuple[Number, Number, Number],
    *,
    origin_m: tuple[Number, Number, Number] = (0, 0, 0),
    quantum_m: Number,
    bits: int,
) -> tuple[int, int, int]:
    """PDB orthogonalization, then quantization onto the lattice.

    Angles are alpha, beta, gamma in degrees. All-90° cells never enter binary64.
    """
    alpha, beta, gamma = (F(v) for v in angles_deg)
    if (alpha, beta, gamma) == (Fraction(90), Fraction(90), Fraction(90)):
        return cif_orthogonal(
            frac,
            length_angstrom,
            origin_m=origin_m,
            quantum_m=quantum_m,
            bits=bits,
        )
    ca = math.cos(math.radians(float(alpha)))
    cb = math.cos(math.radians(float(beta)))
    cg = math.cos(math.radians(float(gamma)))
    sg = math.sin(math.radians(float(gamma)))
    if sg == 0.0:
        raise FigError("degenerate cell")
    volume_term = 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg
    if volume_term <= 0.0:
        raise FigError("degenerate cell")
    height = math.sqrt(volume_term)
    a, b, c = (float(F(v)) for v in length_angstrom)
    xf, yf, zf = (float(wrap_frac(v)) for v in frac)
    cartesian_a = (
        a * xf + b * cg * yf + c * cb * zf,
        b * sg * yf + c * (ca - cb * cg) / sg * zf,
        c * height / sg * zf,
    )
    quantum = float(F(quantum_m))
    if quantum <= 0.0:
        raise FigError("quantum must be positive")
    origin = tuple(float(F(v)) for v in origin_m)
    coords = []
    for component, origin_component in zip(cartesian_a, origin):
        metres = component * 1e-10
        raw = (metres - origin_component) / quantum
        coords.append(_require(round_half_away_float(raw), bits))
    return tuple(coords)
