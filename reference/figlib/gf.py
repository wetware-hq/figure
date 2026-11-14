"""Sample alphabet GF(2^8), AES polynomial x^8 + x^4 + x^3 + x + 1.

Storage is the raw byte. Addition is XOR. Integer addition of the same
bytes is a different operation and is not the field sum. Multiplication
is defined so a figure can store field elements, not so a viewer must
multiply pixels.
"""

from __future__ import annotations

from .errors import FigError

# 0x11B = 283. Bit 8 is the modulus; reduction XORs 0x1B into the low byte.
AES_POLY = 0x11B


def _byte(n: int) -> int:
    if not isinstance(n, int) or isinstance(n, bool) or n < 0 or n > 255:
        raise FigError(f"{n} is not a GF(2^8) byte")
    return n


def add(a: int, b: int) -> int:
    return _byte(a) ^ _byte(b)


def mul(a: int, b: int, poly: int = AES_POLY) -> int:
    a = _byte(a)
    b = _byte(b)
    if poly != AES_POLY:
        raise FigError("v1 multiplies only in the AES field")
    product = 0
    for _ in range(8):
        if b & 1:
            product ^= a
        a <<= 1
        if a & 0x100:
            a ^= poly
        a &= 0xFF
        b >>= 1
    return product
