"""Failures a conforming reader must report, not repair."""


class FigError(ValueError):
    """The object breaks a normative rule of the `.fig` coordinate format."""


class OutOfLattice(FigError):
    """A coordinate or code does not lie in the declared finite lattice."""


class DuplicateCode(FigError):
    """Two samples claim the same Morton code."""


class NonCanonical(FigError):
    """Field bytes are not the sorted headerless record stream."""


class CoordinateCollision(FigError):
    """Two source positions quantize to one lattice cell."""


class UserInputRequired(FigError):
    """A frame sets the GIF user-input flag; an agent must not invent the input."""
