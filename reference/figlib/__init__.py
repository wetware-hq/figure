"""Reference checker for the `.fig` coordinate format."""

from .compose import Frame, Moment, play
from .embed import cif_cell, cif_orthogonal, dicom_voxel, quantize, rescale
from .errors import FigError
from .field import canonical_bytes, cid_raw, lookup, merkle_root, parse
from .gf import add, mul
from .link import parse_document
from .zorder import coarsen, code_width, decode, encode

__all__ = [
    "FigError",
    "Frame",
    "Moment",
    "add",
    "canonical_bytes",
    "cid_raw",
    "cif_cell",
    "cif_orthogonal",
    "coarsen",
    "code_width",
    "decode",
    "dicom_voxel",
    "encode",
    "lookup",
    "merkle_root",
    "mul",
    "parse",
    "parse_document",
    "play",
    "quantize",
    "rescale",
]
