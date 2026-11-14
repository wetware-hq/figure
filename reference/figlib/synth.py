"""SYNTH-FIG-0001, the synthetic lattice named by the specification.

Invented samples on a 4 by 4 lattice. Not a scan, not a patient, not a cell
refined from a crystal. The bytes are produced by the same functions a reader
checks, then written out so the example file cannot drift from the library.
"""

from __future__ import annotations

from .field import canonical_bytes, cid_raw
from .gf import add, mul
from .zorder import code_width, encode

BITS = 2
DIMS = 2


def _word(n: int) -> bytes:
    return bytes([n])


def frame_samples() -> tuple[tuple[tuple[tuple[int, int], bytes], ...], ...]:
    a = 0x57
    b = 0x83
    first = (
        ((0, 0), _word(a)),
        ((1, 0), _word(b)),
        ((0, 1), _word(add(a, b))),
        ((1, 1), _word(mul(a, b))),
    )
    second = (
        ((1, 1), _word(0x00)),
        ((2, 2), _word(0x11)),
    )
    return first, second


def field_blob(samples: tuple[tuple[tuple[int, int], bytes], ...]) -> bytes:
    width = code_width(DIMS, BITS)
    pairs = [(encode(coord, BITS), value) for coord, value in samples]
    return canonical_bytes(pairs, width)


def render() -> str:
    first, second = frame_samples()
    blobs = (field_blob(first), field_blob(second))
    cids = tuple(cid_raw(blob) for blob in blobs)
    hexes = tuple(blob.hex() for blob in blobs)
    return f"""::meta
title: SYNTH-FIG-0001
version: 1.0.0
media-type: text/vnd.fig+objectgraph
field-type: application/vnd.fig-field
::end

::schema
node-types: [figure, lattice, alphabet, frame, field, embedding, meta]
edge-types: [contains, precedes, supersedes]
::end

::index
# id | type | keywords
fig | figure | root, loop, background
lat | lattice | bits, quantum, axes
alph | alphabet | gf2-8, raw byte
fr0 | frame | delay 10, disposal 1
fld0 | field | headerless, cid
fr1 | frame | delay 5, disposal 2, transparent
fld1 | field | headerless, cid
emb-dicom | embedding | pixel centre, orientation
emb-cif | embedding | fractional site, cell
__view | meta | frame order
::end

::node [id=fig type=figure loop=1 background="00"]
::dense synthetic 4x4 lattice | two GIF frames | not a patient
->[:contains] lat, alph, fr0, fr1, emb-dicom, emb-cif
->[:supersedes] dicom-pixel-coordinate, cif-atom-site-coordinate
::end

::node [id=lat type=lattice bits=2 dims=2 axes=x,y quantum-num=1 quantum-den=10000000000]
::dense each cell is 1 angstrom | axis 0 is the low bit of each Morton group
::end

::node [id=alph type=alphabet name=gf2-8 value-bytes=1 irreducible=283]
::dense GF(2^8) | AES polynomial 283 | addition is XOR | the byte is the element
::end

::node [id=fr0 type=frame seq=0 delay=10 disposal=1 user-input=0]
::dense hold 10 centiseconds | leave the graphic in place
->[:contains] fld0
->[:precedes] fr1
::end

::node [id=fld0 type=field cid="{cids[0]}"]
::dense four records | code byte then sample byte | no header
```hex
{hexes[0]}
```
::end

::node [id=fr1 type=frame seq=1 delay=5 disposal=2 user-input=0 transparent="00"]
::dense hold 5 centiseconds | transparent hole on the product cell | then restore the graphic to background
->[:contains] fld1
::end

::node [id=fld1 type=field cid="{cids[1]}"]
::dense two records | 0x00 is stored and transparent | absence would be no record
```hex
{hexes[1]}
```
::end

::node [id=emb-dicom type=embedding kind=dicom]
::dense Image Position is the centre of column 0 row 0 | first orientation row increases the column | rescale is an integer affine and is not field addition
::end

::node [id=emb-cif type=embedding kind=cif]
::dense wrap the fractional site into [0,1) | orthogonal cells stay rational | other cells use binary64 then round half away from zero | no symmetry expansion
::end

::node [id=__view type=meta]
::dense fr0 draws four field elements | fr1 punches a transparent hole and writes one new cell | disposal 2 then clears both of fr1's stored cells
```mermaid
flowchart LR
  fr0["fr0 disposal 1"] --> fr1["fr1 disposal 2"]
  fr1 --> rest["product cell and (2,2) return to background"]
```
::end
"""
