"""Falsifiers for the `.fig` coordinate rules.

Each test names the rule it tries to break. A passing run is the checker
failing to break that rule.
"""

from __future__ import annotations

import hashlib
import unittest
from fractions import Fraction
from pathlib import Path

from figlib.compose import Frame, play
from figlib.embed import (
    cif_cell,
    cif_orthogonal,
    dicom_voxel,
    quantize,
    require_injective,
    rescale,
    round_half_away,
    wrap_frac,
)
from figlib.errors import (
    CoordinateCollision,
    DuplicateCode,
    FigError,
    NonCanonical,
    OutOfLattice,
    UserInputRequired,
)
from figlib.field import (
    canonical_bytes,
    cid_raw,
    lookup,
    merkle_root,
    parse,
    require_cid,
    split_chunks,
)
from figlib.gf import add, mul
from figlib.link import parse_document
from figlib.synth import render
from figlib.zorder import coarsen, code_width, decode, encode

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "synth-0001.fig"

# Hand-interleaved for bits=2: bit 0 is x, bit 1 is y, bit 2 is x, bit 3 is y.
MORTON_4X4 = {
    (0, 0): 0,
    (1, 0): 1,
    (0, 1): 2,
    (1, 1): 3,
    (2, 0): 4,
    (3, 0): 5,
    (2, 1): 6,
    (3, 1): 7,
    (0, 2): 8,
    (1, 2): 9,
    (0, 3): 10,
    (1, 3): 11,
    (2, 2): 12,
    (3, 2): 13,
    (2, 3): 14,
    (3, 3): 15,
}

EMPTY_CID = "bafkreihdwdcefgh4dqkjv67uzcmw7ojee6xedzdetojuzjevtenxquvyku"
HELLO_CID = "bafkreibm6jg3ux5qumhcn2b3flc3tyu6dmlb4xa7u5bf44yegnrjhc4yeq"


def _poly_mul(a: int, b: int) -> int:
    """Independent GF(2) polynomial multiply, reduced by 0x11B."""
    product = 0
    for bit in range(8):
        if (b >> bit) & 1:
            product ^= a << bit
    modulus = 0x11B
    for bit in range(14, 7, -1):
        if (product >> bit) & 1:
            product ^= modulus << (bit - 8)
    return product & 0xFF


class ZOrderTests(unittest.TestCase):
    def test_hand_table_and_roundtrip(self):
        for coord, code in MORTON_4X4.items():
            self.assertEqual(encode(coord, 2), code)
            self.assertEqual(decode(code, 2, 2), coord)

    def test_three_axis_low_bits(self):
        self.assertEqual(encode((1, 0, 0), 1), 0b001)
        self.assertEqual(encode((0, 1, 0), 1), 0b010)
        self.assertEqual(encode((0, 0, 1), 1), 0b100)
        self.assertEqual(encode((1, 1, 1), 1), 0b111)
        self.assertEqual(decode(0b111, 3, 1), (1, 1, 1))

    def test_exhaustive_roundtrip(self):
        for bits, dims in ((6, 2), (4, 3), (3, 4)):
            limit = 1 << bits
            for index in range(limit**dims):
                coord = tuple((index // (limit**axis)) % limit for axis in range(dims))
                code = encode(coord, bits)
                self.assertEqual(decode(code, dims, bits), coord)
                self.assertLess(code, 1 << (bits * dims))

    def test_aligned_block_is_contiguous(self):
        bits, dims, levels = 5, 3, 2
        step = 1 << levels
        span = 1 << (levels * dims)
        for origin_index in (0, 1, 3):
            origin = tuple(origin_index * step for _ in range(dims))
            start = encode(origin, bits)
            self.assertEqual(start % span, 0)
            found = []
            for offset in range(span):
                local = []
                rest = offset
                for _ in range(dims):
                    local.append(rest % step)
                    rest //= step
                coord = tuple(origin[axis] + local[axis] for axis in range(dims))
                found.append(encode(coord, bits))
            self.assertEqual(sorted(found), list(range(start, start + span)))

    def test_coarsen_is_the_prefix(self):
        code = encode((3, 2), 4)
        self.assertEqual(coarsen(code, 2, 4, 0), code)
        self.assertEqual(coarsen(code, 2, 4, 1), encode((1, 1), 3))
        self.assertEqual(coarsen(code, 2, 4, 1), code >> 2)
        bond = encode((154, 0, 0), 48)
        self.assertEqual(coarsen(bond, 3, 48, 10), encode((154 >> 10, 0, 0), 38))
        self.assertEqual(coarsen(bond, 3, 48, 48), 0)

    def test_fine_cell_sits_inside_the_coarse_cell(self):
        quantum = 1
        levels = 10
        fine = 154
        coarse = fine >> levels
        fine_lo, fine_hi = fine * quantum, (fine + 1) * quantum
        step = quantum << levels
        coarse_lo, coarse_hi = coarse * step, (coarse + 1) * step
        self.assertLessEqual(coarse_lo, fine_lo)
        self.assertLessEqual(fine_hi, coarse_hi)

    def test_one_lattice_holds_a_bond_and_a_body(self):
        bits = 48
        bond = quantize(Fraction(154, 10**12), Fraction(0), Fraction(1, 10**12))
        body = quantize(Fraction(18, 10), Fraction(0), Fraction(1, 10**12))
        self.assertEqual(bond, 154)
        self.assertEqual(body, 1_800_000_000_000)
        self.assertGreaterEqual(body, 1 << 40)
        self.assertLess(body, 1 << 41)
        self.assertEqual(decode(encode((bond, 0, 0), bits), 3, bits), (bond, 0, 0))
        self.assertEqual(decode(encode((body, 0, 0), bits), 3, bits), (body, 0, 0))

    def test_rejects_out_of_range_and_overwide_codes(self):
        with self.assertRaises(OutOfLattice):
            encode((-1, 0), 4)
        with self.assertRaises(OutOfLattice):
            encode((16, 0), 4)
        with self.assertRaises(FigError):
            encode((0, 0, 0), 86)
        with self.assertRaises(OutOfLattice):
            decode(1 << 4, 2, 2)


class FieldTests(unittest.TestCase):
    def test_known_cids(self):
        self.assertEqual(cid_raw(b""), EMPTY_CID)
        self.assertEqual(cid_raw(b"hello"), HELLO_CID)
        require_cid(b"hello", HELLO_CID)
        with self.assertRaises(FigError):
            require_cid(b"hello", EMPTY_CID)

    def test_little_endian_sort_is_not_byte_sort(self):
        blob = canonical_bytes([(1, b"\xaa"), (256, b"\xbb")], 2)
        self.assertEqual(blob, bytes.fromhex("0100aa0001bb"))
        self.assertNotEqual(blob, bytes.fromhex("0001bb0100aa"))
        self.assertEqual(parse(blob, 2, 1), [(1, b"\xaa"), (256, b"\xbb")])
        with self.assertRaises(NonCanonical):
            parse(bytes.fromhex("0001bb0100aa"), 2, 1)

    def test_order_does_not_change_the_cid_and_a_byte_does(self):
        forward = canonical_bytes([(1, b"\x10"), (2, b"\x20")], 1)
        reverse = canonical_bytes([(2, b"\x20"), (1, b"\x10")], 1)
        self.assertEqual(forward, reverse)
        self.assertEqual(cid_raw(forward), cid_raw(reverse))
        changed = canonical_bytes([(1, b"\x11"), (2, b"\x20")], 1)
        self.assertNotEqual(cid_raw(forward), cid_raw(changed))

    def test_empty_duplicate_width_and_tail(self):
        self.assertEqual(canonical_bytes([], 1), b"")
        self.assertEqual(parse(b"", 1, 1), [])
        with self.assertRaises(DuplicateCode):
            canonical_bytes([(1, b"\x01"), (1, b"\x02")], 1)
        with self.assertRaises(FigError):
            canonical_bytes([(1, b"\x01"), (2, b"\x02\x03")], 1)
        with self.assertRaises(NonCanonical):
            parse(b"\x01\x02\x03", 1, 1)

    def test_header_is_absent_and_sparse_code_is_not_an_offset(self):
        samples = [(3, b"\x00"), (12, b"\x77")]
        blob = canonical_bytes(samples, 1)
        self.assertEqual(len(blob), 4)
        self.assertFalse(blob.startswith(b"GIF"))
        self.assertFalse(blob.startswith(b"DICM"))
        self.assertEqual(blob[0], 3)
        record = 2
        self.assertNotEqual(12 * record, blob.index(b"\x77"))
        self.assertEqual(lookup(blob, 12, 1, 1), b"\x77")
        self.assertIsNone(lookup(blob, 4, 1, 1))
        self.assertIsNone(lookup(b"", 0, 1, 1))

    def test_merkle_root_is_not_the_cid(self):
        blob = canonical_bytes([(0, b"a"), (1, b"b"), (2, b"c")], 1)
        chunks = split_chunks(blob, 1, 2)
        self.assertEqual(chunks, [blob[0:2], blob[2:4], blob[4:6]])

        def leaf(chunk: bytes) -> bytes:
            return hashlib.sha256(b"\x00" + chunk).digest()

        pair = hashlib.sha256(b"\x01" + leaf(chunks[0]) + leaf(chunks[1])).digest()
        expected = hashlib.sha256(b"\x01" + pair + leaf(chunks[2])).digest()
        self.assertEqual(merkle_root(chunks), expected)
        other = merkle_root(split_chunks(blob, 2, 2))
        self.assertNotEqual(other, expected)
        self.assertEqual(cid_raw(blob), cid_raw(b"".join(split_chunks(blob, 2, 2))))


class FiniteFieldTests(unittest.TestCase):
    def test_aes_vector_is_not_integer_arithmetic(self):
        self.assertEqual(mul(0x57, 0x83), 0xC1)
        self.assertEqual(mul(0x57, 0x83), _poly_mul(0x57, 0x83))
        self.assertEqual(add(0x57, 0x83), 0xD4)
        self.assertNotEqual(add(0x57, 0x83), (0x57 + 0x83) & 0xFF)
        self.assertEqual(mul(0x57, 1), 0x57)
        self.assertEqual(mul(0x57, 0), 0)
        self.assertEqual(mul(add(0x57, 0x83), 0x13), add(mul(0x57, 0x13), mul(0x83, 0x13)))

    def test_rescale_is_the_integer_representative(self):
        self.assertEqual(rescale(b"\x00", 1, -1024), -1024)
        self.assertEqual(rescale(b"\x01", 1, -1024), -1023)
        self.assertEqual(rescale(bytes([0x90, 0x01]), 1, -1024), 400 - 1024)


class TimeTests(unittest.TestCase):
    def test_disposal_two_punches_a_transparent_hole_and_three_does_not(self):
        first = Frame((((1, 1), b"\xc1"),), delay=10, disposal=1)
        hole = (((1, 1), b"\x00"), ((2, 2), b"\x11"))
        wiped = play(
            [first, Frame(hole, delay=5, disposal=2, transparent=b"\x00")],
            bits=4,
            background=b"\x00",
            loop=1,
        )
        kept = play(
            [first, Frame(hole, delay=5, disposal=3, transparent=b"\x00")],
            bits=4,
            background=b"\x00",
            loop=1,
        )
        self.assertEqual(wiped[1].visible[(1, 1)], b"\xc1")
        self.assertEqual(wiped[1].visible[(2, 2)], b"\x11")
        self.assertNotIn((1, 1), wiped[1].rested)
        self.assertNotIn((2, 2), wiped[1].rested)
        self.assertEqual(kept[1].rested[(1, 1)], b"\xc1")
        self.assertNotIn((2, 2), kept[1].rested)

    def test_region_disposal_restores_cells_the_frame_did_not_store(self):
        frames = [
            Frame((((0, 0), b"\x57"), ((1, 0), b"\x83")), disposal=1),
            Frame(
                (((0, 0), b"\x11"),),
                disposal=2,
                region=((0, 2), (0, 1)),
            ),
        ]
        moments = play(frames, bits=2, background=b"\x00", loop=1)
        self.assertEqual(moments[1].visible[(1, 0)], b"\x83")
        self.assertEqual(moments[1].rested, {})

    def test_disposal_zero_matches_leave_in_place(self):
        samples = (((0, 0), b"\x57"),)
        left = play([Frame(samples, disposal=0)], bits=2, background=b"\x00", loop=1)
        same = play([Frame(samples, disposal=1)], bits=2, background=b"\x00", loop=1)
        self.assertEqual(left[0].rested, same[0].rested)
        self.assertEqual(left[0].rested[(0, 0)], b"\x57")

    def test_draw_background_is_not_a_stored_zero_on_the_canvas(self):
        drawn = play(
            [Frame((((0, 0), b"\x00"),), disposal=1)],
            bits=2,
            background=b"\x00",
            loop=1,
        )
        self.assertEqual(drawn[0].visible, {})
        ink = play(
            [Frame((((0, 0), b"\x00"),), disposal=1)],
            bits=2,
            background=b"\x11",
            loop=1,
        )
        self.assertEqual(ink[0].visible[(0, 0)], b"\x00")

    def test_loop_count_and_user_input(self):
        frame = Frame((((0, 0), b"\x01"),), delay=7, disposal=1)
        twice = play([frame], bits=2, background=b"\x00", loop=2)
        self.assertEqual([moment.play_index for moment in twice], [0, 1])
        self.assertEqual(twice[0].delay, 7)
        with self.assertRaises(FigError):
            play([frame], bits=2, background=b"\x00", loop=0)
        capped = play([frame], bits=2, background=b"\x00", loop=0, plays=2)
        self.assertEqual(len(capped), 2)
        with self.assertRaises(UserInputRequired):
            play([Frame((((0, 0), b"\x01"),), user_input=True)], bits=2, background=b"\x00", loop=1)

    def test_sample_outside_region_is_rejected(self):
        frame = Frame((((2, 0), b"\x01"),), region=((0, 2), (0, 2)))
        with self.assertRaises(FigError):
            play([frame], bits=4, background=b"\x00", loop=1)


class EmbeddingTests(unittest.TestCase):
    def test_dicom_axes_and_half_away(self):
        common = dict(
            ipp_mm=(0, 0, 0),
            column_direction=(1, 0, 0),
            row_direction=(0, 1, 0),
            slice_direction=(0, 0, 1),
            column_spacing_mm=2,
            row_spacing_mm=5,
            slice_spacing_mm=4,
            origin_m=(0, 0, 0),
            quantum_m=Fraction(1, 1000),
            bits=16,
        )
        self.assertEqual(dicom_voxel(1, 0, 0, **common), (2, 0, 0))
        self.assertEqual(dicom_voxel(0, 1, 0, **common), (0, 5, 0))
        self.assertEqual(dicom_voxel(0, 0, 3, **common), (0, 0, 12))
        stretched = dict(common, column_direction=(2, 0, 0))
        self.assertEqual(dicom_voxel(1, 0, 0, **stretched), (4, 0, 0))
        half = dict(common, ipp_mm=(Fraction(3, 2), 0, 0))
        self.assertEqual(dicom_voxel(0, 0, 0, **half), (2, 0, 0))

    def test_negative_patient_space_needs_an_origin(self):
        args = dict(
            ipp_mm=(-10, 0, 0),
            column_direction=(1, 0, 0),
            row_direction=(0, 1, 0),
            slice_direction=(0, 0, 1),
            column_spacing_mm=1,
            row_spacing_mm=1,
            slice_spacing_mm=1,
            quantum_m=Fraction(1, 1000),
            bits=16,
        )
        with self.assertRaises(OutOfLattice):
            dicom_voxel(0, 0, 0, origin_m=(0, 0, 0), **args)
        self.assertEqual(
            dicom_voxel(0, 0, 0, origin_m=(Fraction(-1, 100), 0, 0), **args),
            (0, 0, 0),
        )

    def test_coarse_quantum_collides(self):
        args = dict(
            ipp_mm=(0, 0, 0),
            column_direction=(1, 0, 0),
            row_direction=(0, 1, 0),
            slice_direction=(0, 0, 1),
            column_spacing_mm=1,
            row_spacing_mm=1,
            slice_spacing_mm=1,
            origin_m=(0, 0, 0),
            quantum_m=Fraction(1, 100),
            bits=16,
        )
        landed = [dicom_voxel(column, 0, 0, **args) for column in (0, 1)]
        self.assertEqual(landed[0], landed[1])
        with self.assertRaises(CoordinateCollision):
            require_injective(landed)

    def test_cif_orthogonal_is_exact_and_wraps(self):
        quantum = Fraction(1, 10**10)
        self.assertEqual(
            cif_orthogonal((Fraction(1, 2), Fraction(1, 5), 0), (10, 10, 10), quantum_m=quantum, bits=8),
            (5, 2, 0),
        )
        self.assertEqual(
            cif_orthogonal((Fraction(1, 2), 0, 0), (5, 1, 1), quantum_m=quantum, bits=8),
            (3, 0, 0),
        )
        self.assertEqual(wrap_frac(Fraction(-1, 5)), Fraction(4, 5))
        self.assertEqual(wrap_frac(Fraction(6, 5)), Fraction(1, 5))
        self.assertEqual(
            cif_orthogonal((Fraction(-1, 5), 0, 0), (10, 1, 1), quantum_m=quantum, bits=8)[0],
            8,
        )

    def test_cif_gamma_60_vector(self):
        # Fractional 1 is the next lattice point, so it wraps to the origin.
        wrapped = cif_cell(
            (0, 1, 0),
            (2, 2, 2),
            (90, 90, 60),
            quantum_m=Fraction(1, 10**12),
            bits=16,
        )
        self.assertEqual(wrapped, (0, 0, 0))
        # Site halfway along b. Y = 50*sqrt(3) pm = sqrt(7500), nearer to 87 than to 86.
        coord = cif_cell(
            (0, Fraction(1, 2), 0),
            (2, 2, 2),
            (90, 90, 60),
            quantum_m=Fraction(1, 10**12),
            bits=16,
        )
        self.assertEqual(coord, (50, 87, 0))
        self.assertLess(87**2 - 7500, 7500 - 86**2)

    def test_round_half_away_from_zero(self):
        self.assertEqual(round_half_away(5, 2), 3)
        self.assertEqual(round_half_away(-3, 2), -2)
        self.assertEqual(round_half_away(1, 2), 1)
        self.assertEqual(round_half_away(1, 3), 0)
        self.assertEqual(round_half_away(2, 3), 1)


class AlphabetTests(unittest.TestCase):
    def test_raw_word_keeps_dicom_rescale(self):
        blob = canonical_bytes([(2, bytes([0x90, 0x01]))], 1)
        text = f"""::node [id=fig type=figure loop=1 background="0000"]
::end
::node [id=lat type=lattice bits=8 dims=1 axes=x quantum-num=1 quantum-den=1000]
::end
::node [id=alph type=alphabet name=raw value-bytes=2]
::end
::node [id=fr0 type=frame seq=0 delay=0 disposal=1 user-input=0]
->[:contains] fld0
::end
::node [id=fld0 type=field cid="{cid_raw(blob)}"]
```hex
{blob.hex()}
```
::end
"""
        document = parse_document(text)
        self.assertEqual(document.alphabet, "raw")
        self.assertEqual(document.frames[0].field.samples, (((2,), bytes([0x90, 0x01])),))
        self.assertEqual(rescale(document.frames[0].field.samples[0][1], 1, -1024), -624)

    def test_axes_and_precedes_must_agree(self):
        blob = canonical_bytes([(0, b"\x01")], 1)
        cid = cid_raw(blob)
        base = f"""::node [id=fig type=figure loop=1 background="00"]
::end
::node [id=lat type=lattice bits=2 dims=2 axes=x quantum-num=1 quantum-den=1]
::end
::node [id=alph type=alphabet name=gf2-8 value-bytes=1 irreducible=283]
::end
::node [id=fr0 type=frame seq=0 delay=0 disposal=1 user-input=0]
->[:contains] fld0
::end
::node [id=fld0 type=field cid="{cid}"]
```hex
{blob.hex()}
```
::end
"""
        with self.assertRaises(FigError):
            parse_document(base)
        swapped = base.replace("axes=x ", "axes=x,y ")
        swapped = swapped.replace(
            "->[:contains] fld0",
            "->[:contains] fld0\n->[:precedes] fr0",
        )
        with self.assertRaises(FigError):
            parse_document(swapped)

    def test_malformed_edge_is_rejected(self):
        text = '::node [id=fig type=figure loop=1 background="00"]\n->[:contains\n::end\n'
        with self.assertRaises(FigError):
            parse_document(text)

    def test_gf2_8_rejects_a_different_polynomial(self):
        blob = canonical_bytes([(0, b"\x01")], 1)
        text = f"""::node [id=fig type=figure loop=1 background="00"]
::end
::node [id=lat type=lattice bits=2 dims=2 axes=x,y quantum-num=1 quantum-den=1]
::end
::node [id=alph type=alphabet name=gf2-8 value-bytes=1 irreducible=1]
::end
::node [id=fr0 type=frame seq=0 delay=0 disposal=1 user-input=0]
->[:contains] fld0
::end
::node [id=fld0 type=field cid="{cid_raw(blob)}"]
```hex
{blob.hex()}
```
::end
"""
        with self.assertRaises(FigError):
            parse_document(text)


class SynthDocumentTests(unittest.TestCase):
    def test_example_matches_the_generator_and_plays(self):
        text = render()
        self.assertEqual(EXAMPLE.read_text(encoding="utf-8"), text)
        document = parse_document(text)
        self.assertEqual(document.bits, 2)
        self.assertEqual(document.dims, 2)
        self.assertEqual(document.axes, ("x", "y"))
        self.assertEqual(document.quantum_num, 1)
        self.assertEqual(document.quantum_den, 10**10)
        self.assertEqual(document.alphabet, "gf2-8")
        self.assertEqual(document.time_num, 1)
        self.assertEqual(document.time_den, 100)
        self.assertEqual(document.irreducible, 0x11B)
        self.assertEqual(document.value_bytes, 1)
        self.assertEqual(code_width(document.dims, document.bits), 1)
        self.assertEqual(len(document.frames), 2)
        self.assertEqual(document.frames[0].delay, 10)
        self.assertEqual(document.frames[0].disposal, 1)
        self.assertEqual(document.frames[1].delay, 5)
        self.assertEqual(document.frames[1].disposal, 2)
        self.assertEqual(document.frames[1].transparent, b"\x00")
        product = document.frames[0].field.samples[3]
        self.assertEqual(product, ((1, 1), bytes([mul(0x57, 0x83)])))
        field_sum = document.frames[0].field.samples[2][1]
        self.assertEqual(field_sum, bytes([add(0x57, 0x83)]))
        self.assertNotEqual(field_sum, bytes([(0x57 + 0x83) & 0xFF]))
        supersedes = [dst for _, kind, dst in document.edges if kind == "supersedes"]
        self.assertEqual(supersedes, ["dicom-pixel-coordinate", "cif-atom-site-coordinate"])
        moments = play(
            [
                Frame(frame.field.samples, delay=frame.delay, disposal=frame.disposal, transparent=frame.transparent)
                for frame in document.frames
            ],
            bits=document.bits,
            background=document.background,
            loop=document.loop,
        )
        self.assertEqual(
            moments[0].visible,
            {
                (0, 0): b"\x57",
                (1, 0): b"\x83",
                (0, 1): b"\xd4",
                (1, 1): b"\xc1",
            },
        )
        self.assertEqual(moments[1].visible[(2, 2)], b"\x11")
        self.assertEqual(moments[1].visible[(1, 1)], b"\xc1")
        self.assertNotIn((1, 1), moments[1].rested)
        self.assertNotIn((2, 2), moments[1].rested)
        self.assertEqual(moments[1].rested[(0, 1)], b"\xd4")
        for frame in document.frames:
            self.assertEqual(cid_raw(frame.field.blob), frame.field.cid)
            self.assertFalse(frame.field.blob.startswith(b"GIF89a"))

    def test_spec_pins_the_same_identifiers(self):
        text = (ROOT / "fig.md").read_text(encoding="utf-8")
        document = parse_document(EXAMPLE.read_text(encoding="utf-8"))
        figure_cid = cid_raw(EXAMPLE.read_bytes())
        for frame in document.frames:
            self.assertIn(frame.field.cid, text)
        self.assertIn(figure_cid, text)
        self.assertIn(EMPTY_CID, text)
        self.assertIn(HELLO_CID, text)
        self.assertIn(document.frames[0].field.blob.hex(), text)
        self.assertIn(document.frames[1].field.blob.hex(), text)


if __name__ == "__main__":
    unittest.main()
