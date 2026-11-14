"""Reader for the restricted ObjectGraph a `.fig` link document uses.

The field bytes are headerless. This document is the link: lattice, alphabet,
GIF frame fields, and the CID of each field. Fetching a URL is the caller's
job; the reader checks an inlined hex field against the CID written on the node.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from .errors import FigError
from .field import cid_raw, parse
from .zorder import code_width, decode

_ATTR = re.compile(r'([A-Za-z_][A-Za-z0-9_-]*)=("(?:[^"\\]|\\.)*"|[^\s\]]+)')
_NODE = re.compile(r"^::node\s+\[(.*)\]\s*$")
_EDGE = re.compile(r"^->\[:([a-z0-9-]+)\]\s+(.+)$")


@dataclass(frozen=True)
class LinkedField:
    node_id: str
    cid: str
    blob: bytes
    samples: tuple[tuple[tuple[int, ...], bytes], ...]


@dataclass(frozen=True)
class LinkedFrame:
    node_id: str
    seq: int
    delay: int
    disposal: int
    transparent: bytes | None
    user_input: bool
    field: LinkedField


@dataclass(frozen=True)
class Document:
    node_ids: tuple[str, ...]
    edges: tuple[tuple[str, str, str], ...]
    loop: int
    background: bytes
    time_num: int
    time_den: int
    bits: int
    dims: int
    axes: tuple[str, ...]
    quantum_num: int
    quantum_den: int
    alphabet: str
    value_bytes: int
    irreducible: int | None
    frames: tuple[LinkedFrame, ...]


def _attrs(body: str) -> dict[str, str]:
    found = _ATTR.findall(body)
    if not found:
        raise FigError("node has no attributes")
    attrs: dict[str, str] = {}
    for key, value in found:
        if key in attrs:
            raise FigError(f"duplicate attribute {key}")
        if value.startswith('"'):
            value = value[1:-1]
        attrs[key] = value
    return attrs


def parse_document(text: str) -> Document:
    nodes: dict[str, dict[str, str]] = {}
    hex_blobs: dict[str, str] = {}
    edges: list[tuple[str, str, str]] = []
    order: list[str] = []
    current: str | None = None
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("```"):
            lang = line[3:].strip()
            block: list[str] = []
            i += 1
            while i < len(lines) and not lines[i].startswith("```"):
                block.append(lines[i])
                i += 1
            if i >= len(lines):
                raise FigError("unclosed fence")
            if lang == "hex" and current is not None:
                hex_blobs[current] = "".join(block)
            i += 1
            continue
        node = _NODE.match(line)
        if node:
            attrs = _attrs(node.group(1))
            node_id = attrs.get("id")
            if not node_id or node_id in nodes:
                raise FigError("node id missing or repeated")
            nodes[node_id] = attrs
            order.append(node_id)
            current = node_id
            i += 1
            continue
        edge = _EDGE.match(line)
        if line.startswith("->"):
            if edge is None or current is None:
                raise FigError(f"malformed edge: {line}")
            for target in edge.group(2).split(","):
                name = target.strip()
                if not name:
                    raise FigError("empty edge target")
                edges.append((current, edge.group(1), name))
            i += 1
            continue
        if line.startswith("::end"):
            current = None
        i += 1

    figures = [n for n in nodes.values() if n.get("type") == "figure"]
    lattices = [n for n in nodes.values() if n.get("type") == "lattice"]
    alphabets = [n for n in nodes.values() if n.get("type") == "alphabet"]
    if len(figures) != 1 or len(lattices) != 1 or len(alphabets) != 1:
        raise FigError("a link document needs one figure, one lattice, and one alphabet")
    figure, lattice, alphabet = figures[0], lattices[0], alphabets[0]
    bits = int(lattice["bits"])
    dims = int(lattice["dims"])
    value_bytes = int(alphabet["value-bytes"])
    name = alphabet["name"]
    if name == "gf2-8":
        if value_bytes != 1 or int(alphabet["irreducible"]) != 0x11B:
            raise FigError("gf2-8 is one byte in the AES field")
    elif name == "raw":
        if "irreducible" in alphabet:
            raise FigError("raw has no polynomial")
    else:
        raise FigError(f"unknown alphabet {name}")
    width = code_width(dims, bits)
    background = bytes.fromhex(figure["background"])
    time_num = int(figure.get("time-num", "1"))
    time_den = int(figure.get("time-den", "100"))
    if time_num < 1 or time_den < 1:
        raise FigError("time quantum must be a positive rational")
    if int(figure["loop"]) < 0:
        raise FigError("loop must be a non-negative integer")
    if int(lattice["quantum-num"]) < 1 or int(lattice["quantum-den"]) < 1:
        raise FigError("quantum must be positive")
    axes = tuple(lattice["axes"].split(","))
    if len(axes) != dims or any(part == "" for part in axes):
        raise FigError("axes must match dims")
    if len(background) != value_bytes:
        raise FigError("background width does not match the alphabet")

    frame_nodes = [n for n in nodes.values() if n.get("type") == "frame"]
    frame_nodes.sort(key=lambda n: int(n["seq"]))
    linked: list[LinkedFrame] = []
    for frame in frame_nodes:
        targets = [dst for src, kind, dst in edges if src == frame["id"] and kind == "contains"]
        fields = [t for t in targets if nodes.get(t, {}).get("type") == "field"]
        if len(fields) != 1:
            raise FigError(f"frame {frame['id']} must contain one field")
        field = nodes[fields[0]]
        raw_hex = hex_blobs.get(field["id"])
        if raw_hex is None:
            raise FigError(f"field {field['id']} has no inlined bytes")
        try:
            blob = bytes.fromhex(raw_hex)
        except ValueError as exc:
            raise FigError(f"field {field['id']} hex is not whole bytes") from exc
        digest = cid_raw(blob)
        if digest != field["cid"]:
            raise FigError(f"field {field['id']} cid does not match its bytes")
        samples = []
        for code, value in parse(blob, width, value_bytes):
            samples.append((decode(code, dims, bits), value))
        if int(frame["disposal"]) not in (0, 1, 2, 3):
            raise FigError("disposal is a GIF method 0..3")
        if int(frame["delay"]) < 0:
            raise FigError("delay must be a non-negative integer")
        flag = frame.get("user-input", "0")
        if flag not in ("0", "1"):
            raise FigError("user-input is 0 or 1")
        transparent = None
        if "transparent" in frame:
            transparent = bytes.fromhex(frame["transparent"])
        linked.append(
            LinkedFrame(
                node_id=frame["id"],
                seq=int(frame["seq"]),
                delay=int(frame["delay"]),
                disposal=int(frame["disposal"]),
                transparent=transparent,
                user_input=frame.get("user-input", "0") == "1",
                field=LinkedField(field["id"], digest, blob, tuple(samples)),
            )
        )
    if [frame.seq for frame in linked] != list(range(len(linked))):
        raise FigError("frame seq is not 0..n-1")
    seq_of = {frame.node_id: frame.seq for frame in linked}
    for src, kind, dst in edges:
        if kind != "precedes":
            continue
        if src not in seq_of or dst not in seq_of or seq_of[dst] != seq_of[src] + 1:
            raise FigError("precedes does not follow seq")
    return Document(
        node_ids=tuple(order),
        edges=tuple(edges),
        loop=int(figure["loop"]),
        background=background,
        time_num=time_num,
        time_den=time_den,
        bits=bits,
        dims=dims,
        axes=axes,
        quantum_num=int(lattice["quantum-num"]),
        quantum_den=int(lattice["quantum-den"]),
        alphabet=alphabet["name"],
        value_bytes=value_bytes,
        irreducible=int(alphabet["irreducible"]) if name == "gf2-8" else None,
        frames=tuple(linked),
    )
