"""Frame time borrowed from GIF89a, on a sparse finite field.

Visible order for each frame:

1. Draw every stored sample that is not the transparent word.
2. Hold for `delay` time quanta (the default quantum is 1/100 s, the GIF centisecond).
3. Dispose the graphic, then start the next frame.

The graphic area is the set of stored codes, or the declared region if one
is set. Disposal method 2 restores that whole area to the background,
including cells whose word was transparent. Disposal method 3 restores only
cells the frame overwrote, which is the GIF wording "area overwritten".
Method 0 is defined as method 1 so the next canvas is a function, not a
decoder choice. Method 0 in GIF89a leaves the action unspecified.

The canvas is a total function. A cell with no entry displays the background
word. Absence in the field is not the zero sample, and a transparent word is
not the background unless they happen to be the same byte string.
"""

from __future__ import annotations

from dataclasses import dataclass

from .errors import FigError, UserInputRequired
from .zorder import encode


@dataclass(frozen=True)
class Frame:
    samples: tuple[tuple[tuple[int, ...], bytes], ...]
    delay: int = 0
    disposal: int = 1
    transparent: bytes | None = None
    user_input: bool = False
    region: tuple[tuple[int, int], ...] | None = None


@dataclass(frozen=True)
class Moment:
    play_index: int
    seq: int
    delay: int
    visible: dict[tuple[int, ...], bytes]
    rested: dict[tuple[int, ...], bytes]


def _check_frame(frame: Frame, dims: int, bits: int, value_bytes: int) -> None:
    if frame.disposal not in (0, 1, 2, 3):
        raise FigError(f"disposal {frame.disposal} is not a GIF method")
    if not isinstance(frame.delay, int) or isinstance(frame.delay, bool) or frame.delay < 0:
        raise FigError("delay must be a non-negative integer of time quanta")
    if frame.transparent is not None and len(frame.transparent) != value_bytes:
        raise FigError("transparent word has the wrong width")
    seen: set[tuple[int, ...]] = set()
    for coord, value in frame.samples:
        if len(coord) != dims:
            raise FigError("sample dimension does not match the lattice")
        if coord in seen:
            raise FigError(f"duplicate coordinate {coord}")
        seen.add(coord)
        if len(value) != value_bytes:
            raise FigError("sample word has the wrong width")
        encode(coord, bits)
        if frame.region is not None and not _in_region(coord, frame.region):
            raise FigError(f"{coord} lies outside the frame region")
    if frame.region is not None:
        if len(frame.region) != dims:
            raise FigError("region dimension does not match the lattice")
        for lo, hi in frame.region:
            if lo < 0 or hi < lo or hi > (1 << bits):
                raise FigError("region is outside the lattice")


def _in_region(coord: tuple[int, ...], region: tuple[tuple[int, int], ...]) -> bool:
    return all(lo <= c < hi for c, (lo, hi) in zip(coord, region))


def _dims_of(frames: list[Frame]) -> int:
    for frame in frames:
        if frame.samples:
            return len(frame.samples[0][0])
        if frame.region is not None:
            return len(frame.region)
    raise FigError("lattice dimension is unknown")


def play(
    frames: tuple[Frame, ...] | list[Frame],
    *,
    bits: int,
    background: bytes,
    loop: int,
    plays: int | None = None,
) -> list[Moment]:
    """Return one moment per displayed frame.

    `loop` is the Netscape play count: a positive count is how many times the
    frame list is shown, and 0 means unbounded. An unbounded list needs an
    explicit `plays` cap so a checker cannot loop forever. `plays`, when set,
    is that cap and wins over `loop`.
    """
    if not frames:
        raise FigError("a figure has no frames")
    if not isinstance(loop, int) or isinstance(loop, bool) or loop < 0:
        raise FigError("loop must be a non-negative integer")
    if any(frame.user_input for frame in frames):
        raise UserInputRequired("GIF user-input flag is set")
    dims = _dims_of(list(frames))
    value_bytes = len(background)
    if value_bytes < 1:
        raise FigError("background word is empty")
    for frame in frames:
        _check_frame(frame, dims, bits, value_bytes)
    count = plays if plays is not None else loop
    if not isinstance(count, int) or isinstance(count, bool) or count < 1:
        raise FigError("unbounded loop needs an explicit positive play count")

    canvas: dict[tuple[int, ...], bytes] = {}
    moments: list[Moment] = []
    for play_index in range(count):
        for seq, frame in enumerate(frames):
            before = dict(canvas)
            drawn: list[tuple[int, ...]] = []
            for coord, value in frame.samples:
                if frame.transparent is not None and value == frame.transparent:
                    continue
                drawn.append(coord)
                if value == background:
                    canvas.pop(coord, None)
                else:
                    canvas[coord] = value
            visible = dict(canvas)
            disposal = 1 if frame.disposal == 0 else frame.disposal
            if disposal == 1:
                pass
            elif disposal == 2:
                if frame.region is None:
                    for coord, _ in frame.samples:
                        canvas.pop(coord, None)
                else:
                    for coord in list(canvas):
                        if _in_region(coord, frame.region):
                            canvas.pop(coord, None)
            else:
                for coord in drawn:
                    if coord in before:
                        canvas[coord] = before[coord]
                    else:
                        canvas.pop(coord, None)
            moments.append(
                Moment(play_index, seq, frame.delay, visible, dict(canvas))
            )
    return moments
