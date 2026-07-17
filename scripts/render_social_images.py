#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = [
    ROOT / "public" / "dotmatch-og.png",
    ROOT / "public" / "dotmatch-twitter.png",
]

W, H = 1200, 630
INK = (16, 21, 19)
MUTED = (72, 88, 81)
LINE = (215, 228, 222)
BG = (247, 251, 249)
WHITE = (255, 255, 255)
GREEN = (15, 107, 87)
BLUE = (23, 79, 134)
AMBER = (153, 107, 19)
PURPLE = (102, 80, 157)
RED = (168, 61, 54)


def font(size: int, *, weight: str = "regular") -> ImageFont.FreeTypeFont:
    candidates = [
        Path("/System/Library/Fonts/Supplemental/Arial.ttf"),
        Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
        Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Regular.ttf"),
        Path("/System/Library/Fonts/SFNS.ttf"),
        Path("/System/Library/Fonts/HelveticaNeue.ttc"),
        Path("/Library/Fonts/Arial.ttf"),
    ]
    bold_candidates = [
        Path("/System/Library/Fonts/Supplemental/Arial Bold.ttf"),
        Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
        Path("/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf"),
        Path("/System/Library/Fonts/SFNS.ttf"),
        Path("/System/Library/Fonts/HelveticaNeue.ttc"),
        Path("/Library/Fonts/Arial Bold.ttf"),
    ]
    for path in bold_candidates if weight == "bold" else candidates:
        if path.exists():
            return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default(size=size)


def text_block(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    *,
    font_obj: ImageFont.FreeTypeFont,
    fill: tuple[int, int, int],
    width: int,
    line_gap: int,
) -> int:
    x, y = xy
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        candidate = f"{current} {word}".strip()
        if current and font_obj.getlength(candidate) > width:
            lines.append(current)
            current = word
        else:
            current = candidate
    if current:
        lines.append(current)

    bbox = font_obj.getbbox("Ag")
    line_height = bbox[3] - bbox[1]
    for line in lines:
        draw.text((x, y), line, font=font_obj, fill=fill)
        y += line_height + line_gap
    return y


def rounded_box(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int, int, int],
    *,
    fill: tuple[int, int, int],
    outline: tuple[int, int, int] | None = None,
    radius: int = 8,
    width: int = 2,
) -> None:
    draw.rounded_rectangle(xy, radius=radius, fill=fill, outline=outline, width=width if outline else 1)


def render_card() -> Image.Image:
    card = Image.new("RGB", (W, H), BG)
    draw = ImageDraw.Draw(card)

    # Subtle technical grid.
    for x in range(0, W, 32):
        draw.line((x, 0, x, H), fill=(232, 241, 236), width=1)
    for y in range(0, H, 32):
        draw.line((0, y, W, y), fill=(232, 241, 236), width=1)

    rounded_box(draw, (64, 58, 1136, 572), fill=WHITE, outline=LINE, radius=10)

    mark_x, mark_y = 96, 90
    rounded_box(draw, (mark_x, mark_y, mark_x + 62, mark_y + 62), fill=WHITE, outline=INK, radius=12, width=3)
    draw.line((mark_x + 31, mark_y + 8, mark_x + 31, mark_y + 54), fill=(15, 107, 87), width=5)
    draw.line((mark_x + 8, mark_y + 31, mark_x + 54, mark_y + 31), fill=(23, 79, 134), width=5)
    draw.text((178, 88), "DotMatch", font=font(38, weight="bold"), fill=INK)
    draw.text(
        (180, 132),
        "Known-target read assignment",
        font=font(23, weight="bold"),
        fill=GREEN,
    )

    heading_end = text_block(
        draw,
        (96, 210),
        "See which reads match—and which do not.",
        font_obj=font(62, weight="bold"),
        fill=INK,
        width=610,
        line_gap=7,
    )
    support_end = text_block(
        draw,
        (100, heading_end + 22),
        "Unique, ambiguous, unmatched, and invalid outcomes stay visible.",
        font_obj=font(26),
        fill=MUTED,
        width=610,
        line_gap=7,
    )
    if support_end > 540:
        raise ValueError(f"Social card copy exceeds the safe area: y={support_end}")

    panel = (760, 128, 1088, 520)
    rounded_box(draw, panel, fill=(250, 253, 251), outline=LINE, radius=8)
    draw.text((790, 158), "ONE OUTCOME PER READ", font=font(16, weight="bold"), fill=GREEN)
    draw.text((790, 187), "known target window", font=font(22, weight="bold"), fill=INK)

    outcomes = [
        ("unique", GREEN, (236, 248, 241)),
        ("ambiguous", AMBER, (255, 248, 230)),
        ("none", PURPLE, (244, 241, 251)),
        ("invalid", RED, (255, 242, 240)),
    ]
    y = 232
    for label, color, fill in outcomes:
        rounded_box(draw, (790, y, 1058, y + 58), fill=fill, outline=color, radius=8)
        draw.rectangle((813, y + 19, 831, y + 39), fill=color)
        draw.text((852, y + 14), label, font=font(22, weight="bold"), fill=INK)
        y += 69

    return card


def main() -> None:
    card = render_card()
    for out in OUTPUTS:
        card.save(out, "PNG", optimize=True)
        print(out)


if __name__ == "__main__":
    main()
