#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from textwrap import wrap

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "public" / "dotmatch-social-art.png"
OUTPUTS = [
    ROOT / "public" / "dotmatch-og.png",
    ROOT / "public" / "dotmatch-twitter.png",
]

W, H = 1200, 630
INK = (16, 21, 19)
MUTED = (73, 88, 81)
GREEN = (14, 124, 90)
GREEN_2 = (35, 176, 130)
BLUE = (29, 102, 209)
AMBER = (218, 151, 15)
LINE = (213, 226, 219)
WHITE = (250, 252, 250)


def font(size: int, *, weight: str = "regular") -> ImageFont.FreeTypeFont:
    candidates = [
        Path("/System/Library/Fonts/SFNS.ttf"),
        Path("/System/Library/Fonts/HelveticaNeue.ttc"),
        Path("/Library/Fonts/Arial.ttf"),
    ]
    bold_candidates = [
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
    avg_char = max(font_obj.getlength("abcdefghijklmnopqrstuvwxyz") / 26, 1)
    chars = max(int(width / avg_char), 10)
    for line in wrap(text, chars):
        draw.text((x, y), line, font=font_obj, fill=fill)
        y += font_obj.size + line_gap
    return y


def rounded_box(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int, int, int],
    *,
    fill: tuple[int, int, int],
    outline: tuple[int, int, int] | None = None,
    radius: int = 14,
) -> None:
    draw.rounded_rectangle(xy, radius=radius, fill=fill, outline=outline, width=2 if outline else 1)


def render_card() -> Image.Image:
    source = Image.open(SOURCE).convert("RGB")
    source_ratio = source.width / source.height
    target_ratio = W / H
    if source_ratio > target_ratio:
        new_h = H
        new_w = round(source.width * (H / source.height))
    else:
        new_w = W
        new_h = round(source.height * (W / source.width))
    bg = source.resize((new_w, new_h), Image.Resampling.LANCZOS)
    left = max((new_w - W) // 2, 0)
    top = max((new_h - H) // 2, 0)
    card = bg.crop((left, top, left + W, top + H)).convert("RGBA")

    veil = Image.new("RGBA", (W, H), (255, 255, 255, 0))
    veil_draw = ImageDraw.Draw(veil)
    for x in range(W):
        if x < 650:
            alpha = 242
        elif x < 1000:
            alpha = int(242 - ((x - 650) / 350) * 140)
        else:
            alpha = 86
        veil_draw.line([(x, 0), (x, H)], fill=(WHITE[0], WHITE[1], WHITE[2], alpha))
    card.alpha_composite(veil)

    draw = ImageDraw.Draw(card)
    title_font = font(112, weight="bold")
    lede_font = font(43, weight="bold")
    body_font = font(30)
    small_font = font(25, weight="bold")
    chip_font = font(23, weight="bold")

    # Brand mark.
    rounded_box(draw, (72, 66, 140, 134), fill=(255, 255, 255), outline=INK, radius=16)
    draw.line((103, 72, 103, 128), fill=(14, 124, 90, 88), width=5)
    draw.line((78, 102, 134, 102), fill=(29, 102, 209, 66), width=5)
    draw.text((162, 72), "DotMatch", font=font(34, weight="bold"), fill=INK)
    draw.text((164, 114), "guide counts and barcode QC", font=font(24), fill=MUTED)

    draw.text((72, 178), "DotMatch", font=title_font, fill=INK)
    lede_bottom = text_block(
        draw,
        (75, 300),
        "Count guides. Split barcodes. See what failed.",
        font_obj=lede_font,
        fill=INK,
        width=700,
        line_gap=6,
    )
    text_block(
        draw,
        (75, lede_bottom + 20),
        "FASTQ assignment for known short-DNA targets, with ambiguous reads kept visible.",
        font_obj=body_font,
        fill=MUTED,
        width=690,
        line_gap=10,
    )

    chips = [
        ("87,437 guides", GREEN),
        ("ambiguous reads visible", BLUE),
        ("331k reads/s", GREEN_2),
        ("HTML + TSV QC", AMBER),
    ]
    chip_positions = [(75, 494), (306, 494), (75, 552), (306, 552)]
    for (label, color), (x, y) in zip(chips, chip_positions):
        text_w = int(draw.textlength(label, font=chip_font))
        rounded_box(draw, (x, y, x + text_w + 36, y + 48), fill=(255, 255, 255), outline=LINE, radius=10)
        draw.rounded_rectangle((x + 15, y + 15, x + 28, y + 33), radius=4, fill=color)
        draw.text((x + 38, y + 12), label, font=chip_font, fill=INK)

    return card.convert("RGB")


def main() -> None:
    if not SOURCE.exists():
        raise SystemExit(f"Missing source art: {SOURCE}")
    card = render_card()
    for out in OUTPUTS:
        card.save(out, "PNG", optimize=True)
        print(out)


if __name__ == "__main__":
    main()
