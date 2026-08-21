from pathlib import Path
from PIL import Image, ImageOps, ImageDraw


ROOT = Path(__file__).parent
PREFIXES = [
    "mult-main",
    "mult-methods",
    "mult-supp",
    "fs-main",
    "fs-methods",
    "fs-supp",
]


for prefix in PREFIXES:
    paths = sorted(ROOT.glob(f"{prefix}-*.png"))
    if not paths:
        continue
    thumbs = []
    for path in paths:
        image = Image.open(path).convert("RGB")
        target_width = 850
        target_height = round(image.height * target_width / image.width)
        image = image.resize((target_width, target_height), Image.Resampling.LANCZOS)
        canvas = Image.new("RGB", (target_width + 20, target_height + 50), "white")
        canvas.paste(image, (10, 40))
        ImageDraw.Draw(canvas).text((10, 10), path.stem, fill="black")
        thumbs.append(canvas)

    columns = 2
    rows = (len(thumbs) + columns - 1) // columns
    cell_width = max(image.width for image in thumbs)
    cell_height = max(image.height for image in thumbs)
    sheet = Image.new("RGB", (columns * cell_width, rows * cell_height), "#d8d8d8")
    for index, image in enumerate(thumbs):
        x = (index % columns) * cell_width
        y = (index // columns) * cell_height
        sheet.paste(image, (x, y))
    output = ROOT / f"{prefix}-contact.png"
    sheet.save(output, optimize=True)
    print(output)
