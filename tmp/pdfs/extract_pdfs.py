from pathlib import Path
from pypdf import PdfReader


PDFS = [
    Path(r"C:\Users\willi\Downloads\MultifSuSiE_Manuscript-7.pdf"),
    Path(r"C:\Users\willi\Downloads\2025.08.17.670732v1.full.pdf"),
]


for pdf_path in PDFS:
    reader = PdfReader(pdf_path)
    output_path = Path(__file__).parent / f"{pdf_path.stem}.txt"
    with output_path.open("w", encoding="utf-8") as output:
        for page_number, page in enumerate(reader.pages, start=1):
            output.write(f"\n===== PDF PAGE {page_number} =====\n")
            output.write(page.extract_text() or "")
            output.write("\n")
    print(f"{pdf_path.name}: {len(reader.pages)} pages -> {output_path}")
