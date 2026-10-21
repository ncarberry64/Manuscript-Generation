"""Render every final page and compact contact sheets for visual review."""
from pathlib import Path
import hashlib
import json
import pypdfium2 as pdfium
from PIL import Image, ImageOps, ImageDraw

root = Path(__file__).resolve().parents[3]
pdf = root / 'manuscript/strong_lensing_main.pdf'
out = Path(__file__).resolve().parent
cache = root / 'data/strong_lensing_response_completion_20261021/native_cache/pdf_review'
cache.mkdir(parents=True, exist_ok=True)
doc = pdfium.PdfDocument(str(pdf))
thumbs = []
for i, page in enumerate(doc):
    im = page.render(scale=1.7).to_pil().convert('RGB')
    im.save(cache / f'page_{i+1:02d}.png')
    thumb = ImageOps.contain(im, (360, 485))
    cell = Image.new('RGB', (384, 520), 'gray')
    cell.paste(thumb, ((384-thumb.width)//2, 18))
    ImageDraw.Draw(cell).text((12, 500), f'Page {i+1}', fill='white')
    thumbs.append(cell)
for start in range(0, len(thumbs), 12):
    sheet = Image.new('RGB', (384*4, 520*3), 'gray')
    for j, thumb in enumerate(thumbs[start:start+12]):
        sheet.paste(thumb, ((j%4)*384, (j//4)*520))
    sheet.save(out / f'contact_{start+1:02d}_{min(start+12,len(thumbs)):02d}.jpg', quality=90)
(out/'render_receipt.json').write_text(json.dumps({
    'pdf_sha256': hashlib.sha256(pdf.read_bytes()).hexdigest(),
    'pages': len(doc), 'all_pages_rendered': True, 'scale': 1.7,
    'page_directory': str(cache), 'visual_review_complete': False,
}, indent=2)+'\n')
print(f'Rendered {len(doc)} pages, 4 contact sheets')
