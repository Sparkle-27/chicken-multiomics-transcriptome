#!/usr/bin/env python3
from pathlib import Path
from pypdf import PdfReader,PdfWriter
import csv,json,sys
root=Path(sys.argv[1]).resolve()
files=[Path(x) for x in (root/'figures/pdf_order.txt').read_text().splitlines() if x]
assert len(files)==29 and len(set(files))==29
writer=PdfWriter();index=[]
for n,p in enumerate(files,1):
    reader=PdfReader(p)
    assert len(reader.pages)==1,(p,len(reader.pages))
    assert len((reader.pages[0].extract_text() or '').strip())>20,p
    writer.add_page(reader.pages[0])
    writer.add_outline_item(p.stem,n-1)
    index.append(dict(page=n,file=p.name,title=p.stem.replace('_',' ')))
writer.add_metadata({'/Title':'Chicken tissue-specific differential expression and QC','/Author':'Analysis generated from verified study counts','/Subject':'DESeq2 1.46.0; TS vs CB; main and sensitivity analyses'})
dest=root/'figures/gene_expression_figures.pdf'
with dest.open('wb') as f:writer.write(f)
assert len(PdfReader(dest).pages)==29
with (root/'figures/figure_index.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=['page','file','title']);w.writeheader();w.writerows(index)
print('Combined PDF:',dest,'pages:',len(index))
