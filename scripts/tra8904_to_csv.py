#!/usr/bin/env python3
"""Convert the TRA8904 ODS sheet to a rectangular CSV (uniform width, no trimming)."""
import zipfile
import csv
import xml.etree.ElementTree as ET

NS = {"table": "urn:oasis:names:tc:opendocument:xmlns:table:1.0",
      "office": "urn:oasis:names:tc:opendocument:xmlns:office:1.0",
      "text": "urn:oasis:names:tc:opendocument:xmlns:text:1.0"}
CELL = f"{{{NS['table']}}}table-cell"
ROW = f"{{{NS['table']}}}table-row"
REP_COLS = f"{{{NS['table']}}}number-columns-repeated"
REP_ROWS = f"{{{NS['table']}}}number-rows-repeated"
VAL = f"{{{NS['office']}}}value"
DATE = f"{{{NS['office']}}}date-value"


def cell_text(c):
    for a in (VAL, DATE):
        v = c.get(a)
        if v is not None:
            return v
    return " ".join("".join(p.itertext())
                    for p in c.iter(f"{{{NS['text']}}}p")
                    if "".join(p.itertext()).strip())


def main():
    src = "data/raw/traffic/tra8904-km-by-la.ods"
    dst = "data/raw/traffic/tra8904_full.csv"
    with zipfile.ZipFile(src) as z:
        root = ET.fromstring(z.read("content.xml"))
    rows = []
    for tbl in root.iter(f"{{{NS['table']}}}table"):
        if tbl.get(f"{{{NS['table']}}}name", "") != "TRA8904":
            continue
        for row in tbl.iter(ROW):
            rrep = int(row.get(REP_ROWS, "1") or 1)
            vals = []
            for cell in row.iter(CELL):
                crep = int(cell.get(REP_COLS, "1") or 1)
                vals.extend([cell_text(cell)] * crep)
            for _ in range(rrep):
                rows.append(vals)
    trimmed = []
    for r in rows:
        while r and r[-1] == "":
            r.pop()
        trimmed.append(r)
    width = max(len(r) for r in trimmed)
    with open(dst, "w", newline="") as fh:
        w = csv.writer(fh)
        for r in trimmed:
            w.writerow(r + [""] * (width - len(r)))
    print("rows:", len(rows), "width:", width)


if __name__ == "__main__":
    main()
