#!/usr/bin/env python3
"""Convert an ODS spreadsheet to CSV using only the Python stdlib.

Usage: ods_to_csv.py input.ods output.csv [sheet_name]
If sheet_name is omitted, all sheets are concatenated with a sheet column.
"""
import sys
import zipfile
import csv
import xml.etree.ElementTree as ET

NS = {
    "office": "urn:oasis:names:tc:opendocument:xmlns:office:1.0",
    "table": "urn:oasis:names:tc:opendocument:xmlns:table:1.0",
    "text": "urn:oasis:names:tc:opendocument:xmlns:text:1.0",
}

CELL = f"{{{NS['table']}}}table-cell"
ROW = f"{{{NS['table']}}}table-row"
TABLE = f"{{{NS['table']}}}table"
REP_COLS = f"{{{NS['table']}}}number-columns-repeated"
REP_ROWS = f"{{{NS['table']}}}number-rows-repeated"
VALUE_ATTR = f"{{{NS['office']}}}value"
DATE_ATTR = f"{{{NS['office']}}}date-value"


def cell_text(cell):
    for attr in (VALUE_ATTR, DATE_ATTR):
        v = cell.get(attr)
        if v is not None:
            return v
    parts = []
    for p in cell.iter(f"{{{NS['text']}}}p"):
        parts.append("".join(p.itertext()))
    return " ".join(x for x in parts if x.strip())


def main():
    src, dst = sys.argv[1], sys.argv[2]
    want_sheet = sys.argv[3] if len(sys.argv) > 3 else None
    with zipfile.ZipFile(src) as z:
        root = ET.fromstring(z.read("content.xml"))
    tables = root.iter(TABLE)
    out_rows = []
    for tbl in tables:
        name = tbl.get(f"{{{NS['table']}}}name", "")
        if want_sheet is not None and name != want_sheet:
            continue
        for row in tbl.iter(ROW):
            rrep = int(row.get(REP_ROWS, "1") or 1)
            vals = []
            for cell in row.iter(CELL):
                crep = int(cell.get(REP_COLS, "1") or 1)
                v = cell_text(cell)
                vals.extend([v] * crep)
            for _ in range(rrep):
                out_rows.append([name] + vals if want_sheet is None else vals)
    if not out_rows:
        sys.exit(f"sheet not found: {want_sheet}")
    trimmed = []
    for r in out_rows:
        r = list(r)
        while r and r[-1] == "":
            r.pop()
        trimmed.append(r)
    width = max(len(r) for r in trimmed)
    with open(dst, "w", newline="") as fh:
        w = csv.writer(fh)
        for r in trimmed:
            w.writerow(r + [""] * (width - len(r)))
    print(f"wrote {len(out_rows)} rows to {dst}")


if __name__ == "__main__":
    main()
