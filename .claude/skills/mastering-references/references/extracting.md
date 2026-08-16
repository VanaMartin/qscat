# Extracting text, and getting the page numbers right

## Extract

```bash
uv run --with pypdf python -c "
from pypdf import PdfReader
src = 'reference/literature/<stem>.pdf'
r = PdfReader(src)
print('pages:', len(r.pages))
open(src.replace('.pdf', '.txt'), 'w').write(
    ''.join(f'\n===== PAGE {i+1} =====\n' + (p.extract_text() or '')
            for i, p in enumerate(r.pages)))
"
```

The `.txt` is gitignored, so it can sit beside the PDF as a working file.

## The page-offset check — do this before writing any locator

**The extractor's page index is not the printed page number.** Front matter,
cover sheets and offprint banners shift it, and a note full of locators that are
all off by three is worse than a note with none: the reader trusts it and lands
in the wrong place.

Find the offset once:

```bash
grep -n "===== PAGE" reference/literature/<stem>.txt | head -3
sed -n '/===== PAGE 1 =====/,/===== PAGE 2 =====/p' reference/literature/<stem>.txt | head -20
```

Read the printed page number off the first page of body text and compare it to
its extractor index. Then state the mapping at the top of your working notes and
apply it to every locator.

Two common cases in this collection:

- **APS journal articles** print per-article pages like `032721-1` … `032721-17`.
  Extractor page 1 is usually printed page `-1`, so the offset is zero and the
  locator is `p. 032721-6`.
- **Charles University theses** carry several unnumbered front pages, so printed
  page 1 is typically extractor page 9-11. Locators use the **printed** number.

Spot-check the offset again near the end of the document — a mid-document insert
can shift it.

## Finding what to extract

Work from the repository inward, not from the paper outward. The note carries
what the repo depends on, so start by asking what that is:

```bash
# who already cites this source, and for what
grep -rn "Houfek\|Rescigno\|McCurdy" --include="*.py" --include="*.md" \
  libs apps docs/physics CLAUDE.md | head -20
```

Then locate each of those claims in the source and anchor it. Anything the repo
does not lean on goes in "Not used here" as a one-liner, not in the body.

## Transcribing equations

Keep the paper's own numbering — that is the locator. Transcribe in plain-text
notation matching the repo's docstring style rather than LaTeX:

```
H_el(r; R) = -1/2 d^2/dr^2 + V(R, r)                    p. 032721-3, Eq. (5)
V_res(R)   = E_res(R) - (i/2) Gamma(R)                  p. 022714-8, Eq. (41)
```

State the units the paper uses, and whether they are the repo's atomic units.

## Checking parity with the repo

Before writing "matches", run the check:

```bash
grep -rn "12766.36\|0.75102\|1.1535" libs/qscat/qscat/model/library.py
```

Report the outcome either way. Two disagreements found this way so far — the
H₂⁺ reduced mass, and the thesis grid tables versus the eMoScat decks — were
both real, and both mattered.
