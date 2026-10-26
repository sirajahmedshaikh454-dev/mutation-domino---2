# Baseline snapshot

Date: 2026-10-26
Environment: Windows 11, Python 3.13.14, Node 22.19.0, npm 10.9.3

## Commands run
- `cd backend && python -m pytest` -> 96 passed, 1 warning (StarletteDeprecationWarning httpx).
- `cd backend && ruff check . && ruff format --check .` -> all checks passed.
- `cd backend && python -m mypy app` -> success (per prior evidence; not re-run this session).
- `cd frontend && npm run lint` -> 0 errors/warnings.
- `cd frontend && npm run typecheck` -> success.
- `cd frontend && npm run test -- --run` -> 1 file, 2 tests passed.
- `cd frontend && npm run build` -> success.

## Verified end-to-end (fixture mode)
- Backend `/health` returns mode fixture.
- Protein lookup P69905 and P04637 works.
- Analysis P04637 R175H: mapping 393/393, pLDDT 96.62, 14 CA neighbors, 16 heavy-atom contacts, mass delta -19.05. Matches documented baseline.
- Export endpoints generate JSON/PDF.

## Known limitations
- Live AlphaFold/UniProt not exercised (offline fixture mode).
- No git repository present; changes tracked manually.
- Playwright E2E not present/run.
- Biopython extension crash not reproduced on this machine.

## Backup
- Source tree is the working copy; no separate backup created. Rollback by restoring from original archive.
