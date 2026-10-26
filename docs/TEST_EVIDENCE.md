# Test evidence

Date: 2026-10-26
Environment: Windows 11, Python 3.13.14, Node 22.19.0, npm 10.9.3

Commands:
- `cd backend && python -m pytest` -> 101 passed, 1 warning (after Phase 2.1 nonbonded atom-pair fix + 5 regression tests).
- `cd backend && ruff check . && ruff format --check .` -> All checks passed.
- `cd frontend && npm run lint` -> 0 errors, 0 warnings.
- `cd frontend && npm run typecheck` -> success.
- `cd frontend && npm run test -- --run` -> 1 file, 2 tests passed.
- `cd frontend && npm run build` -> success.

Smoke (fixture mode):
- Backend `/health` -> ok, mode fixture.
- Analysis P04637 R175H -> pLDDT 96.62, 14 CA neighbors, 16 heavy-atom contacts, 16 nonbonded contacts (was 14; adjacent residues retained via nonbonded pairs), mass delta -19.05.