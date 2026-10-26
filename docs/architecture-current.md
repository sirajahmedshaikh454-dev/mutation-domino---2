# Current architecture

Frontend: React 19 + TypeScript + Vite + Tailwind v4 + TanStack Query + Mol*.
Backend: FastAPI + Pydantic v2 + SQLAlchemy 2 + Alembic + SQLite + httpx + Biopython + NumPy + ReportLab.

Data flow:
accession -> UniProt adapter (or fixture) -> AlphaFold adapter (or fixture) -> bounded disk cache
-> sequence/structure mapping validation -> deterministic analysis engine (domain/) -> SQLite persistence
-> JSON/PDF exports -> Mol* viewer (wild-type only).

Key modules:
- backend/app/domain/*: pure scientific core (variant, structure, mapping, contacts, properties, analysis).
- backend/app/integrations/*: bounded HTTP, UniProt, AlphaFold, cache, fixtures.
- backend/app/services/*: orchestration.
- backend/app/persistence/*: SQLite tables/repositories/migrations.
- backend/app/exports/*: JSON and PDF.
- frontend/src: pages, components, API client, types, Mol* viewer.

Scientific boundary: wild-type structural context only; no mutant structure; pLDDT not mutation-effect confidence; no composite score.
