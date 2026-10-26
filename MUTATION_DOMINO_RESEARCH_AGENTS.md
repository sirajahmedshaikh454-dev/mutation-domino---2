# Mutation Domino 2.0 — Research-Grade Agent Execution Specification

> Mission: Extend the EXISTING Mutation Domino application into a reliable, evidence-aware, reproducible protein-variant research workbench. Preserve every working behavior. Do not claim biological predictions or scientific novelty without independent validation.

## 0. Non-negotiable operating contract

1. **Inspect before editing.** Read existing `AGENTS.md`, README, backend routes/models/services, frontend components, tests, migration history, and the most recent real report. Build a capability inventory and baseline test results. Never assume a feature is missing merely because it is not in this document.
2. **No destructive rewrite.** Do not scaffold a new project over the existing one, rename API contracts casually, delete data, overwrite SQLite databases, clear cached structure files, remove functioning UI, or replace working algorithms with approximations. Make incremental changes behind feature flags.
3. **Backward compatibility.** Existing endpoints, response fields, report fields, data records, user workflows, and export formats remain functional. Add schema versions, additive migrations, adapters, and deprecation notices if changes are needed. Existing reports must still be readable.
4. **Preserve the scientific boundary.** Existing calculations describe *wild-type/reference structure* and coarse residue-property differences. Do NOT claim mutant structure, altered contacts, ΔΔG, disease pathogenicity, drug binding, or clinical effects unless independently computed by a documented, validated method; label each prediction and its uncertainty. Never equate pLDDT with pathogenicity.
5. **No invented evidence.** No fabricated references, citations, clinical labels, confidence values, ground truth, assay results, benchmark scores, or successful test claims. Missing information must be displayed as missing.
6. **Zero-dollar, CPU-first.** The default local workflow must work on Windows/Git Bash with a 16 GB RAM, no GPU, and no paid APIs. Optional heavyweight or cloud services must be disabled by default. No mandatory authentication or LLM for local single-user operation.
7. **Safe engineering.** Do not bypass Windows Application Control, disable antivirus, or change security policy to load a blocked DLL. Diagnose, provide supported alternatives, and ask for approval for impactful changes.
8. **No silent failures.** Distinguish network errors, upstream 502, parse errors, unavailable structures, mapping failures, missing annotations, and unsupported scientific conclusions.
9. **Every phase is gated.** Implement, test, review, document, and only then proceed. Stop or report blockers instead of marking unverified tasks done.
10. **Human oversight.** Agent-generated scientific prose is informational, citation-linked, and never clinical advice. No real patient data in demos or public logs.

## 1. Baseline acceptance fixture — preserve EXACT behavior

Use existing P04637 human p53 R175H report as an immutable historical fixture, not as an unquestionable scientific gold standard. It documents:

- UniProt P04637, human p53, TP53; canonical R175H; structure chain A residue 175; 393/393 alignment; mapping status `mapped`.
- UniProt REST sequence and AlphaFold DB structure `AF-P04637-F1-model_v6.cif` (model v6); source URL, retrieved time, SHA-256, algorithm/schema versions.
- Mutation-site pLDDT 96.62 (very high); clearly local *reference structure* confidence.
- Arginine vs histidine coarse descriptors: basic/basic, positive/positive **as existing coarse labels**, aromatic false/true, Kyte–Doolittle -4.5/-3.2 (delta +1.3), mass 156.19/137.14 Da (delta -19.05), side-chain heavy atoms 7/6.
- 14 CA neighbors within 8 Å; 16 other residues with a heavy-atom pair within 5 Å; closest CA CYS176 3.81 Å, ARG174 3.82 Å, GLY245 5.03 Å; closest heavy-atom distances CYS176 1.33 Å, ARG174 1.34 Å.
- Existing PDF/JSON exports, scientific disclaimers, provenance, and UI workflows where implemented.

**Scientific improvement without breaking output:** Keep legacy raw heavy-atom distances and counts under versioned fields. Add a separate `nonbonded_contacts` view that excludes covalently bonded atom pairs, and mark sequence-adjacent residues. Do not silently change historical numbers. Flag histidine charge as pH-dependent in explanatory text while retaining historical coarse classification in backward-compatible data.

Regression tests should check the source/version and algorithm version before expecting exact numeric fixtures; allow source updates to create new versioned fixtures, never rewrite history. Independent coordinate checks must validate the fixture and reveal any prior bug rather than freezing it forever.

## 2. Intended users and real-life jobs-to-be-done

**A. Structural biologist:** upload/select a UniProt accession, inspect mutation site and nearby structural features, annotate evidence, compare variants, export reproducible research dossier.

**B. Industrial enzyme engineer:** analyze a *non-clinical* enzyme, annotate active sites, compare candidate substitutions, produce an experiment shortlist with suggested wet-lab assays (activity, thermostability, pH tolerance), not claims of improved yield.

**C. University instructor/student:** guided protein/mutation lesson with 3D visualization, definitions, and quiz; distinguish observed, predicted, and unknown.

**D. Research assistant:** batch screen a curated list of variants, inspect data provenance, filter by annotation coverage, download machine-readable outputs.

**E. Clinical genetics researcher:** optional public variant evidence exploration with prominent non-diagnostic labeling; no patient-specific treatment or pathogenicity classification generated by the tool.

Do not add generic features without a demonstrated user need. Prioritize A/B/C.

## 3. Architecture — extend existing stack

- Frontend: existing React + TypeScript + Vite + Tailwind/shadcn, Mol* viewer; use TanStack Query for server state only if already present or justified, and React context/CSS variables for light/dark/system theme.
- Backend: existing Python/FastAPI, Pydantic, Biopython, NumPy/SciPy, HTTPX, SQLAlchemy, Alembic.
- Database: SQLite local with WAL, foreign keys, indexed unique identifiers; optional PostgreSQL adapter for deployed multi-user workloads. Never require a hosted database for local use.
- External data: UniProt REST, AlphaFold DB; optionally curated InterPro, PDBe/SIFTS, RCSB PDB, ClinVar, Ensembl VEP, AlphaMissense, and PubMed/Crossref only through well-defined provider adapters, respecting licenses, terms, rate limits and coverage. Do not assume every service supports every protein or that all variant identifiers are directly interoperable.
- Scientific modules: `sequence_mapping`, `structure_parsing`, `contacts`, `confidence`, `chemistry`, `functional_annotations`, `evidence_graph`, `variant_comparison`, `reporting`, `validation`.
- Storage: source structures in content-addressed cache with atomic writes, checksums, timestamps, provenance; reports separate from DB; no binaries stored unnecessarily in SQLite.
- Background jobs: local bounded worker/thread/process queue for batch analyses, with persisted status and cancellation; avoid Celery/Redis unless scale requires it.
- Reproducibility: every analysis records accession, isoform, sequence checksum, structure accession/version/checksum, model type, parser/algorithm versions, parameters, source timestamps, warnings and status.

Data path: `React → FastAPI → validated provider adapters/cache → sequence/structure mapping → scientific feature engine → evidence/uncertainty layer → SQLite + reports → Mol* and charts`.

## 4. Scientific modules to implement

### 4.1 Verified mapping and structural quality

- Support canonical accession, selected isoform, chain, insertion codes, altlocs, missing residues, nonstandard residues, and multi-chain structures where supported.
- Require reference residue identity match; never map by index alone when numbering differs. Record ambiguous/unmapped results and reject invalid substitutions.
- Parse PDB and mmCIF with explicit model/chain selection; validate coordinates, residue uniqueness, and finite values.
- Distinguish AlphaFold predictions from experimentally determined structures. Clearly label structure origin, resolution if experimental, pLDDT if AlphaFold, and any missing confidence.
- Verify mapping with independent fixtures including indels, missing loops, multiple chains, mismatched reference residue, and structure revisions.

### 4.2 Contact science, not just distances

- Preserve legacy CA radius and minimum heavy-atom contacts.
- Add nonbonded contacts with documented exclusions (covalently bonded atoms; peptide bonds), sequence-separation filters, per-atom pair evidence, contact type *only where chemically supported*, and sensitivity across configurable thresholds.
- Distinguish contact geometry from hydrogen bonds, salt bridges, disulfide bonds, pi interactions, or energetic coupling; implement chemically valid rules with references or do not label them.
- Add residue solvent-accessible surface area (SASA) and buried/exposed status using a tested library or algorithm; document settings and limitations.
- Optional structural-site proximity: curated catalytic/binding/PTM sites from provenance-backed annotations; missing annotation != no site.
- Optional multi-structure comparison: experimental/reference models only when alignment and numbering can be validated; report RMSD and coverage honestly, no claim of mutation-caused change without mutant structure.

### 4.3 Mutation property descriptors

- Physicochemical deltas with source and scale (mass, charge under stated pH assumption, hydropathy, volume, polarity, aromaticity, proline/glycine special cases).
- Distinguish residue mass from free amino-acid mass; distinguish coarse categorical labels from computed protonation.
- No unsupported 'damage score'. If an exploratory priority index is added, label it **non-clinical heuristic**, document weights, sensitivity, missing-data policy, and lack of validation; never disguise as pathogenicity probability.

### 4.4 Evidence-aware variant graph — core research hypothesis

Define typed nodes: `variant`, `protein`, `residue`, `structural_region`, `functional_site`, `annotation`, `publication`, `structure`, `assay`.

Define typed edges: `maps_to`, `within_distance_of`, `annotated_as`, `supported_by`, `contradicted_by`, `derived_from`.

Every edge stores: `source_id`, source version/date, URL/DOI when available, evidence category (`experimental`, `curated`, `computational`, `geometric`, `hypothesis`), extraction method, confidence/quality *if meaningful*, and caveats. Geometric proximity is not proof of functional effect. Do not convert missing data into negative evidence.

Graph view must support click-through provenance, filtering by evidence type, contradiction visibility, and exporting graph JSON/GraphML. Keep graph computations deterministic. Use relational tables/JSON initially; no graph database is required.

### 4.5 Uncertainty and claim-audit engine — core research hypothesis

- Distinguish uncertainty from structure confidence, residue mapping, annotation completeness, conflicting publications, and model calibration.
- Define rules: claim types and minimum evidence; examples: `nearby` requires coordinates + threshold, `experimentally observed` requires assay/source, `pathogenic` requires external curated clinical interpretation and must be attributed rather than predicted.
- Before exporting prose, scan generated claims for unsupported causal language; require explicit evidence links or downgrade to hypotheses.
- Present concise natural-language explanations generated from deterministic templates. Optional local LLM explanations must be disabled by default, cite underlying records, and never invent sources.
- Add explanation audit logs and reviewer override with reason.

### 4.6 Variant comparison and batch analysis

- Side-by-side comparison (same protein and optionally different proteins with carefully normalized features); shared 3D camera where meaningful; per-variant evidence coverage and missing-data flags.
- Batch CSV import with strict schema, preview, deduplication, cancellation, retry, progress, per-row errors, and export; protect against giant inputs and formula injection in CSV exports.
- Precomputed AlphaMissense or other scores only where legitimately available and license-compliant; label as **third-party prediction**, separate from Mutation Domino's geometry and clinical labels.

### 4.7 Experiment planning — real-world impact

- For industrial enzymes: user-defined objective (temperature, pH, activity, substrate specificity), annotate known sites and suggest **tests** rather than asserting outcomes; example outputs: variant list, hypothesis, control, measurement, replicates, source rationale, failure modes.
- Offer prioritized *review queues* based on evidence completeness and user criteria, not fake predictions of enzyme improvement.
- Allow researchers to record their own measured assay results and compare hypotheses with observations. Use units, conditions, assay protocols, provenance, and immutable revisions.
- Export experiment worksheet CSV/PDF. No autonomous lab execution or hazardous wet-lab protocols.

### 4.8 Teaching mode and accessibility

- Explain DNA → gene → protein → folding → mutation with tooltips, Roman Urdu/English language switch, and contextual lessons; preserve technical terms.
- Interactive 3D site selection, guided R175H lesson, clear 'what this DOES / DOES NOT show' panels, colorblind-friendly palettes, keyboard navigation, WCAG contrast, and reduced motion.
- Responsive light/dark/system theme across Mol*, tables, exports, and graphs. Persist preference, respect OS changes, do not reset viewer or refetch structure on theme change.

## 5. Real engineering quality

### External API reliability

- Timeouts, bounded exponential backoff for transient 429/5xx with Retry-After, per-provider rate limiting, circuit breaker when appropriate, atomic cache, ETag/Last-Modified where supported, stale-cache warning, offline analysis of cached structures, and actionable errors.
- Avoid hardcoded AlphaFold file URLs where API metadata provides a canonical downloadable link; validate status, type, file size, coordinates, and checksums.
- Never turn upstream 502 into fabricated successful analysis.

### 3D viewer quality

- Mol* loads real CIF/PDB coordinates; shows visible cartoon, mutation highlight, neighbors, confidence coloring, residue labels, legend, selection sync, camera focus/reset, error/loading states, accessibility, and cleanup.
- Test rendering in a real browser with screenshots; a canvas existing in DOM is insufficient.

### Security and privacy

- Validate inputs, size limits, safe paths, upload sanitization, dependency pinning, license checks, request timeouts, safe URL allowlist for upstream fetches, CORS origin restriction, no secrets in logs, safe report file names.
- Optional accounts/teams only after local MVP is stable; privacy notice and deletion/export controls. No protected health information by default.

### Database evolution

Additive schema candidates: `source_records`, `structure_artifacts`, `analysis_runs`, `analysis_parameters`, `residue_features`, `evidence_nodes`, `evidence_edges`, `annotations`, `assay_observations`, `claim_audits`, `batch_jobs`, `report_artifacts`. Add Alembic migrations and backups; keep old IDs and reports valid. Use indexes on accession, variant key, source ID, run ID, status, created_at. Avoid duplicating entire large structures per run.

### API extensions (illustrative; adapt to existing conventions)

- `GET /api/v1/proteins/{accession}`
- `GET /api/v1/proteins/{accession}/structures`
- `POST /api/v1/analyses` and `GET /api/v1/analyses/{id}`
- `GET /api/v1/analyses/{id}/evidence`
- `GET /api/v1/analyses/{id}/claims`
- `POST /api/v1/comparisons`
- `POST /api/v1/batches`, `GET /api/v1/batches/{id}`
- `POST /api/v1/experiments`, `GET /api/v1/experiments/{id}`
- `GET /api/v1/reports/{id}.{format}`

Use OpenAPI/Pydantic models, versioning, pagination, errors with stable codes, and contract tests. Do not break existing routes.

## 6. Research contribution: hypotheses, not novelty claims

Proposed title: **Mutation Domino: An Evidence-Traceable and Uncertainty-Aware Workbench for Protein Missense Variant Interpretation**.

Candidate research questions:
- RQ1: Does provenance-linked evidence graph improve factual traceability and reduce unsupported biological claims compared with a conventional structural-neighborhood report?
- RQ2: Does explicit uncertainty labeling improve interpretation accuracy for nonexpert users?
- RQ3: Does the workflow reduce time-to-answer for researchers doing curated variant context analysis without reducing scientific correctness?

Do not claim the first-ever evidence graph, first AlphaFold mutation viewer, or first interpretable variant tool. Perform a literature review of AlphaFold DB, Missense3D, ProtVar, AlphaMissense, FoldX/Rosetta, and relevant evidence/variant platforms, with a dated feature-comparison matrix and verifiable sources. Exact novelty must be narrowed after literature search.

### Study design and benchmarks

- Curate versioned public protein/variant examples across multiple proteins, plus a held-out evaluation set. Record accession, isoform, mutation, source, annotations, experimental support, licenses, and exclusions. Avoid training/test leakage and duplicate variants across splits.
- Compare **B0** current Mutation Domino structural report, **B1** B0 + evidence graph, **B2** B1 + uncertainty/claim audit. Where meaningful, compare workflows with available existing tools using reproducible protocols and documented differences.
- Primary outcomes: unsupported-claim rate, evidence-link correctness, factual accuracy of generated explanations, missing-evidence handling, and time-to-completion of standardized tasks. Define operational scoring rubrics *before* testing.
- Secondary outcomes: mapping success rate, contact-calculation agreement against independent reference implementations, provenance completeness, API success/latency, user-rated comprehensibility, and usability.
- User study: consenting participants, predefined tasks, randomized/counterbalanced order, no real patient cases, anonymized responses, informed consent and institutional review/ethics check where required. Report sample size, expertise, uncertainty, effect sizes and confidence intervals.
- Statistics: pre-register hypotheses where possible, bootstrap confidence intervals, appropriate paired comparisons, multiple-testing correction when applicable; publish negative results too.
- Reproducibility: release code, environment lockfiles, test fixtures, data manifest and license-compatible evaluation scripts. Keep raw proprietary or restricted datasets out of repo.
- Scientific paper sections: Abstract, Introduction, Related Work, Research Questions, Methods, Data, Baselines, Experiments, Results, Error Analysis, Limitations, Ethics, Reproducibility, Conclusion.
- No guarantee of publication, accuracy, novelty, or superiority. Make no numeric performance claims until measured.

## 7. Agent roles and ownership

**Orchestrator / Release Engineer**: repo inventory, baseline, dependency graph, task planning, feature flags, backwards compatibility, merge gates, release notes.

**Structural Bioinformatics Agent**: sequence mapping, PDB/mmCIF, contact validation, chemistry, pLDDT, independent fixtures, scientific caveats.

**Evidence & Literature Agent**: provider adapters, data provenance, literature feature matrix, license review, evidence graph schema, contradiction handling.

**Backend & Data Agent**: FastAPI, migrations, caching, resilient HTTP, batch jobs, reports, OpenAPI, data integrity.

**Frontend & Mol* Agent**: interactive protein viewer, compare screens, graph, teaching mode, dark/light/system, accessibility and responsive behavior.

**Research Methods Agent**: hypotheses, dataset protocol, baselines, annotation rubric, experiments, statistics, reproducibility artifacts.

**QA / Security Agent**: regression suite, browser screenshots, API contract, property-based tests, security checks, performance benchmarks, release gate.

Each agent must record files changed, why, tests run, results, unresolved issues, and evidence. Orchestrator merges only passing, reviewed increments. Parallel work requires agreed API/data contracts first.

## 8. Phased execution (complete each gate before next)

**Phase 0 — Snapshot and protect (mandatory):** Git branch, backups, inspect project, record baseline run and R175H PDF/JSON, list passing tests and known failures. Deliver `docs/baseline.md`, `docs/architecture-current.md`. Gate: baseline reproducible, rollback documented.

**Phase 1 — Fix core reliability:** resolve blocked Biopython extension through approved environment/package approach (no security bypass), upstream 502 handling, structure caching, Mol* visible rendering, camera reset, pLDDT extraction, JSON/PDF tests. Gate: end-to-end P04637/R175H verified, error cases tested, no data loss.

**Phase 2 — Scientific integrity:** independent residue mapping and contacts validation, bonded vs nonbonded separation, confidence warnings, chemistry notes, versioned results. Gate: golden fixtures + independent calculations pass.

**Phase 3 — Evidence platform:** source adapters, graph, evidence filters, claim-audit engine, source-linked explanations. Gate: all edges traceable, missing data honest, adversarial unsupported claims rejected.

**Phase 4 — User workflows:** variant comparison, batch jobs, experiment planning, teaching mode, responsive UI and accessibility. Gate: representative researcher, engineer and student journeys pass browser E2E.

**Phase 5 — Research study:** literature review, frozen benchmark, baselines/ablations, blinded manual review where needed, user study with ethics review, stats and error analysis. Gate: results reproducible from scripts, no invented metrics.

**Phase 6 — Release:** Docker optional, local Git Bash one-command startup script, docs, demo, changelog, data license notices, reproducibility package and paper draft. Gate: clean install on supported Windows machine and CPU-only smoke test.

## 9. Definition of done

- Existing working features preserved and regression-tested; previous reports remain readable.
- P04637 R175H end-to-end works on the actual local stack: correct mapping, actual structure rendered, confidence correctly extracted, scientifically interpretable neighborhood, valid JSON/PDF.
- Added evidence graph and claim-audit are transparent, source-linked, and never silently infer causal biological effects.
- Real-world workflows (structural research, enzyme experiment planning, teaching) can be completed without fabricated data or paid services.
- Theme and 3D visualization work in both modes without state loss.
- Tests include unit, scientific property, integration, browser E2E, regression, upstream outage, malformed data, and migration/rollback scenarios.
- Research claims are limited to what controlled evaluation supports; novelty verified by literature comparison; paper draft contains real results or explicitly marked pending sections.
- Agent handoff includes exact changed files, run commands, test commands, logs/results, screenshots, scientific limitations, and known blockers.

## 10. Immediate command to the coding-agent team

> Read this specification AND the repository's existing AGENTS.md. Inspect the entire project and current tests. Begin with Phase 0; do not change code until baseline, backups, and existing API/report contracts are documented. Prioritize the existing Python/Biopython startup crash, any AlphaFold 502 failures, Mol* visibility, and confidence mapping. Then implement phased improvements behind feature flags, preserve all verified working features, and stop making success claims unless you ran the corresponding tests. Propose a research question and benchmark before claiming novelty. At each phase report the exact diff, tests, screenshots, blockers, and next steps. Continue as far as the environment allows without paid services, fake scientific results, or unsafe security bypasses.
