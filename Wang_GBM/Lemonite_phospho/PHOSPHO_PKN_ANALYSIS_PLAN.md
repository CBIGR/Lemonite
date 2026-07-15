# Plan — Phosphosite prior-knowledge evaluation against modules

**Goal.** Do for **phosphosite regulators** what the pipeline already does for metabolites:
test whether a phosphosite assigned to a module has *known* regulatory connections to
that module's genes (and, specifically, to the module's **TF regulators**), using a
prior-knowledge network (PKN). Analogous to `calculate_metabolite_gene_enrichment` in
`nextflow/scripts/evaluate_against_PKN.py`.

Status legend: ☐ todo · ◐ in progress · ☑ done

---

## 0. Context established (investigation done)

- **The current PKN has NO phosphosite edges.** `nextflow/PKN/Lemonite_PKN.tsv`
  (`Node1, Node2, Source, Type`) contains only `PPI` (2,540,152 edges) and
  `metabolite-gene` (381,227 edges). So the phospho→gene prior knowledge must be
  **built and integrated first** — this is the main new piece of work.
- **How the pipeline validates metabolites** (the template to mirror,
  `evaluate_against_PKN.py:499`): build `known_interactions = {(regulator_id, gene)}`
  from the PKN; per module, count observed (regulator × module-gene) pairs that are
  known vs. expected from global density; hypergeometric test + BH-FDR. Output per-module
  fold-enrichment / p / FDR, and draw subnetworks with the supported edges.
- **How a new edge Type is added to the PKN** (template: `build_PKN/…/step2b_hPTM/` +
  `step3_final/combiner.py`): a self-contained module fetches raw prior knowledge, writes
  `<type>_network.tsv`, and the combiner standardizes it to `Node1/Node2/Source/Type` and
  concatenates. hPTM there uses `Type = histone-modification` with `Node1=Enzyme,
  Node2=Mark`. We add phospho the same way.
- **Phosphosite ID structure:** `GENE:RefSeqProtein:Residue[;Residue…][.n]`, e.g.
  `ANK2:NP_001139.3:S1858`. The **phosphoprotein gene symbol = token before the first
  `:`**. (Already used by `prepare_phospho_regulator.py`.)
- **Resources available:**
  - `.sif` has **OmnipathR + decoupleR + dorothea** (R) → kinase–substrate and enzyme–PTM
    networks (needs network access at build time — see risk R1).
  - `nextflow/PKN/CollecTRI_network.txt` → `source(TF) → target(gene)` with mode-of-reg
    (1,175 TFs). Local, no network needed.
  - `nextflow/PKN/Lemonite_PKN.tsv` PPI layer → phosphoprotein–gene physical interactions.
  - Selected phospho regulators live in each variant's
    `run/Networks/Phospho2targets_*.txt` (`Regulator`, pipe-separated `Targets`) and
    `run/ModuleViewer_files/Phospho.selected_regs_list.txt`.
- **Signal is already visible without any new PKN:** of the 209 phosphoproteins selected
  in the top-2000 variant, **11 are TFs in CollecTRI** (SOX2, SOX10, PML, WWTR1/TAZ, NFIA,
  MYT1, MEOX2, IRX1, KCNIP3, SFPQ, SOX8) — several are core GBM TFs. So the TF-connection
  question is likely to yield interpretable hits immediately.

---

## 1. What "a known phosphosite→gene edge" means

A phosphosite sits on a phosphoprotein `P` (gene symbol). We define three
biologically-motivated edge types to test, each with a source (independent, can be run/
scored separately, and combined):

| Edge type | Node1 → Node2 | Meaning | Source |
|---|---|---|---|
| **A. kinase→substrate** | kinase `P` → substrate gene | phospho-event on a kinase that phosphorylates a module gene | **OmniPath enzsub** via `decoupleR::get_ksn_omnipath()` (aggregates PhosphoSitePlus, SIGNOR, DEPOD, dbPTM, phospho.ELM, HPRD, NetworKIN… — see PHOSPHO_PKN_RESOURCES.md) |
| **B. phospho-TF→target** | TF `P` → target gene | phospho-event on a TF that transcriptionally regulates a module gene | CollecTRI (`source=P`), local |
| **C. phosphoprotein PPI** | `P` → interacting gene | phospho-event on a protein that physically interacts with a module gene | existing PPI layer of `Lemonite_PKN.tsv`, local |
| *(optional)* D. Kinase Library | kinase → phosphosite@residue | motif-predicted kinase for our exact phosphosites | Johnson atlas (PhosphoSitePlus); site-resolved; **non-commercial DL — do not auto-scrape** |

**Resource research done 2026-07-03 → see `PHOSPHO_PKN_RESOURCES.md`.** Key finding: pulling
the **OmniPath enzyme-substrate KSN** (edge type A) in one shot subsumes ~18 curated+predicted
databases, so it is the highest-value first fetch and is license-clean + already in the `.sif`.
Only the Johnson **Kinase Library** atlas is genuinely additive beyond OmniPath, but its bulk
download is non-commercial-gated (needs the user's manual academic download — flagged, not
auto-fetched).

Notes / decisions to make (see §6 open questions):
- Edge is defined at the **phosphoprotein (gene) level**, not the exact residue — public
  substrate/TF resources are gene-level, not site-resolved for arbitrary substrates.
  (Optionally, kinase-substrate resources ARE site-resolved on the *substrate* side; we
  can keep residue info as metadata but match at gene level for enrichment.)
- **B (phospho-TF→target)** is the most direct answer to the user's "connections between
  phosphosites and TF regulators" question, and needs only local CollecTRI. Start here.

---

## 2. Build the phospho prior-knowledge network  ☐

Mirror `step2b_hPTM`. New helper(s) in this folder (kept separate from the pipeline):

- `build_phospho_pkn.py` → writes `phospho_pkn.tsv` with columns
  `Node1  Node2  Source  Type  [Residue]  [url]`, where:
  - **B (CollecTRI):** for each phosphoprotein `P` that is a CollecTRI `source`, emit
    `P → target` rows, `Source=CollecTRI`, `Type=phospho-TF-target`. *(local, do first)*
  - **A (OmniPath enzsub):** `import_omnipath_enzsub()` in the `.sif`; keep
    `enzyme_genesymbol → substrate_genesymbol`, `Type=kinase-substrate`,
    `Source=OmniPath:<resources>`, residue = `substrate_residue_type+offset`.
  - **C (PPI):** not rebuilt — read directly from `Lemonite_PKN.tsv` (`Type=PPI`) at
    evaluation time, restricted to phosphoprotein nodes.
- Cache the OmniPath pull to `phospho_pkn_cache/` so it is reproducible offline once
  fetched (network is the only external dependency — see R1).

**Integration — DECIDED (user, 2026-07-03): standalone supplementary PKN.**
Keep `phospho_pkn.tsv` (this folder) separate and load it only in the phospho evaluation.
The shared `nextflow/PKN/Lemonite_PKN.tsv` is NOT modified. Same `Node1/Node2/Source/Type`
schema as the main PKN so the evaluation code can treat it identically.

Rejected alternative (recorded for context):
- *Merge into the PKN* — append to `Lemonite_PKN.tsv` with the new Types via the
  same standardize-and-concat pattern as `combiner.py`. More faithful to the pipeline but
  edits a large shared, gitignored artifact.

---

## 3. Evaluate phospho regulators against modules  ☐

New `evaluate_phospho_against_pkn.py` (adapts `calculate_metabolite_gene_enrichment`):

1. Load, per variant, module→genes (`tight_clusters.txt` filtered to the 46 network
   modules / `specific_modules.txt`) and module→phospho-regulators
   (`Phospho2targets_*.txt` / `Phospho.selected_regs_list.txt`).
2. Map each phospho regulator → phosphoprotein gene `P`.
3. Build `known_interactions = {(P, gene)}` from `phospho_pkn.tsv` (+ PPI subset).
4. Per module: observed = #(P, module-gene) pairs in `known_interactions`; expected from
   global density; **hypergeometric SF + BH-FDR**; report fold-enrichment.
5. Emit `Phospho_PKN_evaluation.{csv,txt}`: per-module table + supported (P → gene) edge
   list, tagged by edge type (A/B/C), plus a `.mvf` for ModuleViewer (mirror the metabolite
   `metabolite_LemoniteKG_interactions.mvf`).
6. **TF-specific view:** for each module, list phospho-TF→target support where the target
   is a module gene AND (bonus) the TF's phospho is a selected regulator of *that* module —
   i.e. "this module's phospho-TF regulator has known targets among the module's genes."
   Also cross-link to the module's existing **TF (Lovering) regulators**: does a phospho
   event sit on a protein known to regulate / be regulated by the module's assigned TFs?

Run for **both variants** (`top2000_variable`, `all_phosphosites`), scoring each edge type
separately and combined, so we can see which layer drives the signal.

---

## 3b. Visualize phospho interactions on heatmaps + subnetworks  ☐  (user request)

Show the *known* phospho interactions (from the phospho-PKN we integrate) on the module
figures, exactly the way metabolite/PPI/HumanNet interactions are already shown. Mechanism
traced from the existing pipeline:

- **`.mvf` interaction files** drive the heatmaps. `evaluate_against_PKN.py` writes three:
  `metabolite_LemoniteKG_interactions.mvf`, `PPI_interactions.mvf`,
  `HumanNet_interactions.mvf` into `ModuleViewer_files/`. Format (verified):
  ```
  ::TYPE=Lemonite_KG
  ::TITLE:Lemonite_KG
  ::OBJECT=GENES
  ::COLOR=YELLOW
  <module_id>\t<gene1>|<gene2>|…\t<regulator_name>
  ```
  Each line = "in module M, regulator R has known interactions with these module genes."
  `module_viewer.py` loads each by **explicit filename** (`load_metabolite_interactions`,
  `load_ppi_interactions`) and draws them as annotation tracks / arcs.
- **Subnetwork figures:** `draw_subnetwork()` (`evaluate_against_PKN.py:205`) builds a
  per-module graph, adding regulator→target and PPI/HumanNet edges categorized by source
  with per-category colors (see `get_edge_source_and_type` + the color map ~line 349).

**Phospho additions (mirror the above, do NOT duplicate metabolite logic):**
1. `evaluate_phospho_against_pkn.py` writes **`phospho_LemoniteKG_interactions.mvf`** (same
   header/format; `COLOR` = a new phospho color, e.g. `ORANGE`), one line per
   (module, phosphosite regulator) with the module genes it has a known phospho-PKN edge to.
   Also split by edge type if useful (`phospho_KSN_…`, `phospho_TF_…`).
2. **Heatmaps:** add a `load_*` + draw call for the phospho `.mvf` in `module_viewer.py`
   (new track/arc set, own color + legend entry). Small, additive edit — keep the existing
   metabolite/PPI/HumanNet blocks untouched.
3. **Subnetworks:** feed the phospho-PKN edges into `draw_subnetwork` as a new edge category
   (`phospho` / `kinase-substrate` / `phospho-TF-target`) with its own color + legend, so
   phospho regulator→target edges render alongside metabolite/PPI edges.
4. Re-run the heatmap + subnetwork steps for both variants so the figures include phospho.

**Decision to confirm:** edit the shared pipeline scripts (`module_viewer.py`,
`evaluate_against_PKN.py`) on this branch, vs. keep phospho-viz in separate wrapper scripts
in `Lemonite_phospho/`. Leaning **separate wrappers** (or a thin monkey-patch) to stay
consistent with "standalone / don't touch shared pipeline"; but a clean, guarded edit to the
shared scripts is acceptable on the branch since it's isolated. → default: separate wrappers
that reuse the shared functions, falling back to a guarded edit if the functions aren't
importable.

---

## 4. Report / interpret  ☐

- Per-variant summary: # phospho regulators with ≥1 known edge to their module; # modules
  with FDR<0.05 phospho-PKN enrichment; breakdown by edge type A/B/C.
- Highlight the GBM-relevant phospho-TFs (SOX2/SOX10/WWTR1/PML…) whose module targets are
  PKN-supported — these are the mechanistically interpretable findings.
- Append a "Phospho–PKN evaluation" section to `README.md`.

---

## 5. Sanity checks before trusting results  ☐

- Gene-symbol namespace: phosphoprotein symbols vs. module gene symbols vs. CollecTRI/
  OmniPath symbols must be reconciled (aliases). Spot-check overlap counts.
- Null/background: confirm the hypergeometric background (all module genes × all selected
  phosphoproteins) matches the metabolite convention so numbers are comparable.
- Compare against a shuffled-regulator control to confirm enrichment isn't trivial.

---

## 6. Open questions for the user

- **Q1. Integration — RESOLVED (2026-07-03): standalone supplementary phospho-PKN.** Do NOT
  modify `Lemonite_PKN.tsv`.
- **Q2. Edge types:** start with B (phospho-TF→target, local/CollecTRI) only, or build all
  three (A kinase-substrate needs an OmniPath network fetch)? *(still open — leaning: build
  B first since it's local and answers the TF question, then attempt A.)*
- **Q3. Matching granularity:** gene-level (phosphoprotein↔gene) — recommended, matches the
  resources — or attempt residue-resolved where available? *(still open)*
- **Q4. Direction of the TF question:** (a) phospho *on a TF* whose targets are module genes,
  and/or (b) phospho on a kinase that regulates the module's assigned TF regulators? *(still
  open)*

## Risks

- **R1 (network):** OmniPath enzsub (edge type A) requires internet at build time.
  **CONFIRMED reachable 2026-07-03** — `curl -I https://omnipathdb.org/` → 200 from the
  sandbox; curl/wget present. B+C are fully local regardless. Kinase Library (D) is
  non-commercial-gated → not auto-fetched (user downloads manually if wanted).
- **R2 (symbol drift):** RefSeq-based phosphosite gene symbols may differ from HGNC symbols
  used elsewhere; needs an alias map (the PKN build already ships `ensembl_mapping_*`).

---

---

# ▶ RESUME HERE (read this first next week)

## ✅ STATUS 2026-07-04: phospho-PKN build + evaluation + visualization COMPLETE for both
variants. Full results in **`PHOSPHO_PKN_REPORT.md`** (read that first). Q1–Q4 all resolved
(Q3=gene-level, Q4=phospho-on-TF, both cross- and same-module reported). Remaining OPTIONAL
next steps only: Johnson Kinase Library atlas (manual academic DL), shuffled-regulator null,
promote builder to a `step2c_phospho/` pipeline package. Nothing committed/pushed.

**Branch:** `phospho-regulator-analysis` (created 2026-07-03 off `main`; nothing committed
or pushed — all work is uncommitted/untracked in the working tree). Check with
`git -C /home/borisvdm/repo/LemonIte branch --show-current`.

**Decisions locked:**
- Standalone supplementary phospho-PKN — do NOT touch `nextflow/PKN/Lemonite_PKN.tsv`.
- Same `Node1/Node2/Source/Type` schema as the main PKN.

**Resource research done → `PHOSPHO_PKN_RESOURCES.md`.** Build plan: **A (OmniPath KSN) +
B (CollecTRI) + C (PPI)** into one standalone `phospho_pkn.tsv`. OmniPath's KSN alone
subsumes ~18 curated+predicted DBs (PhosphoSitePlus/SIGNOR/DEPOD/dbPTM/phospho.ELM/HPRD/…),
so it is the main external fetch and is license-clean. Kinase Library atlas (D) is
non-commercial-gated — leave for a manual user download, don't auto-scrape.

**Still need user input (non-blocking):** Q3 (granularity), Q4 (TF direction) — §6. Sensible
defaults chosen if no answer: gene-level matching (Q3), TF-question direction (a) (Q4).

**Immediate next action — `build_phospho_pkn.py` (this folder) producing standalone
`phospho_pkn.tsv`** with `Node1,Node2,Source,Type[,Residue]`:
- **A (OmniPath KSN):** in the `.sif`,
  `decoupleR::get_ksn_omnipath()` (or `OmnipathR::enzyme_substrate()` filtered to
  phosphorylation) → `Node1=enzyme_genesymbol Node2=substrate_genesymbol
  Source=OmniPath:<resources> Type=kinase-substrate` (keep substrate residue as metadata).
  Cache the raw pull to `phospho_pkn_cache/` for offline reproducibility.
- **B (CollecTRI):** `nextflow/PKN/CollecTRI_network.txt` (`source,target,mor`) →
  `Node1=TF Node2=target Source=CollecTRI Type=phospho-TF-target`.
- **C (PPI):** read from `Lemonite_PKN.tsv` (`Type=PPI`) at eval time, restricted to
  phosphoprotein nodes (not copied into phospho_pkn.tsv).
Then `evaluate_phospho_against_pkn.py` adapting `evaluate_against_PKN.py:499`
(hypergeometric per module), run for both variants, score A/B/C separately + combined.

**Optional (bigger, pipeline-native):** add a `step2c_phospho/` package to
`build_PKN/PKN_build_pipeline_opus4.8/` mirroring `step2b_hPTM/` (resumable via
`utils/cache.py`, 429-aware via `utils/api_retry.py`) so the phospho-PKN build is a
first-class, monthly-refreshable pipeline step. Standalone output only (do not merge into
the shared PKN).

**Overnight-run readiness (2026-07-03):** `.claude/settings.local.json` now allowlists
WebSearch/WebFetch, curl/wget, archive tools, and git branch/checkout/add (commit & push
still DENIED). Network to omnipathdb confirmed. ⚠️ Permission-rule changes only take effect
after a config reload — if commands prompt mid-run, open `/hooks` once or restart before
leaving it unattended, else an unattended build could stall on a prompt.

**Key paths (all verified this session):**
- Metabolite-eval template: `nextflow/scripts/evaluate_against_PKN.py:499`
  (`calculate_metabolite_gene_enrichment`).
- PKN-build / new-edge-type template: `build_PKN/PKN_build_pipeline_opus4.8/step2b_hPTM/`
  + `step3_final/combiner.py`.
- Main PKN (PPI + metabolite-gene only, NO phospho): `nextflow/PKN/Lemonite_PKN.tsv`.
- TF resource: `nextflow/PKN/CollecTRI_network.txt` (1,175 TFs).
- Selected phospho regulators per variant:
  `…/noProteomics_percentile2_divide_by_sum_phospho/<variant>/run/Networks/Phospho2targets_*.txt`
  and `…/run/ModuleViewer_files/Phospho.selected_regs_list.txt`; variants =
  `top2000_variable`, `all_phosphosites`.
- Run tools inside `.sif`: `nextflow/lemontree-pipeline_v1.0.0.sif` (has OmnipathR,
  decoupleR, dorothea, pandas). A local python env with pandas/numpy/openpyxl:
  `/home/borisvdm/Software/miniconda3/envs/LemonIte/bin/python3`.

**Facts to reuse (don't re-derive):**
- Phosphosite ID = `GENE:RefSeq:Residue…`; phosphoprotein gene = token before first `:`.
- top-2000 variant: 209 unique selected phosphoproteins, of which **11 are CollecTRI TFs**:
  SOX2, SOX10, PML, WWTR1, NFIA, MYT1, MEOX2, IRX1, KCNIP3, SFPQ, SOX8.
- Modules in the phospho-augmented network: **46** (coherence_threshold 0.5).

---

## Progress log

- 2026-07-03 (session 1) — Investigation complete: confirmed no phospho edges in PKN;
  identified metabolite-eval template (`evaluate_against_PKN.py:499`) and PKN-build template
  (`step2b_hPTM` + `combiner.py`); confirmed OmnipathR/decoupleR/dorothea in `.sif` and
  CollecTRI locally; found 11/209 selected phosphoproteins are CollecTRI TFs. Permissions
  scoped to repo + GBM data/results (no git commit/push) in `.claude/settings.local.json`.
  Plan drafted.
- 2026-07-03 (session 1) — User decisions: **standalone phospho-PKN** (§2 updated, Q1
  resolved). Created branch **`phospho-regulator-analysis`** to isolate this work (nothing
  committed/pushed). Added the "RESUME HERE" section above.
- 2026-07-03 (session 1) — **Online resource review** → new `PHOSPHO_PKN_RESOURCES.md`.
  Decided A+B+C build: **OmniPath KSN** (`decoupleR::get_ksn_omnipath()`) is the key add
  (subsumes ~18 DBs incl. PhosphoSitePlus/SIGNOR/DEPOD/dbPTM/phospho.ELM/HPRD/NetworKIN);
  CollecTRI (B) + PPI (C) local; Kinase Library atlas (D) is non-commercial-gated → manual
  DL only. Network to omnipathdb **confirmed 200** from sandbox. Extended
  `.claude/settings.local.json` for an **unattended overnight run** (WebSearch/WebFetch,
  curl/wget, tar/gunzip/unzip, git branch/checkout/add; commit+push still denied). Plan +
  RESUME updated. **Next:** implement `build_phospho_pkn.py` (A via `.sif`, B, C) → run
  `evaluate_phospho_against_pkn.py` on both variants. ⚠️ new perm rules need a config reload
  (open `/hooks` or restart) to be guaranteed live before leaving it running.
