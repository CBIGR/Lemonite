# TFA sample mismatch in `LemonPreprocessed_complete.txt`

## Symptom

In a sepsis run, the ATF5 row differs between the two preprocessed expression files, and the
values in `complete` are the TFA values but placed under the wrong sample names:

```
LemonPreprocessed_expression.txt   symbol  ensembl_gene_id  C1       C2        C3
                                   ATF5    ENSG00000169136  0.11866  -1.17917  -0.5241

LemonPreprocessed_complete.txt     symbol  ensembl_gene_id  C1       C2        C3
                                   ATF5    ENSG00000169136  0.57868  0.56967   0.52725

TFA_consensus.txt                          C1       C10      C11
                                   ATF5    0.57868  0.56967  0.52725
```

`complete.txt` stores samples in the count-matrix order (`C1 C2 C3 C4 C5 P1 …`), while
`TFA_consensus.txt` stores them alphabetically (`C1 C10 C11 C12 …`). The first three TFA
values were copied into the first three columns of `complete.txt` regardless of those labels,
so `C2` holds C10's activity and `C3` holds C11's.

## Root cause

`scripts/Preprocessing_TFA_RNA.R`, in the block that replaces TF expression with TF activity:

```r
RNA_preprocessed_withTFA[lovering_tfs_in_hvg, ] <- TFA_df_lovering[lovering_tfs_in_hvg, , drop=FALSE]
```

Both sides are named data frames, so this looks name-safe. It is not: R's
`df[rows, ] <- other_df[rows, ]` copies **by column position** and ignores the column names of
the right-hand side.

```r
a <- data.frame(C1=1, C2=2, C3=3, row.names="ATF5")
b <- data.frame(C1=10, C10=20, C11=30, row.names="ATF5")
a["ATF5", ] <- b["ATF5", , drop=FALSE]
#      C1 C2 C3
# ATF5 10 20 30    <- C2 got b's C10 value, C3 got b's C11 value
```

decoupleR returns its conditions sorted alphabetically, which for sample IDs like `C1…C14`,
`P1…P9` is a different order than the count matrix. Every column where the two orders diverge
gets the wrong sample's activity — for the sepsis layout that is nearly all of them.

`TFA_consensus.txt` itself is correct. It is built by `pivot_wider` from decoupleR's long
`(condition, source, score)` output, so labels stay attached to their values; only the merge
into `RNA_preprocessed_withTFA` breaks the link.

## Scope

- Affects **TFs that are also highly variable genes**, which take the assignment above.
- TFs added by `rbind()` a few lines later are **not** affected: `rbind()` on data frames does
  match columns by name.
- Only `LemonPreprocessed_complete.txt` is corrupted. `LemonPreprocessed_expression.txt` is
  built from `RNA_preprocessed_noTFA`, which never touches the TFA matrix.
- Downstream: tight clustering reads `LemonPreprocessed_expression.txt`, so **modules are
  unaffected**. Regulator assignment reads `complete.txt`, so **TF regulator assignment and
  anything derived from it are affected**. Metabolite and protein regulators are unaffected.

## Verification

Reproduced with decoupleR 2.12.0 on a synthetic matrix using the real sepsis sample order, by
running the pipeline's own code path:

```
expression matrix column order : C1 C2 C3 C4 C5 P1 ...
TFA_df        column order      : C1 C10 C11 C12 C13 C14 ...
TFA_df sorted alphabetically?   : TRUE

WITHOUT FIX -> TF row equals name-matched TFA:
    ATF5  (HVG, positional assignment) : FALSE
    FOXN3 (added via rbind)            : TRUE
WITH FIX:  all TRUE
```

The corrupted row shows the same off-by-label shift as the real data: the value written for
`C2` is the activity belonging to `C10`.

## Fix

Align the TFA columns to the expression matrix by name once, where `TFA_df_lovering` is built:

```r
TFA_df_lovering <- TFA_df[rownames(TFA_df) %in% Lovering_TF_list, colnames(variable_genes), drop=FALSE]
```

This makes the positional assignment correct and leaves the `rbind()` path unchanged. Runs
made before this fix need their TF regulator assignment (and TF2target enrichment) redone;
clustering, metabolite and protein regulators can be reused.
