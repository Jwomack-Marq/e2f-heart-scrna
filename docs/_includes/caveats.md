::: {.callout-important}
## The five things that constrain every result in this book

**1. n = 1 animal per genotype × timepoint.** The two lanes per sample are the same
library sequenced twice, not biological replicates. Cell-level Wilcoxon p-values
therefore treat technical replicates as biological ones (pseudoreplication): they
under-estimate variance and return extreme values, many so small they underflow to
zero in double precision and are floored at `1e-300`. **Use the p-axis to rank
candidates, never to claim significance.** Tables in the app sort by effect size for
this reason.

**2. The KO and WT animals are different sexes.** Y-linked genes top the KO-up list.
Seven genes — `Eif2s3y`, `Kdm5d`, `Uty`, `Ddx3y`, `Xist`, `Tsix`, `Gt(ROSA)26Sor` —
are flagged as `sex/construct` in every table and excluded from every enrichment
input. They are left *in* the DE tables, visible, because deleting a row is how a
confound becomes invisible.

**3. The knockout is not confirmed at the transcript level.** `E2f7` and `E2f8` mRNA
are not reduced in the KO cells — most likely a conditional allele that a 3'-biased
assay cannot see. Everything called a "KO effect" here is an effect of *this animal*
versus *that animal*.

**4. P7 was FACS cycling-enriched and P0 essentially was not.** The project notes put
the enrichment at 4.5–5.2×, so the raw cycling fraction *rises* from P0 to P7 in this
data while it *falls* in reality. Any P0-vs-P7 contrast on all cells partly reads out
the sort rather than development, which is why every temporal contrast exists in a
phase-matched (G1-only) stratum as well as a raw one. Only within-timepoint
comparisons (P7 KO vs WT) are free of it. Measured *inside* the cardiomyocyte
compartment of this bundle the ratio is milder than the headline figure — see
[Confounds and reproducibility](11-confounds-repro.qmd).

**5. Mitochondrially-encoded genes shift as a block between the libraries.** `mt-`
genes are up in the KO in all seven tested subclusters and down in none. Seven
independently clustered populations do not agree that cleanly by biology; this is a
read-fraction difference. It matters because those genes *are* the OXPHOS /
electron-transport terms that otherwise headline the KO-up enrichment.

Taken together: **this is a descriptive pilot.** It is well suited to generating
hypotheses and ranking candidates, and not suited to any claim of the form "the
knockout causes X". Valid inference needs a replicated, sex-matched cohort of at least
three animals per condition.
:::
