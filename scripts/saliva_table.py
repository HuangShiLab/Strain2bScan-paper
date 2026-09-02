"""Read `results/saliva_strain_long.tsv` by column NAME.

Reading this table positionally is how the saliva community matrices came to be built on
`support` -- a raw marker count -- instead of an abundance. `support` scales with how many
markers a cluster happens to carry in the database and with sequencing depth, so two clusters
at the same true abundance get different values. This module exists so the column is an
explicit, checked choice rather than an index that silently means something else after a
format change.

Which column to use:

  sample_fraction   fraction of the sample's tag observations assigned to this cluster.
                    The only one comparable ACROSS species, so the one a community matrix
                    (Bray-Curtis, PERMANOVA, ordination) must be built on.
  global_abundance  share of the strain-resolved part of the sample, in cells.
  within_abund      sums to 1.0 WITHIN each species. Fine for "which strain of this species",
                    meaningless once several species are put in one matrix.
  support, coverage, depth   diagnostics, not abundances.
"""

COMMUNITY_ABUNDANCE = "sample_fraction"


def read_long(path, require=()):
    """Rows of `path` as dicts keyed by header name.

    `require` names columns the caller cannot work without; a missing one is a hard error,
    because the alternative -- falling back to whatever sits at that index -- is the bug this
    module was written to prevent. A table written before those columns existed must be
    regenerated with `profile_saliva.py`, which needs the reads.
    """
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, l.rstrip("\n").split("\t"))) for l in f if l.strip()]
    missing = [c for c in require if c not in header]
    if missing:
        raise SystemExit(
            f"{path}: missing column(s) {missing}; header is {header}.\n"
            f"This table predates them -- regenerate it with scripts/profile_saliva.py "
            f"(needs the saliva reads). Do NOT substitute another column: `support` is a "
            f"marker count and `within_abund` is normalised within a species, and neither "
            f"is a cross-species abundance."
        )
    return rows
