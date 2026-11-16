"""Writing result tables in the layout of the CSV files of apfc.m / rpfc.m."""


def write_csv(table, path):
    """Write ``table`` (index = row labels) as CSV: header ``Line,<column labels>``, one line per row with
    the label followed by the values with 17 significant digits (``%.17g``), lines ending in ``\\n``.
    This is the layout of the files written under Octave with tests/octave_shims/writetable.m."""
    with open(path, "w", newline="\n") as fh:
        fh.write(",".join(["Line"] + [str(c) for c in table.columns]) + "\n")
        for label, row in zip(table.index, table.to_numpy(float)):
            fh.write(",".join([str(label)] + ["%.17g" % v for v in row]) + "\n")
