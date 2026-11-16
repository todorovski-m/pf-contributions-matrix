"""Row and column labels, as in apfc.m (and rpfc.m)."""


def disambiguate_parallel(from_labels, to_labels):
    """Parallel-branch scheme of apfc.m: the second and later branches whose from-end label repeats get the
    suffix _2, _3, ... on the from-end label and the same suffix on the to-end label."""
    from_labels, to_labels = list(from_labels), list(to_labels)
    seen = {}
    for i, nm in enumerate(from_labels):
        if nm in seen:
            seen[nm] += 1
            from_labels[i] = f"{nm}_{seen[nm]}"
            to_labels[i] = f"{to_labels[i]}_{seen[nm]}"
        else:
            seen[nm] = 1
    return from_labels, to_labels


def unique_labels(labels):
    """unique_labels of apfc.m: the second and later occurrences of a label get the suffix _2, _3, ...;
    a suffix is skipped if the resulting label is already in use."""
    labels = list(labels)
    used = set(labels)
    seen = {}
    for i, nm in enumerate(labels):
        if nm in seen:
            k = seen[nm] + 1
            while f"{nm}_{k}" in used:
                k += 1
            seen[nm] = k
            labels[i] = f"{nm}_{k}"
            used.add(labels[i])
        else:
            seen[nm] = 1
    return labels
