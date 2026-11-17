"""Options, errors and warnings of the Python port (the counterparts of opts, error() and warning() of the Matlab
code). The option names and defaults are those of pfc_checks('options') (tol, tol_div, write, solved, name), plus
``checks``, which exists only in Python."""

_DEFAULTS = {"tol": 1e-6, "tol_div": 1e-9, "write": True, "solved": False, "name": "", "checks": True}


class PfcError(ValueError):
    """An error of the port, with the identifier of the Matlab code (for example 'pfc:rpfc:negativeSelfTerm') in
    ``identifier``. Identifiers of the form 'pfc:input:...' that have no Matlab counterpart are listed in the README."""

    def __init__(self, identifier, message):
        super().__init__(message)
        self.identifier = identifier


class PfcWarning(UserWarning):
    """A warning of the port; ``identifier`` holds the Matlab identifier (for example 'pfc:apfc:cyclic')."""

    def __init__(self, identifier, message):
        super().__init__(message)
        self.identifier = identifier


class Options:
    """Options of apfc and rpfc, with the names and defaults of opts of apfc.m and rpfc.m:

    tol      deletion and classification threshold (default 1e-6, in the units of the case: MW, Mvar)
    tol_div  guard against division by (near) zero, negligible self-term threshold and attenuation slack (1e-9)
    write    accepted for compatibility; no effect: apfc and rpfc never write files (use write_csv)
    solved   accepted for compatibility; no effect: the input of the port is always a solved case
    name     accepted for compatibility; stored in info["opts"] (the case name is Case.name)
    checks   Python only (default True): False switches off the optional checks (attenuation range,
             cyclic-digraph warning, negative self-term guard); the guards that protect the solve
             (zero inflow, spectral radius at or above 1) always run.

    Unknown names raise PfcError (a ValueError) 'pfc:options:unknownField'; invalid values raise
    'pfc:options:value', as in pfc_checks('options')."""

    __slots__ = tuple(_DEFAULTS)

    def __init__(self, **kwargs):
        for k in kwargs:
            if k not in _DEFAULTS:
                raise PfcError("pfc:options:unknownField",
                               f"Unknown option '{k}'. Known options: {', '.join(_DEFAULTS)}.")
        for k, v in _DEFAULTS.items():
            object.__setattr__(self, k, kwargs.get(k, v))
        for k in ("tol", "tol_div"):
            v = getattr(self, k)
            ok = isinstance(v, (int, float)) and not isinstance(v, bool) and v > 0 and v != float("inf") and v == v
            if not ok:
                raise PfcError("pfc:options:value", f"opts.{k} must be a positive finite scalar.")
            object.__setattr__(self, k, float(v))
        if not isinstance(self.name, str):
            raise PfcError("pfc:options:value", "opts.name must be a character vector.")
        for k in ("write", "solved", "checks"):
            object.__setattr__(self, k, bool(getattr(self, k)))

    def __setattr__(self, key, value):
        raise AttributeError("Options is immutable; use replace()")

    def replace(self, **kwargs):
        d = self.as_dict(); d.update(kwargs)
        return Options(**d)

    def as_dict(self):
        return {k: getattr(self, k) for k in _DEFAULTS}

    def __eq__(self, other):
        return isinstance(other, Options) and self.as_dict() == other.as_dict()

    def __repr__(self):
        return "Options(" + ", ".join(f"{k}={getattr(self, k)!r}" for k in _DEFAULTS) + ")"


def resolve_options(options=None, tol=None, tol_div=None):
    """The options of one call: ``options`` (an Options, a dict, or None for the defaults), with the keyword
    arguments tol and tol_div of T13a and T13b taking precedence when given."""
    if options is None:
        options = Options()
    elif isinstance(options, dict):
        options = Options(**options)
    elif not isinstance(options, Options):
        raise PfcError("pfc:options:type", "opts must be a struct.")
    over = {k: v for k, v in (("tol", tol), ("tol_div", tol_div)) if v is not None}
    return options.replace(**over) if over else options
