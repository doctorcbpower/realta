import logging

from realta.imf.base import IMF
from realta.imf.chabrier import ChabrierIMF
from realta.imf.kroupa import KroupaIMF
from realta.imf.salpeter import SalpeterIMF

logger = logging.getLogger("realta")


def get_imf(imf_type: int, slope: float | None = None) -> IMF:
    """Construct the IMF selected by config.imf_type.

    1=Salpeter, 2=Kroupa (default), 3=Chabrier (see each IMF class's own
    docstring for its specific provenance). Falls back to SalpeterIMF
    with a logged warning for any other value, rather than raising --
    this is a permissive default, not itself a documented Power et al.
    (2009) behaviour.

    `slope` (config.imf_slope) overrides SalpeterIMF's own alpha
    (default 2.35) when imf_type=1 -- a continuously-sweepable single
    power-law slope, needed for Figure 4's (alpha_IMF, f_bin)
    degeneracy grid (docs/science/paper1-detailed-work-breakdown.md,
    item A4). Deliberately Salpeter-only: Kroupa's multi-segment break
    structure has no single slope to sweep in the same unambiguous
    way, and Chabrier's log-normal low-mass form has no power-law
    slope at all -- this is an additive extension to the existing
    Salpeter preset, not a new IMF family, per A4's own scope note.
    Ignored (with no warning) for imf_type != 1.
    """
    if imf_type == 1:
        return SalpeterIMF() if slope is None else SalpeterIMF(alpha=slope)
    elif imf_type == 2:
        return KroupaIMF()
    elif imf_type == 3:
        return ChabrierIMF()
    else:
        logger.warning(f"Unknown IMF type {imf_type}, defaulting to SalpeterIMF")
        return SalpeterIMF() if slope is None else SalpeterIMF(alpha=slope)
