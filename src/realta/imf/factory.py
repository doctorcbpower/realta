import logging

from realta.imf.base import IMF
from realta.imf.chabrier import ChabrierIMF
from realta.imf.kroupa import KroupaIMF
from realta.imf.salpeter import SalpeterIMF

logger = logging.getLogger("realta")


def get_imf(imf_type: int) -> IMF:
    """Construct the IMF selected by config.imf_type.

    1=Salpeter, 2=Kroupa (default), 3=Chabrier (see each IMF class's own
    docstring for its specific provenance). Falls back to SalpeterIMF
    with a logged warning for any other value, rather than raising --
    this is a permissive default, not itself a documented Power et al.
    (2009) behaviour.
    """
    if imf_type == 1:
        return SalpeterIMF()
    elif imf_type == 2:
        return KroupaIMF()
    elif imf_type == 3:
        return ChabrierIMF()
    else:
        logger.warning(f"Unknown IMF type {imf_type}, defaulting to SalpeterIMF")
        return SalpeterIMF()
