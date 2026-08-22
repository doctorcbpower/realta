import logging

from realta.imf.base import IMF
from realta.imf.chabrier import ChabrierIMF
from realta.imf.kroupa import KroupaIMF
from realta.imf.salpeter import SalpeterIMF

logger = logging.getLogger("realta")


def get_imf(imf_type: int) -> IMF:
    """Factory function to get IMF instance."""
    if imf_type == 1:
        return SalpeterIMF()
    elif imf_type == 2:
        return KroupaIMF()
    elif imf_type == 3:
        return ChabrierIMF()
    else:
        logger.warning(f"Unknown IMF type {imf_type}, defaulting to SalpeterIMF")
        return SalpeterIMF()
