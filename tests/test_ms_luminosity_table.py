"""Unit tests for MSLuminosityTable.get_lbol()'s mass-rescaling formula
-- closes the gap noted in docs/provenance.md Section 4: this was the
site of this session's own mass-normalization bug (the MS luminosity
curve was baked to a fixed fiducial 1e6 Msun cluster and added directly
to lumx_tot from a population of a different actual mass), fixed but
previously only verified manually via notebook re-execution, not by any
automated test.
"""

import pytest

from realta.io.tables import MSLuminosityTable


def test_get_lbol_scales_linearly_with_total_mass():
    """L_bol(age, mass) must scale linearly in mass at fixed age.

    SSP bolometric luminosity scales linearly with total mass formed at
    fixed IMF/metallicity/age -- see MSLuminosityTable's class
    docstring. Checked at several ages across the tabulated range.
    """
    table = MSLuminosityTable(imetal=2)
    assert table.loaded

    for age_myr in (0.5, 5.0, 20.0, 80.0):
        lbol_1x = table.get_lbol(age_myr, table.FIDUCIAL_CLUSTER_MASS_MSUN)
        lbol_2x = table.get_lbol(age_myr, 2.0 * table.FIDUCIAL_CLUSTER_MASS_MSUN)
        lbol_half = table.get_lbol(age_myr, 0.5 * table.FIDUCIAL_CLUSTER_MASS_MSUN)

        assert lbol_1x > 0.0
        assert lbol_2x == pytest.approx(2.0 * lbol_1x, rel=1e-9)
        assert lbol_half == pytest.approx(0.5 * lbol_1x, rel=1e-9)


def test_get_lbol_zero_mass_gives_zero_luminosity():
    table = MSLuminosityTable(imetal=2)
    assert table.get_lbol(10.0, 0.0) == 0.0


def test_get_lbol_outside_tabulated_age_range_returns_zero():
    """No extrapolation outside the tabulated 0.1-100 Myr range (by design
    -- see the class docstring's "Domain of validity" note), regardless
    of the mass passed in.
    """
    table = MSLuminosityTable(imetal=2)
    mass = table.FIDUCIAL_CLUSTER_MASS_MSUN

    assert table.get_lbol(0.05, mass) == 0.0  # below the tabulated minimum
    assert table.get_lbol(200.0, mass) == 0.0  # above the tabulated maximum
    assert table.get_lbol(0.0, mass) == 0.0
    assert table.get_lbol(-5.0, mass) == 0.0


def test_get_lbol_matches_manually_computed_value_at_fiducial_mass():
    """Cross-check against a value captured from an actual run, at the
    fiducial 1e6 Msun mass (no rescaling applied -- ratio is exactly 1).
    """
    table = MSLuminosityTable(imetal=2)
    lbol = table.get_lbol(10.0, table.FIDUCIAL_CLUSTER_MASS_MSUN)
    assert lbol == pytest.approx(7.971425654408602e41, rel=1e-9)


def test_all_three_metallicity_tables_load_and_scale():
    for imetal in (1, 2, 3):
        table = MSLuminosityTable(imetal=imetal)
        assert table.loaded, f"imetal={imetal} table failed to load"
        lbol_1x = table.get_lbol(10.0, table.FIDUCIAL_CLUSTER_MASS_MSUN)
        lbol_2x = table.get_lbol(10.0, 2.0 * table.FIDUCIAL_CLUSTER_MASS_MSUN)
        assert lbol_1x > 0.0
        assert lbol_2x == pytest.approx(2.0 * lbol_1x, rel=1e-9)


def test_missing_data_file_returns_zero_not_crash():
    table = MSLuminosityTable(imetal=2, data_dir="/nonexistent/path")
    assert not table.loaded
    assert table.get_lbol(10.0, table.FIDUCIAL_CLUSTER_MASS_MSUN) == 0.0
