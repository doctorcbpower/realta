"""Unit tests for UVLuminosityTable -- see docs/provenance.md Section 7.

The real fuv_lbol_z*.dat data (generated via
scripts/generate_fuv_luminosities.py, which requires FSPS + SPS_HOME)
now exists in src/realta/data/, so these tests cover: (1) that the real
files load and scale correctly (analogous to
test_ms_luminosity_table.py's test_all_three_metallicity_tables_load_and_scale),
(2) the same graceful missing-data degradation MSLuminosityTable has,
for an explicitly nonexistent data_dir, and (3) the rescaling
arithmetic itself, exercised by writing a small synthetic table to a
temp data_dir rather than depending on the real FSPS values -- this is
the part of the class that is independent of what the real tabulated
values turn out to be.
"""

import pytest

from realta.io.tables import UVLuminosityTable


def test_all_three_metallicity_tables_load_and_scale():
    for imetal in (1, 2, 3):
        table = UVLuminosityTable(imetal=imetal)
        assert table.loaded, f"imetal={imetal} table failed to load"
        luv_1x = table.get_luv(10.0, table.FIDUCIAL_CLUSTER_MASS_MSUN)
        luv_2x = table.get_luv(10.0, 2.0 * table.FIDUCIAL_CLUSTER_MASS_MSUN)
        assert luv_1x > 0.0
        assert luv_2x == pytest.approx(2.0 * luv_1x, rel=1e-9)


def test_missing_data_file_explicit_nonexistent_dir():
    table = UVLuminosityTable(imetal=2, data_dir="/nonexistent/path")
    assert not table.loaded
    assert table.get_luv(10.0, table.FIDUCIAL_CLUSTER_MASS_MSUN) == 0.0


def _write_synthetic_table(tmp_path, filename: str):
    # Two points is the minimum get_luv() needs to interpolate; values
    # are arbitrary but well within the 0.1-100 Myr domain.
    path = tmp_path / filename
    path.write_text(
        "# synthetic test table\n"
        "# Cluster mass: 1e6 Msun\n"
        "# age_myr  log10_lfuv_total_erg_s\n"
        " 1.0  42.0\n"
        "10.0  43.0\n"
    )
    return path


def test_get_luv_scales_linearly_with_total_mass(tmp_path):
    _write_synthetic_table(tmp_path, "fuv_lbol_z8e-3.dat")
    table = UVLuminosityTable(imetal=2, data_dir=str(tmp_path))
    assert table.loaded

    luv_1x = table.get_luv(5.0, table.FIDUCIAL_CLUSTER_MASS_MSUN)
    luv_2x = table.get_luv(5.0, 2.0 * table.FIDUCIAL_CLUSTER_MASS_MSUN)
    luv_half = table.get_luv(5.0, 0.5 * table.FIDUCIAL_CLUSTER_MASS_MSUN)

    assert luv_1x > 0.0
    assert luv_2x == pytest.approx(2.0 * luv_1x, rel=1e-9)
    assert luv_half == pytest.approx(0.5 * luv_1x, rel=1e-9)


def test_get_luv_zero_mass_gives_zero_luminosity(tmp_path):
    _write_synthetic_table(tmp_path, "fuv_lbol_z8e-3.dat")
    table = UVLuminosityTable(imetal=2, data_dir=str(tmp_path))
    assert table.get_luv(5.0, 0.0) == 0.0


def test_get_luv_outside_tabulated_age_range_returns_zero(tmp_path):
    _write_synthetic_table(tmp_path, "fuv_lbol_z8e-3.dat")
    table = UVLuminosityTable(imetal=2, data_dir=str(tmp_path))
    mass = table.FIDUCIAL_CLUSTER_MASS_MSUN

    assert table.get_luv(0.5, mass) == 0.0  # below tabulated minimum (1.0)
    assert table.get_luv(20.0, mass) == 0.0  # above tabulated maximum (10.0)
    assert table.get_luv(0.0, mass) == 0.0
    assert table.get_luv(-5.0, mass) == 0.0
