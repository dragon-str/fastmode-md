"""Unit tests for the pure functions in fastmode.py.

These tests do not need GROMACS.  They cover the parsers, the mass
repartition transform, the mdp builder, the log readers, and the four-check
validator.
"""

import numpy as np
import pytest

import fastmode


# --------------------------------------------------------------- tiny systems

TOP_ONE_MOL = """\
[ moleculetype ]
; name  nrexcl
MOL     3

[ atoms ]
; nr type resnr residue atom cgnr charge mass
1  C  1  MOL  C  1  0.0  12.0000
2  H  1  MOL  H  1  0.0   1.0080

[ bonds ]
1 2
"""

TOP_TWO_H = """\
[ moleculetype ]
MOL     3

[ atoms ]
; nr type resnr residue atom cgnr charge mass
1  C  1  MOL  C  1  0.0  12.0000
2  H  1  MOL  H1 1  0.0   1.0080
3  H  1  MOL  H2 1  0.0   1.0080

[ bonds ]
1 2
1 3
"""

TOP_WATER = """\
[ moleculetype ]
SOL     2

[ atoms ]
; nr type resnr residue atom cgnr charge mass
1  OW 1  SOL  OW 1  0.0  15.9994
2  HW 1  SOL  HW 1  0.0   1.0080

[ bonds ]
1 2
"""


# ---------------------------------------------------------------------- model

def test_read_gro_natoms(tmp_path):
    p = tmp_path / "box.gro"
    p.write_text("comment line\n  123\n1MOL ...\n")
    assert fastmode.read_gro_natoms(str(p)) == 123


def test_parse_topology_counts_atoms_and_bonds():
    mols = fastmode.parse_topology(TOP_ONE_MOL)
    assert list(mols) == ["MOL"]
    mol = mols["MOL"]
    assert mol["order"] == [1, 2]
    assert mol["atoms"][1]["name"] == "C"
    assert mol["atoms"][2]["name"] == "H"
    assert mol["atoms"][2]["mass"] == pytest.approx(1.0080)
    assert mol["bonds"] == [(1, 2)]


def test_is_water():
    assert fastmode.is_water("SOL")
    assert fastmode.is_water("tip3")
    assert not fastmode.is_water("MOL")
    assert not fastmode.is_water("ALA")


def test_bonded_heavy():
    mols = fastmode.parse_topology(TOP_ONE_MOL)
    atoms = mols["MOL"]["atoms"]
    bonds = mols["MOL"]["bonds"]
    assert fastmode.bonded_heavy(atoms, bonds, 2) == 1
    assert fastmode.bonded_heavy(atoms, bonds, 1) is None


def test_is_protein(tmp_path):
    protein = tmp_path / "protein.top"
    protein.write_text(TOP_ONE_MOL)
    assert fastmode.is_protein(str(protein))

    water = tmp_path / "water.top"
    water.write_text(TOP_WATER)
    assert not fastmode.is_protein(str(water))


# ----------------------------------------------------------- hmr mass transform

def test_hmr_transform_conserves_mass(tmp_path):
    out = tmp_path / "hmr.top"
    path, stats = fastmode.hmr_transform(TOP_ONE_MOL, str(out))
    assert path == str(out)
    assert stats["hydrogens"] == 1
    assert stats["rel_mass_change"] <= fastmode.MASS_TOL

    mols = fastmode.parse_topology(out.read_text())
    atoms = mols["MOL"]["atoms"]
    assert atoms[2]["mass"] == pytest.approx(1.0080 * fastmode.HMR_FACTOR)
    assert atoms[1]["mass"] == pytest.approx(12.0 - 1.0080 * (fastmode.HMR_FACTOR - 1.0))
    total = sum(a["mass"] for a in atoms.values())
    assert total == pytest.approx(12.0 + 1.0080, rel=1e-6)


def test_hmr_transform_debits_multiple_hydrogens(tmp_path):
    out = tmp_path / "hmr.top"
    _, stats = fastmode.hmr_transform(TOP_TWO_H, str(out))
    assert stats["hydrogens"] == 2
    mols = fastmode.parse_topology(out.read_text())
    carbon = mols["MOL"]["atoms"][1]["mass"]
    assert carbon == pytest.approx(12.0 - 2 * 1.0080 * (fastmode.HMR_FACTOR - 1.0))
    total = sum(a["mass"] for a in mols["MOL"]["atoms"].values())
    assert total == pytest.approx(12.0 + 2 * 1.0080, rel=1e-6)


def test_hmr_transform_skips_water(tmp_path):
    out = tmp_path / "hmr.top"
    path, stats = fastmode.hmr_transform(TOP_WATER, str(out))
    assert path is None
    assert stats["hydrogens"] == 0


def test_hmr_transform_scales_solute(tmp_path):
    out = tmp_path / "hmr.top"
    _, stats = fastmode.hmr_transform(TOP_ONE_MOL, str(out), scale=1.1)
    assert stats["rel_mass_change"] <= fastmode.MASS_TOL
    mols = fastmode.parse_topology(out.read_text())
    total = sum(a["mass"] for a in mols["MOL"]["atoms"].values())
    assert total == pytest.approx((12.0 + 1.0080) * 1.1, rel=1e-5)


# --------------------------------------------------------------------- mdp I/O

def test_output_interval_rounds_up():
    assert fastmode.output_interval(1000, 1, 10) == 10
    assert fastmode.output_interval(1000, 1, 0) == 1
    assert fastmode.output_interval(1000, 5, 7) == 10
    assert fastmode.output_interval(1000, 5, 1) == 5


def test_build_mdp_plain():
    spec = {"dt": 0.005, "hmr": False, "mts": 0, "nstlist": 10,
            "tol": 0.005, "constraints": "all-bonds"}
    text, nstxout = fastmode.build_mdp(spec, 40000, 200.0)
    assert "constraints     = all-bonds" in text
    assert "mts = yes" not in text
    assert nstxout == 10


def test_build_mdp_mts_is_multiple_of_factor():
    spec = {"dt": 0.005, "hmr": False, "mts": 3, "nstlist": 10, "tol": 0.005}
    text, nstxout = fastmode.build_mdp(spec, 40000, 200.0)
    assert "mts = yes" in text
    assert "mts-level2-factor = 3" in text
    assert nstxout % 3 == 0


def test_setting_label():
    spec = {"dt": 0.005, "hmr": True, "mts": 0, "nstlist": 10, "tol": 0.005,
            "constraints": "all-bonds"}
    assert fastmode.setting_label(spec) == "dt5fs_hmron_mtsoff_nstlist10_tol0.005_allbonds"


def test_vary_does_not_mutate():
    base = {"dt": 0.002, "hmr": False}
    out = fastmode.vary(base, dt=0.005)
    assert out["dt"] == 0.005
    assert base["dt"] == 0.002


def test_grid_specs_size():
    base = dict(fastmode.REF_SPEC)
    assert len(list(fastmode.grid_specs(base, protein=True))) == 4 * 2 * 3 * 4 * 2
    assert len(list(fastmode.grid_specs(base, protein=False))) == 4 * 1 * 3 * 4 * 2


# ------------------------------------------------------------------ log readers

def test_performance_ns_per_day(tmp_path):
    p = tmp_path / "run.log"
    p.write_text("step 0\n\n   Performance:    123.4     ns/day\n")
    assert fastmode.performance_ns_per_day(str(p)) == pytest.approx(123.4)


def test_performance_ns_per_day_missing(tmp_path):
    p = tmp_path / "run.log"
    p.write_text("no performance here\n")
    assert fastmode.performance_ns_per_day(str(p)) is None


def test_conserved_drift(tmp_path):
    p = tmp_path / "run.log"
    p.write_text("Conserved energy drift: 0.0123 kJ/mol/ps per atom\n")
    assert fastmode.conserved_drift(str(p)) == pytest.approx(0.0123)


def test_conserved_drift_missing(tmp_path):
    p = tmp_path / "run.log"
    p.write_text("nothing\n")
    assert fastmode.conserved_drift(str(p)) is None


def test_read_xvg(tmp_path):
    p = tmp_path / "e.xvg"
    p.write_text(
        "# comment\n"
        '@ title "energy"\n'
        '@ s0 legend "Temperature"\n'
        '@ s1 legend "Density"\n'
        "0.0 300.0 1000.0\n"
        "1.0 301.0 1001.0\n"
    )
    cols, legends = fastmode.read_xvg(str(p))
    assert legends == ["Temperature", "Density"]
    assert np.allclose(cols[0], [0.0, 1.0])
    assert np.allclose(cols[1], [300.0, 301.0])
    assert np.allclose(cols[2], [1000.0, 1001.0])


def test_last_error_fatal():
    stderr = "reading mdp\nFatal error:\nThe bond is too long\n----\n"
    assert fastmode.last_error(stderr) == "The bond is too long"


def test_last_error_fallback():
    assert fastmode.last_error("some plain problem") == "some plain problem"


# ------------------------------------------------------------------- validator

def _rec(**kw):
    rec = {"drift": 0.0, "temperature": 300.0, "density": 1.0,
           "rdf_g": [1.0, 2.0], "reasons": []}
    rec.update(kw)
    return rec


def _ref():
    return {"temperature": 300.0, "density": 1.0, "rdf_g": [1.0, 2.0]}


def test_validate_passes():
    rec = _rec()
    assert fastmode.validate(rec, _ref())
    assert rec["pass"] is True
    assert rec["reasons"] == []


def test_validate_flags_drift():
    rec = _rec(drift=0.05)
    assert not fastmode.validate(rec, _ref())
    assert any("drift" in r for r in rec["reasons"])


def test_validate_flags_temperature():
    rec = _rec(temperature=304.0)
    assert not fastmode.validate(rec, _ref())
    assert any("temperature" in r for r in rec["reasons"])


def test_validate_flags_density():
    rec = _rec(density=1.02)
    assert not fastmode.validate(rec, _ref())
    assert any("density" in r for r in rec["reasons"])


def test_validate_flags_missing_rdf():
    rec = _rec()
    del rec["rdf_g"]
    assert not fastmode.validate(rec, _ref())
    assert any("rdf" in r for r in rec["reasons"])


def test_validate_rejected_is_failure():
    rec = _rec(rejected=True)
    assert not fastmode.validate(rec, _ref())
    assert rec["pass"] is False
