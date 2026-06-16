# Zero-Offset Schakel 2011 RT-SE Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a new reactive-transport-driven seismoelectric simulator that keeps the geometry and output contract of `seismoelectric_offset_liu2018_spectral.py`, forces zero finite offset, and replaces the waveform core with a Schakel 2011 Sommerfeld-style solver.

**Architecture:** The new script imports the validated reactive-transport mapping, Schakel 2010 coefficient solver, and plotting/output helpers from `seismoelectric_offset_liu2018_spectral.py`. It adds a Schakel 2011 pressure-source/Sommerfeld waveform synthesizer and a matching peak-amplitude routine, then writes the same result file types under a new output directory.

**Tech Stack:** Python, NumPy, Pandas, Matplotlib, SciPy if available, `unittest`.

---

### Task 1: Tests For The New Output Contract And Solver Hooks

**Files:**
- Create: `tests/test_zerooffset_schakel2011_rt.py`
- Read: `seismoelectric_offset_liu2018_spectral.py`
- Read: `schakel2011_geophysics_reproduction.py`

- [ ] **Step 1: Write failing tests**

Create tests that import `seismoelectric_zerooffset_schakel2011_sommerfeld.py` and verify:

```python
def test_default_config_forces_zero_offset_and_preserves_receiver_geometry():
    cfg = model.ZeroOffsetSchakelConfig()
    base = base_model.SEConfig()
    assert cfg.offset_D == 0.0
    assert cfg.z_s == base.z_s
    assert cfg.receiver_z_min == base.receiver_z_min
    assert cfg.receiver_z_max == base.receiver_z_max
    assert cfg.receiver_spacing == base.receiver_spacing
```

```python
def test_sommerfeld_pressure_coefficient_uses_schakel_pressure_normalization():
    coeff = {"R_E": 6.0 + 8.0j}
    actual = model.pressure_normalized_re_from_coeff(coeff, rho_f=1000.0, omega=2.0)
    assert actual == (6.0 + 8.0j) / (1000.0 * 4.0)
```

```python
def test_waveform_synthesis_returns_existing_output_shape():
    cfg = model.ZeroOffsetSchakelConfig()
    cfg.receiver_z_min = -0.002
    cfg.receiver_z_max = 0.002
    cfg.receiver_spacing = 0.002
    cfg.waveform_nt = 24
    row = pd.Series({"Porosity": 0.24, "Permeability_mD": 100.0, "Tortuosity": 2.0, "OutletHConc": 1.0e-10})
    z, t, u = model.synthesize_waveforms_schakel2011(row, cfg, n_frequencies=3, n_theta=6, integration_method="fixed")
    assert u.shape == (len(z), len(t))
    assert len(z) == 3
    assert len(t) == cfg.waveform_nt
```

```python
def test_pipeline_writes_same_named_result_types(tmp_path):
    input_csv = tmp_path / "rt.csv"
    input_csv.write_text("Time_s,Porosity,Permeability_mD,Tortuosity,OutletHConc\n0,0.24,100,2,1e-10\n60,0.25,120,2.1,1e-10\n")
    outdir = tmp_path / "out"
    model.run_simulation(input_csv, outdir, n_frequencies=3, n_theta=6, peak_n_frequencies=3, peak_n_theta=6, integration_method="fixed")
    for name in [...]:
        assert (outdir / name).exists()
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
python -m unittest tests.test_zerooffset_schakel2011_rt -v
```

Expected: import failure because the new script does not exist yet.

### Task 2: New Schakel 2011 Zero-Offset Simulator

**Files:**
- Create: `seismoelectric_zerooffset_schakel2011_sommerfeld.py`
- Modify: none

- [ ] **Step 1: Implement config and helpers**

Create `ZeroOffsetSchakelConfig(SEConfig)` with `offset_D = 0.0`, Schakel source parameters (`transducer_radius_m`, frequency band, evanescent gamma limit), and helper functions mirroring `schakel2011_geophysics_reproduction.py`: Bessel fallbacks, piston directivity, source FFT spectrum, EM vertical wavenumber, pressure-normalized `R_E`.

- [ ] **Step 2: Implement row-dependent Sommerfeld response**

For each reactive-transport row and frequency, evaluate the Schakel 2011 real-angle and evanescent branches. Compute row-dependent `R_E` with the existing Schakel 2010 boundary solver and pressure-normalize it as `R_E/(rho_f omega^2)`. Use `source_z_m = -cfg.z_s`, `radial_m = abs(cfg.offset_D) = 0`, and the existing receiver `z` positions.

- [ ] **Step 3: Implement waveform and peak routines**

Use positive-frequency synthesis with `exp(i omega t)` and the same `T0`, `waveform_t_before`, `waveform_t_after`, `waveform_nt` window as the existing script. Save the new snapshot as `waveform_snapshot_schakel2011.*` and expose peak columns named with the existing spectral names plus Schakel-specific aliases for compatibility.

- [ ] **Step 4: Implement CLI and run summary**

Provide `--input`, `--outdir`, geometry overrides, frequency/integration controls, and `run_simulation(...)`. Default `--outdir` is `se_results_zerooffset_schakel2011`.

### Task 3: Verification, Results, And Audit

**Files:**
- Create: `se_results_zerooffset_schakel2011/`
- Create: `se_results_zerooffset_schakel2011/formula_audit.md`

- [ ] **Step 1: Run focused tests**

Run:

```bash
python -m unittest tests.test_zerooffset_schakel2011_rt -v
```

- [ ] **Step 2: Run related regression tests**

Run:

```bash
python -m unittest tests.test_waveform_diagnostics tests.test_schakel2011_geophysics_reproduction -v
```

- [ ] **Step 3: Generate results**

Run:

```bash
python seismoelectric_zerooffset_schakel2011_sommerfeld.py --input global_evolution.xlsx --outdir se_results_zerooffset_schakel2011
```

- [ ] **Step 4: Dispatch two subagent audits**

Ask one reviewer to audit formula/theory alignment and one reviewer to audit geometry/output compliance against the user requirements. Fix any material issues.
