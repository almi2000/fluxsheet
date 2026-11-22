# Model verification and validation

This page records the scope of the current checks. A balance check verifies the
implemented equations at a solved point. It does not establish agreement with
experimental data or validate an entire operating range.

| Family | Existing verification | Experimental coverage and remaining work |
| --- | --- | --- |
| Thermodynamic packages and flash | Pure-component phase boundaries, mixture bubble and dew points, binary diagram endpoints, water/ethanol azeotrope and water/butanol phase split in `test_phase_physics.py` and `test_thermo_models.py`. | The azeotrope and liquid-liquid tests use selected literature-derived values. Offline comparisons with three selected ThermoML references are recorded below; broader independent series remain open. |
| Distillation, absorption and extraction | Fenske and Underwood hand calculations, a single-stage Raoult/Rachford-Rice comparison, total and per-stage balances, feed and side streams, pump-around, reactive stages, Murphree efficiency and phase equilibria in `test_stage_validation.py`, `test_columns.py`, `test_column_side_streams.py` and `test_more_units.py`. | Equilibrium reactive stages are compared with a measured methyl acetate column profile below; the comparison fails the screening limits. Independent reactive kinetics remain open. |
| Reactors | First-order CSTR and PFR closed forms, element balances, adiabatic enthalpy and independent integration of illustrative kinetic laws in `test_reactors.py` and `test_validation.py`. | The HDA kinetic constants in `test_validation.py` are illustrative, so the comparison validates the implementation rather than a measured reaction. The Gibbs reactor on fugacities is compared with measured ammonia equilibria below. Experimental kinetics remain open. |
| Dynamic simulation | Analytical first-order lag and mixed-tank response, integrated tank material and energy inventories, reactor comparison with a new steady state, integrated reactive CSTR mass and energy and controller tests in `test_dynamics.py`; a pure-propane boil/cool transient with time-grid refinement in `test_dynamic_phase_boundaries.py`. | The dynamic column is verified against the steady column below. Experimental transient data remain open. |
| Heat transfer and pressure changes | Exchanger balances, constant-heat-capacity effectiveness comparisons for co- and countercurrent flow, pipe grid convergence to an exponential thermal profile, friction and elevation, turbine isentropic work and multistage compressor tests in `test_core.py`, `test_utilities.py` and `test_more_units.py`; pure-propane evaporation and condensation against independent PH quadrature in both flow arrangements in `test_exchanger_phase_boundaries.py`; two-phase pipe grid convergence in `test_pipe_boundaries_audit.py`; Fanno flow of a gas (outlet Mach number, pressure and the choking length), the multistage compressor's discharge pressure with cooler pressure drops, heats of vaporization against the compound correlations at 1, 5 and 10 bar and liquid densities of the equations of state in `test_engineering_audit.py`. | Multicomponent phase-change heat transfer, two-phase choking against measurements and experimental equipment data remain open. |

## Automatic balance diagnostics

A successful solve exposes `balance_checks` as a list with `scope`, `quantity`,
`component`, `residual`, `unit`, `tolerance` and `status`. The status is
`ok`, `warning` or `not_checkable`. Residuals use inlet plus external transfer
minus outlet. Tolerances are 1e-5 of the larger sum of absolute terms, plus
1e-10 kg/s or mol/s for mass and material, and 0.1 W for energy. Warnings are
also sent through the normal warning list.

Reactor element and energy checks use atomic compositions and standard formation
enthalpies. Yield reactors are checked for mass and their energy when formation
data exist, while element checks are marked not checkable. Logic blocks act as
boundaries. The overall energy check includes the reported enthalpy changes at
property-method crossings into physical units. A crossing into a product or
virtual boundary lies outside that balance. Unsupported unit energy models are
marked not checkable.

Batch reactors are checked on their cycle-average input and output streams,
including the vent. For a constant-volume batch, the energy check includes the
change in pressure-volume energy over the cycle. If the vessel's initial
pressure differs from the inlet pressure, the energy balance is marked not
checkable because the filling work is not provided by the batch model.

The pytest suite records every successful `Solver.solve` return and rejects
unexpected balance warnings or uncheckable diagnostics. Deliberately injected
errors and known model limitations are declared by exact scope, quantity,
component and status; a declared exception must actually occur. This test guard
does not change the production solver's success or warning behavior.

## Reference-data campaign

The NIST ThermoML archive provides experimental measurement records. Raw archive
files are not redistributed in this repository, consistent with `data-sources.md`.
The registered cases can be run together with `python scripts/validate_models.py
--manifest docs/validation-campaign.json --offline --cache-dir PATH --output-dir PATH`.
Use `--download` instead of `--offline` to fetch missing records from NIST. The
cache must be outside the repository. JSON and Markdown outputs distinguish
screening failure (exit 1) from missing or invalid data (exit 2). The old
individual validation scripts remain available. The campaign manifest records
the limits used by those validators; mismatched registered limits are rejected.
A candidate reference for benzene/toluene is Lee et al., Fluid Phase Equilibria
280 (2009), DOI [10.1016/j.fluid.2009.03.012](https://doi.org/10.1016/j.fluid.2009.03.012),
at 333.15 K and 353.15 K. The ideal method uses no fitted binary parameter.
Whether the underlying pure-component vapor-pressure correlations include this
article has not been established.

The XML reports point-specific expanded pressure uncertainties of 0.2 -0.5 kPa
at 95% confidence for the two series. Before running the comparison, the
acceptance limits were set to 5% mean absolute relative pressure error and 10%
maximum error for each series. These engineering screening limits are broader
than the reported measurement uncertainty. The raw file is identified by SHA-256
`c343b21482c9693b87a00823b5af2ad5e73eb724d5254f6133ca4c58aaef0dd8`.

| Temperature | Points | Mean relative pressure error | Maximum error | Maximum-error point |
| --- | --- | --- | --- | --- |
| 333.15 K | 11 | 1.071% | 2.215% | x_benzene = 0.202, 25.349 kPa predicted, 24.8 kPa measured |
| 353.15 K | 11 | 0.749% | 2.479% | x_benzene = 0.207, 51.752 kPa predicted, 50.5 kPa measured |

Both series pass the stated screening limits. Reproduce the comparison with
`python scripts/validate_btx_reference.py PATH_TO_XML`. The XML contains
pressure data for these two series but no vapor-composition data in the
current parser, so this result assesses bubble pressure only.

## Ethanol/water vapor-liquid equilibrium

Lai et al., DOI
[10.1016/j.jct.2013.08.020](https://doi.org/10.1016/j.jct.2013.08.020),
report 16 T-x-y points at 101.3 kPa. These data appear in this project's NRTL and
UNIQUAC fit sources. The comparison below uses UNIFAC Dortmund group parameters,
without the project's fitted binary pair. The provenance of all underlying
pure-component and group parameters has not been independently audited.

Before comparison, screening limits were set at mean/max 2/5 K for temperature
and 0.03/0.08 for ethanol vapor mole fraction. The XML reports point-specific
expanded uncertainties at 95% confidence of 0.13-0.45 K and 0.005-0.015 mole
fraction. The screening limits are broader than those measurement uncertainties.

| Quantity | Mean absolute error | Maximum error | Maximum-error point |
| --- | --- | --- | --- |
| Bubble temperature | 0.0982 K | 0.3166 K | x_ethanol = 0.015, predicted 368.9534 K, measured 369.2700 K |
| Ethanol vapor mole fraction | 0.00430 | 0.01295 | x_ethanol = 0.015, predicted 0.15195, measured 0.13900 |

The series passes both screening limits. XML SHA-256
`9ade859fa53eb443bb45a75f14e742056f0fb137bbfdf5fdc93dfa05d614a11a`.
Reproduce with `python scripts/validate_ethanol_water_reference.py PATH_TO_XML`.
A separate pressure-only holdout comparison for those fitted parameters follows.


### Pressure-only holdout for the fitted NRTL and UNIQUAC models

Cristino et al., DOI
[10.1016/j.fluid.2012.12.014](https://doi.org/10.1016/j.fluid.2012.12.014),
report a 12-point isothermal water/ethanol series at 363.3 K, within the
307.44-369.27 K fit range. The committed NRTL and UNIQUAC metadata explicitly
record this series as excluded from the fit because the vapor compositions are
not available. The metadata therefore support its exclusion from the project's
pair fit. Possible reuse of the underlying measurements by other cited articles,
and the provenance of the pure-component vapor-pressure correlations, have not
been audited.

Before calculation, mean and maximum absolute relative pressure error limits
of 5% and 10% were selected. At 95% confidence the XML reports expanded
uncertainties of 0.0027-0.082 in measured water mole fraction for these points.
The XML does not provide a corresponding point-specific pressure uncertainty.
The pressure screening limits are model-assessment thresholds, not propagated
experimental uncertainties. The measured composition is used as the input to
the model's bubble-pressure calculation.

| Model | Mean relative pressure error | Maximum | Maximum-error point |
| --- | --- | --- | --- |
| NRTL | 0.8205% | 1.6998% | x_water = 0.5777, predicted 146.8544 kPa, measured 144.4 kPa |
| UNIQUAC | 0.8268% | 1.7588% | x_water = 0.5777, predicted 146.9397 kPa, measured 144.4 kPa |

Both pass the stated screening thresholds. The measurements do not include vapor
compositions, so this comparison does not validate those predictions.
XML SHA-256
`cb6925095fa52756bbacd13b7a0152e0e7e34df00fcb83a83be46258d555459c`.
Run `python scripts/validate_ethanol_water_holdout.py PATH_TO_XML --method nrtl`
or pass `--method uniquac`.

## Water/1-butanol liquid-liquid equilibrium

Pirahmadi et al., DOI
[10.1016/j.fluid.2010.09.013](https://doi.org/10.1016/j.fluid.2010.09.013),
report two zero-salt tie lines at 101 kPa and 298.15 and 308.15 K within a
ternary salt study. The selected points have reported water and butanol mass
fractions summing to one in both phases at the published precision. Their
reported 95% expanded uncertainties in butanol mass fraction are 0.0011-0.0012
for the aqueous phase and 0.0020-0.0021 for the organic phase.

Before comparison, absolute butanol mass-fraction error limits were set at 0.05
for the mean of the four phase measurements and 0.10 for the maximum. These
screening limits are broader than measurement uncertainty. UNIFAC Dortmund gives:

| Temperature | Phase | Measured | Predicted | Absolute error |
| --- | --- | --- | --- | --- |
| 298.15 K | aqueous | 0.0703 | 0.110105 | 0.039805 |
| 298.15 K | organic | 0.7932 | 0.881785 | 0.088585 |
| 308.15 K | aqueous | 0.0603 | 0.107526 | 0.047226 |
| 308.15 K | organic | 0.7848 | 0.877616 | 0.092816 |

The mean error is 0.067108 and the maximum is 0.092816 at 308.15 K in the
organic phase. **The UNIFAC mean-error screening criterion fails.**

An independent solution of the two liquid-phase activity-equality equations
using the same Dortmund UNIFAC group model agrees with Fluxsheet's flash
compositions within 1e-6 mole fraction at both temperatures, for three overall
compositions and for reversed component order. The test is in
`tests/test_thermo_models.py`. This verifies consistency with that particular
implementation of the activity model, and supports the conclusion that the
reported mismatch is a limitation of this model/parameter set in this test.
It is not an independent validation of the group parameters themselves.

With the same limits, NRTL gives mean 0.047580 and maximum 0.061067 at 298.15 K
in the organic phase. UNIQUAC gives mean 0.032945 and maximum 0.039465 at
298.15 K in the aqueous phase. Both pass these broad screening criteria, while
their maximum deviations exceed the reported measurement uncertainty.
All methods were evaluated at overall 1-butanol mole fraction 0.2, within the
experimentally reported phase split. These two tie lines do not establish
accuracy elsewhere. Fluxsheet reports ChemSep as the parameter source for
both NRTL and UNIQUAC for water/1-butanol. Their original fit datasets, and the
provenance of the UNIFAC group parameters and pure-component properties, have not
been audited. Independence from this measurement is therefore unconfirmed.

XML SHA-256
`52c1b2b4d1175ff67bdf3041b8b89639a0d1d893956a956355c2e42aad35eb13`.
Reproduce with `python scripts/validate_butanol_water_reference.py PATH_TO_XML`.
Use `--method nrtl` or `--method uniquac` for the other models. The script exits
nonzero when the fixed screening limits fail.

### Additional binary liquid-liquid reference at 288.15 and 318.15 K

Chen and Lin, DOI [10.1016/j.fluid.2007.07.011](https://doi.org/10.1016/j.fluid.2007.07.011),
report two binary water/1-butanol phase compositions at each temperature. Before
the calculations the source hash, table selectors and the same 0.05 mean and
0.10 maximum absolute mass-fraction screening limits were fixed in
`docs/chen2007-lle-manifest.json`. Each phase has two points. ThermoML assigns
101 kPa to the aqueous-phase table and 101.325 kPa to the organic-phase table;
each prediction below uses its table's stated pressure, rather than treating the
different pressures as a single exactly matched tie line.

| Method | Mean absolute error | Maximum error | Maximum-error point | Screening |
| --- | --- | --- | --- | --- |
| UNIFAC Dortmund | 0.059933 | 0.082525 | 318.15 K, organic phase | fail |
| NRTL, ChemSep | 0.054925 | 0.079318 | 288.15 K, organic phase | fail |
| UNIQUAC, ChemSep | 0.031141 | 0.049986 | 288.15 K, aqueous phase | pass |

The source reports point-specific expanded uncertainties in mole fraction at
95% confidence. The script converts them to the mass-fraction basis by linear
uncertainty propagation; these experimental uncertainties and the broader
screening limits are different quantities. As for the earlier LLE comparison,
the underlying ChemSep and UNIFAC parameter sources have not been traced to
their original fit data, so independence of this source is not established.
Four phase observations do not validate a continuous operating range.

Reproduce each method from the original XML with
`python scripts/validate_chen_lle.py PATH_TO_XML --method unifac` (or `nrtl`
or `uniquac`). The source hash is checked before parsing. The script returns
1 for an exceeded screening limit and 2 for invalid or unavailable input.

## Reactive columns and non-ideal Gibbs reactors

`docs/validation-campaign-features.json` registers two literature comparisons for
the equilibrium reactive stages and the Gibbs reactor on fugacities. The reference
values are small published tables kept in the manifest, so nothing is downloaded.
Run `python scripts/validate_features.py --output-dir PATH`. It writes
`features.json` and `features.md` and returns 0 when every limit holds, 1 for a
screening failure and 2 for an invalid manifest or a case that cannot be solved.
`tests/test_validation_campaign.py` runs the campaign and checks the numbers on
this page.

The limits below were fixed in the manifest before the script was first run. For
ammonia this is not a blind comparison: the reactor tests had already compared five
of the points with looser tolerances while the feature was being built. No model
parameter was fitted to either source.

### Ammonia synthesis equilibrium

Larson and Dodge, J. Am. Chem. Soc. 45 (1923) 2918, measured equilibrium
ammonia in a 3:1 hydrogen-nitrogen mixture at 10–100 atm; Larson, J. Am. Chem.
Soc. 46 (1924) 367, extended the series to 300–1000 atm. The manifest carries the
standard tabulation of both at 400, 450 and 500 °C. **These values have not yet
been checked against a copy of the original articles.** The Gibbs reactor solves
each point isothermally with `fugacities: property_method` and Peng-Robinson
without fitted binary parameters. The ideal-gas result is shown for contrast and
is not screened.

| Series | Points | Limits mean/max | Mean error | Maximum error | Maximum-error point | Screening |
| --- | --- | --- | --- | --- | --- | --- |
| 10–300 atm, Peng-Robinson | 15 | 1/2 mol-% | 0.26 mol-% | 1.47 mol-% | 450 °C, 300 atm: 34.35 predicted, 35.82 measured | pass |
| 10–300 atm, ideal gas | 15 | — | 1.14 mol-% | 5.10 mol-% | 400 °C, 300 atm: 41.90 predicted, 47.00 measured | — |
| 600–1000 atm, Peng-Robinson | 6 | 3/6 mol-% | 4.40 mol-% | 6.54 mol-% | 450 °C, 1000 atm: 63.15 predicted, 69.69 measured | **fail** |
| 600–1000 atm, ideal gas | 6 | — | 14.0 mol-% | 18.4 mol-% | 400 °C, 1000 atm: 61.46 predicted, 79.82 measured | — |

From 100 atm on, the fugacities are closer to the measurement than the ideal gas
at every point. Above 300 atm Peng-Robinson without binary parameters still
predicts too little ammonia, by up to 6.5 points at 1000 atm. The reactor is
therefore screened for ordinary synthesis pressures, not for the 600–1000 atm
range. Pure-component formation data and the equation of state are both part of
this comparison; their contributions have not been separated.

### Methyl acetate reactive distillation

Mekala, Data in Brief 27 (2019) 104262,
DOI [10.1016/j.dib.2019.104262](https://doi.org/10.1016/j.dib.2019.104262),
Tables 3 and 4, report measured temperatures on eight stages and liquid
compositions on five stages of a 30-stage column with ten catalytic stages. The
article cites the measurements from earlier work by the same group. The values
were checked against the PMC full text. The registered limits are the
preregistered column limits of the equipment campaign: 2/5 K for temperature and
0.02/0.05 mole fraction for composition.

The Fluxsheet case is `examples/methyl_acetate_equilibrium_column.json`: NRTL,
chemical equilibrium on stages 11–20 with the activity-based K(T) of Huss et al.,
reflux ratio 1.9 and equimolar feeds of 10.2 mol/h. The article does not give the
distillate rate or the feed temperatures; the example assumes half the feed and
320 K.

| Quantity | Points | Limits mean/max | Mean error | Maximum error | Maximum-error point | Screening |
| --- | --- | --- | --- | --- | --- | --- |
| Stage temperature | 8 | 2/5 K | 4.51 K | 13.09 K | stage 30: 371.09 K predicted, 358 K measured | **fail** |
| Liquid mole fraction | 20 | 0.02/0.05 | 0.087 | 0.28 | stage 30, water: 0.951 predicted, 0.669 measured | **fail** |

**The comparison fails both screening limits.** The measured column did not reach
chemical equilibrium: its bottoms still hold 0.10 acetic acid and 0.21 methanol,
and stages 13–25 sit within 1 K of 339 K. An equilibrium model cannot represent a
kinetically limited packing, and the unreported distillate rate and heat losses
remain unknown. The comparison therefore establishes that equilibrium stages are
not a substitute for the kinetic model on this column. It does not test whether
the stage equations themselves are correct; that is verified separately in
`tests/test_column_side_streams.py`, where the activities on every reactive stage
are recomputed from the NRTL model and meet K(T) to a relative 1e-8.

## Dynamic column and batch reactor costing

No experimental reference has been registered for these two features. The checks
below are verification only.

- **Dynamic column.** In `tests/test_dynamic_column.py` a 24-stage
  benzene/toluene column starts at its steady MESH solution and stays there
  undisturbed. After a feed-composition step under level control it comes to rest
  again with nothing vented or accumulated (distillate plus bottoms equal to the
  feed to a relative 1e-5). It then matches the steady column with the same reflux,
  distillate and feed within 2e-5 mole fraction and 0.1 % in both duties. The
  candidate column sources in the equipment qualification table below apply here
  as well; experimental transients remain open.
- **Batch reactor costing.** The vessel uses Turton's jacketed agitated reactor,
  Table A.1, the same correlation as the CSTR. `tests/test_validation_campaign.py`
  evaluates that correlation by hand from the tabulated coefficients. The fill
  fraction of 0.8, the default turnaround of 1 h and the parallel-vessel rule are
  engineering assumptions. They have not been compared with vendor quotes or
  plant costs, and Turton's correlations themselves are study-grade (about ±30 %).

## Equipment-data qualification

The [version-2 equipment campaign](equipment-campaign.md) now generates offline
qualification reports from `docs/equipment-campaign.json`. It returns code 2 for
the currently incomplete sources and preserves completed comparisons in a mixed
campaign. Qualified single-liquid water vessel cases can now run the transient
simulator and compare temperature and level against registered measurements.

The expanded [source qualification and data request sheets](equipment-reference-qualification.md)
record additional candidates and the inputs that remain unverified. The numerical
comparison helpers and stirred-vessel adapter are tested with synthetic data.
The vessel fixtures compare the simulator against independently integrated
enthalpy balances for heat input and feed temperature steps.

The following source assessment is a qualification step, not an experimental
validation of these equipment models. No parameter was fitted to a test output.

| Family | Candidate examined | Current disposition |
| --- | --- | --- |
| Dynamic vessel | [Thornhill et al. CSTH](https://zenodo.org/records/10093059) | The published archive contains Simulink models and measured disturbance/noise sequences. It does not identify the outputs as measured step-response records. Simulation outputs driven by measured noise are not independent experimental step responses; the required two measured input/output steps remain open. |
| Reaction kinetics | [Zeki et al.](https://iasj.rdd.edu.iq/journals/uploads/2024/12/16/b18ca2d92132373e3a8bcef9c13164f8.pdf) | The paper reports 50–60 °C batch experiments, 95% ethanol and one drop of sulfuric acid; initial pH differs between feed ratios. Conversion trajectories appear as plots and the tabulated rate constants were fitted to the same experiments. Quantitative catalyst loading, measured trajectory tables and independently held-out runs have not been established, so the paper is insufficient for the specified independent comparison. |
| Stationary column | [Werle et al.](https://doi.org/10.1016/j.cep.2009.10.010) and [continuous distillation dataset](https://zenodo.org/records/17628963) | The Werle abstract reports a 13-tray ethanol/water pilot column; complete stationary boundary conditions and measured efficiencies have not been obtained. The open dataset contains two-column butanol/water operation, but two of the three normal-operation metadata files do not supply numerical feed compositions. A compatible single-column model with independently specified efficiency and heat losses has not been established. |

Fixed-pressure vessel dynamics do not represent column pressure transients. A
study requiring liquid-liquid column behavior, activity-based kinetics, or
unknown apparatus parameters must not be silently treated as a compatible case.
Each family remains experimentally open until its prespecified source and model
qualification criteria are met.

The preregistered engineering screening limits remain: for a column, mean/max
composition error of 0.02/0.05 mole fraction, temperature error of 2/5 K, and
relative duty error of 10/20%; for kinetics, mean/max absolute conversion error
of 0.03/0.08; and for dynamics, trajectory RMSE and final-value error each no
more than 10% of the measured step amplitude. A dynamic step must exceed five
times the baseline noise. These are model-screening limits, not measurement
uncertainties, and cannot be evaluated until suitable measurements are available.
