# Data sources

A simulator is only as good as its property data. This page records where Fluxsheet's data
comes from today and how we intend to grow it, using **openly licensed sources only**.

## Today

All pure-component and mixture data is provided by the
[`chemicals`](https://github.com/CalebBell/chemicals) and
[`thermo`](https://github.com/CalebBell/thermo) libraries (MIT). They compile published,
openly available datasets, including:

| Data | Typical sources inside `chemicals`/`thermo` |
|---|---|
| Identifiers (~76,000 compounds) | PubChem-derived metadata |
| Critical constants, acentric factor | IUPAC, Yaws-type public compilations, CRC-style tables, estimation methods |
| Vapor pressure, heat capacity, liquid density | Published correlations (DIPPR-form equations with public coefficients, Perry's, VDI-style data where licensed) |
| Peng-Robinson binary interaction parameters | ChemSep database |
| NRTL and UNIQUAC binary parameters | Fluxsheet's own fits to NIST ThermoML VLE data (below) for the pairs in `src/fluxsheet/data/binary_parameters/`, otherwise the ChemSep database; pairs in neither are fitted to UNIFAC (Dortmund) at run time and flagged as estimates in the solve warnings |
| Wilson binary parameters | ChemSep database, UNIFAC estimates for gaps (as above) |
| UNIFAC groups and parameters | Published Dortmund UNIFAC tables (public parameter set) |

A compound is marked `simulation_ready` in the editor when at least critical temperature and
pressure are available. Each library documents its sources in detail; consult them before
relying on a specific value.

## Henry constants (dissolved gases)

Activity-coefficient methods treat light gases as Henry components (see the README). Per gas
and solvent the first available source is used:

| Source | Coverage | Notes |
|---|---|---|
| IAPWS G7-04 (2004), Fernández-Prini, Alvarez and Harvey, J. Phys. Chem. Ref. Data 32 (2003) 903 | He, Ne, Ar, Kr, Xe, H2, N2, O2, CO, CO2, H2S, CH4, C2H6, SF6 in water, 0 to about 350 °C | The guideline's coefficients are public; Fluxsheet refits its equation to thermo's `ln H = A + B/T + C ln T + D T` form between 0 and 250 °C (within 1 %) |
| Sander (2015), Atmos. Chem. Phys. 15, 4399 (CC BY 3.0), as shipped with `thermo` | About 950 solutes in water, temperature-dependent | Single van't Hoff term; poor above about 60 °C |
| ChemSep Henry table (Artistic License 2.0, via `thermo`) | 41 gases in water | Taken over Sander only where both agree within 10 % at 25 °C (ChemSep's N2 and H2 values are off by a factor of 5 to 10) |
| Sander (2015) constants | About 4400 solutes in water, 25 °C only | No heat of solution |
| `thermo`'s Peng-Robinson predictions ("PRTranslated Henry Best") | Gases in about 160 organic solvents | Predicted, not measured; never used for water |
| None: the solvent is left out of the mixing rule | Any solvent without a constant | The other solvents' constants apply in liquids rich in it; always warned (an extrapolated vapor pressure would be orders of magnitude off) |

Literature checks in `tests/test_henry_custom.py`: solubilities of O2, N2, CO2, H2, CH4 and Ar
in water at 25 °C within 5 % of Wilhelm, Battino and Wilcock (1977) and Carroll et al. (1991).

## Estimated properties of user-defined compounds

| Property | Method |
|---|---|
| Tb, Tc, Pc, Vc, Tm, ideal-gas Hf and Gf, ideal-gas Cp(T) | Joback and Reid (1987) group contributions (`thermo.group_contribution.joback`); Cp polynomial used from 280 to 1000 K, extended linearly |
| Joback groups from SMILES | RDKit fragmentation in `thermo` (optional dependency) |
| Acentric factor | Lee-Kesler vapor-pressure equation through Tb, Tc and Pc |
| Vapor pressure | Lee-Kesler or Ambrose-Walton from Tc, Pc and omega, unless Antoine constants are given |
| Liquid volume | COSTALD (Hankinson-Thomson) |
| Cp without groups | Lastovka-Shaw from the formula |
| UNIQUAC r, q | Original UNIFAC group volumes and areas |
| Viscosity, conductivity, surface tension | `thermo`'s corresponding-states methods from the constants |

Joback's typical errors are a few percent in Tb and Tc and up to about 10 % in Pc; an
estimated Tb can put vapor pressures off by a factor of two, so results say when it was
estimated. Given constants always win over estimates.

## Binary parameters from NIST ThermoML

`src/fluxsheet/data/binary_parameters/nrtl.json` and `uniquac.json` hold NRTL and UNIQUAC
parameters that Fluxsheet regressed itself to vapor-liquid equilibrium (VLE) data from the
[NIST TRC ThermoML Archive](https://trc.nist.gov/ThermoML/). The archive holds the
experimental data of articles in J. Chem. Eng. Data, J. Chem. Thermodyn., Fluid Phase
Equilib., Thermochim. Acta and Int. J. Thermophys., extracted by NIST TRC staff in the IUPAC
ThermoML format.

### Terms of use (as read on 2026-10-29)

- The data.nist.gov record of the archive ([doi:10.18434/mds2-2422](https://doi.org/10.18434/mds2-2422))
  lists the [NIST open license](https://www.nist.gov/open/license): NIST data are generally
  not subject to copyright in the United States; derived works may be made and distributed,
  should state that the data were changed, and should acknowledge NIST as the source.
- The [archive page](https://www.nist.gov/mml/acmd/trc/thermoml/thermoml-archive) says the
  files are available "with permission of the journal publishers", and the record carries a
  disclaimer: neither the publishers nor NIST/TRC warrant the correctness of the data or
  accept liability for their use; users should check them against other sources.

Because the files reproduce data from copyrighted journal articles and the publishers'
permission is not stated to extend to redistribution, **this repository does not contain any
ThermoML file or raw data point.** It stores only derived results (fitted parameters, fit
statistics, data ranges and point counts) together with the DOI and citation of every source
article, so each value can be traced to the original measurement. Acknowledgement: *data from
the NIST Thermodynamics Research Center ThermoML Archive, NIST, U.S. Department of Commerce;
parameters derived by Fluxsheet (changed from the original data by regression).*

### How the parameters are made

`fluxsheet fit-binary A B --model nrtl|uniquac` (code in `src/fluxsheet/fitting/`):

1. **Find and download.** The archive's search API (`https://trc.nist.gov/ThermoML-API/objects`,
   the service behind its search page) is queried with both compounds' standard InChIKeys;
   every article reporting binary phase-equilibrium compositions is downloaded from
   `https://trc.nist.gov/ThermoML/<DOI>.xml`. Downloads are cached outside the repository
   (`%LOCALAPPDATA%\fluxsheet\thermoml`, `~/.cache/fluxsheet/thermoml`, or
   `FLUXSHEET_THERMOML_CACHE`); only the two official hosts are contacted.
2. **Parse safely.** Files are untrusted input: `defusedxml` (no DTDs or entities), a size
   limit, validated DOIs before they become file names, range checks on every number. Binary
   data sets with a liquid and a vapor phase and T, P, x (and y) are extracted; ThermoML often
   stores one T-P-x-y set as two tables (T-x and y-x), which are joined row by row. Compounds
   are mapped to CAS numbers by InChIKey (`fluxsheet.compounds.lookup`), checked by formula.
3. **Select data.** Complete T-P-x-y sets up to 300 kPa (the property packages assume an ideal
   vapor), at least 3 mixture points; pure-component points are left out.
4. **Consistency.** Each set must pass the *point test* (Van Ness et al. 1973, Fredenslund
   criterion): a four-term Redlich-Kister G<sup>E</sup> model fitted to T (isobaric) or P
   (isothermal) alone by Barker's method must predict the measured vapor compositions within
   a mean |Δy| of 0.01. The Redlich-Kister/Herington *area test* is computed and stored too,
   but does not decide (it is unreliable for nearly ideal mixtures).
5. **Fit.** Modified Raoult's law with thermo's vapor-pressure correlations, the same model the
   property packages use. Least squares (`scipy.optimize.least_squares`) on bubble T
   (isobaric) or bubble P (isothermal) and on y, scaled by 0.2 K, 0.5 % and 0.005. Always
   fitted: b<sub>12</sub>, b<sub>21</sub>; also a<sub>12</sub>, a<sub>21</sub> when the data
   span at least 30 K. NRTL's α is fixed at 0.3 (0.2 if UNIFAC predicts a liquid-liquid split).
   Several starting points (zero, the UNIFAC estimate, ChemSep's set); the best is kept.
6. **Outliers.** A point is dropped when one of its scaled residuals exceeds three times its
   robust spread (1.4826 × median absolute residual) and three measurement uncertainties, or
   when its liquid lies inside the model's liquid-liquid gap (the reported x is then an overall
   composition). The fit is repeated once without them; the count is stored.
7. **Report.** AAD and RMS deviations in T (K), P (%) and y on the points used, and the same
   statistics for the ChemSep set and for the UNIFAC estimate on the same points.

Each pair in the JSON files records: parameters (`tau_ij = a_ij + b_ij/T` for NRTL,
`tau_ij = exp(a_ij + b_ij/T)` for UNIQUAC, the `thermo` conventions), α, which parameters were
fitted, the T and P range, number of points used and rejected, fit statistics, the comparison,
the DOI and citation of every source article with the data sets used, the data sets that were
not used and why, and the date of the fit.

### How they are used

`PropertyPackage(..., prefer_fitted=True)` (the default) takes a pair's ThermoML fit over
ChemSep's set; with `prefer_fitted=False`, ChemSep comes first and fits only fill its gaps.
Every solve result lists the source of each pair in `parameter_sources` (`ThermoML fit` with
references and data range, `ChemSep`, `UNIFAC estimate`; units with their own method add
entries tagged with that `method`), and so does the phase-diagram view. Outside the stored
T range a fit is an extrapolation, like any regressed parameter set.

### Seed set (fitted 2026-10-29)

Water/ethanol, water/methanol, water/1-propanol, water/2-propanol, acetone/methanol,
ethanol/ethyl acetate, methanol/ethyl acetate and benzene/toluene, each for NRTL and UNIQUAC
(`scripts/fit_seed_pairs.py` refits all of them and prints the comparison with ChemSep and
the UNIFAC estimate on the same points). ChemSep has no NRTL set for water/1-propanol and
water/2-propanol, so the fits replace UNIFAC estimates there.

Pairs looked at but not stored:

- **water/1-butanol**: only one consistent VLE set (13 points). The fit describes it far
  better than ChemSep, but predicts about 40 mol% water in the butanol-rich liquid at 25 °C
  (measured: about 51 %). Partially miscible pairs need vapor-liquid and liquid-liquid data
  fitted together; `fit-binary` prints a note when UNIFAC predicts a split.
- **methanol/ethanol** (only T-P-x data, no vapor compositions), **acetone/water** (only
  Henry's-law and activity-coefficient sets), **toluene/methylcyclohexane** (T-x-y without
  pressure): no complete T-P-x-y set in the archive that the selection rules accept.

### Refitting and adding pairs

```bash
pip install -e ".[dev]"
fluxsheet fit-binary water ethanol --model nrtl            # prints data sets, tests, deviations
fluxsheet fit-binary water ethanol --model nrtl --write    # updates nrtl.json
```

`--max-files` limits how many articles are read (most binary data first; 0 = all),
`--offline` reuses the cache, `--alpha` sets NRTL's α, `--no-consistency` skips the point-test
filter. Review the diff of the JSON file like code: check the deviations against ChemSep, the
rejected data sets, and the sources.

## Plan

1. **Make provenance visible.** Show, per compound and property, which correlation/source was
   used and its valid temperature range, and warn when a flowsheet extrapolates.
2. **Binary parameters from experimental data.** *Started:* NRTL/UNIQUAC parameters regressed
   to [NIST ThermoML](https://www.nist.gov/mml/acmd/trc/thermoml) data, see above. Next: more
   pairs, liquid-liquid data, and per-pair validity ranges checked during a solve.
3. **Estimation as a fallback.** *Started:* Joback and Lee-Kesler estimates for user-defined
   compounds (above). Next: Constantinou–Gani, and UNIFAC groups from SMILES.
4. **Bring your own data.** *Started:* compounds with their own constants and Antoine
   equation in the flowsheet file. Next: binary parameters in a plain file, and plugging in
   commercial databases users hold licenses for (e.g. DIPPR, DDB) without redistributing them.

## What we do not do

We do not copy data out of commercial simulators or licensed databases. Contributions of data
must state their source and license.
