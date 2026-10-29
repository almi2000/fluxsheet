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
| UNIFAC groups and parameters | Published Dortmund UNIFAC tables (public parameter set) |

A compound is marked `simulation_ready` in the editor when at least critical temperature and
pressure are available. Each library documents its sources in detail; consult them before
relying on a specific value.

## Plan

1. **Make provenance visible.** Show, per compound and property, which correlation/source was
   used and its valid temperature range, and warn when a flowsheet extrapolates.
2. **Binary parameters from experimental data.** Use [NIST ThermoML](https://www.nist.gov/mml/acmd/trc/thermoml)
   (open experimental data from journal articles) to regress NRTL/UNIQUAC parameters, stored
   as versioned, reviewable data files in this repository together with the fit quality.
3. **Estimation as a fallback.** Group-contribution methods (Joback, Constantinou–Gani, UNIFAC)
   for compounds without measured data, clearly flagged as estimated.
4. **Bring your own data.** Let users add compounds and parameters in a plain file, and plug in
   commercial databases they hold licenses for (e.g. DIPPR, DDB) without redistributing them.

## What we do not do

We do not copy data out of commercial simulators or licensed databases. Contributions of data
must state their source and license.
