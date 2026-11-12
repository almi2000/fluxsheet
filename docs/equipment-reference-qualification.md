# Equipment reference qualification and data requests

This is a source audit, not an experimental validation. No model parameter was
fitted. A blank field means not established from the inspected material, not that
the information necessarily does not exist. No raw measurements are distributed.

## Dynamic vessel

* [Thornhill CSTH archive](https://zenodo.org/records/10093059/files/CSTH_ZenodoPage_v1.pdf):
  Its introduction identifies a first-principles simulation, measured instrument
  nonlinearities, and experimental noise/disturbance sequences. The inspected
  archive did not identify measured input/output step-response records. Simulated
  responses are not measured validation trajectories. Missing: measured process
  inputs and outputs for two steps, their uncertainty, and evidence of wall heat
  capacity and heat losses. Run identifiers: not established. The inspected
  `CSTHSimulation.zip` SHA-256 is
  `4133f128fbfa2db3776317c572c02f6a9e61cb9e53b242945dfc722560bc2a36`.
* [Hermawan 2011](https://www.davidpublisher.com/Public/uploads/Contribute/55fa0fe6b0f06.pdf),
  sections 2–3, table 1, figures 2–3: a 10 L stirred apparatus has a reported
  5 L liquid inventory, two inlet streams and measured temperature responses.
  Table 1 gives inlet flows 14.870 and 22.634 mL/s, inlet temperatures 26 °C,
  outlet temperature 31 °C and electrical heat input 420 W. Figure 2 reports a
  heat step from 420 to 500 W; figure 3 shows an inlet-temperature disturbance.
  Graphical traces include fitted/model curves; their pointwise measured values,
  baseline noise and sensor uncertainty are not established. Wall heat capacity,
  independently measured heat losses, pressure, and measured inputs through both
  transients are not established. Moreover, the tabulated flows and 5 K temperature
  difference imply roughly 785 W sensible heat for water, substantially above the
  stated 420 W; the discrepancy requires clarification before using the table as
  an energy-balanced boundary condition. No digitization or prediction was attempted.
  Source hash: not established.

**Status: open.** Required data sheet: [ ] run identifiers and full timestamps;
[ ] two independently measured input/output steps; [ ] baseline noise and sensor
uncertainties; [ ] all flows and inlet temperatures versus time; [ ] pressure,
initial liquid and vapor inventories and geometry; [ ] measured heat supplied;
[ ] wall capacity and heat losses or independent evidence they are negligible.

## Kinetics

* [Zeki et al.](https://iasj.rdd.edu.iq/journals/uploads/2024/12/16/b18ca2d92132373e3a8bcef9c13164f8.pdf):
  the inspected publication reports a drop of catalyst, differing initial pH and
  plotted conversion trajectories; its tabulated kinetic constants were fitted to
  the same experiments. Quantitative catalyst loading, numerical trajectories,
  uncertainty and independently held-out runs were not established. Source hash:
  not established. Digitization cannot supply catalyst loading or independence.
* [Beula and Sai](https://doi.org/10.1080/00194506.2014.975761): publisher page
  returned HTTP 403 in this audit. Accessible bibliographic abstracts describe an
  isothermal experiment and fitted kinetics, but do not establish the complete
  numerical trajectories, full initial state or independent parameter provenance.
  Version/hash, qualifying runs and uncertainty remain unverified. The existence
  of a fitted kinetic model alone does not establish a holdout prediction.
* [Misikov et al.](https://rjsvd.com/0044-4537/article/view/693639): the
  available abstract describes measurements at 323.15 K. As presented it does
  not supply the required two temperatures; complete input conditions, numerical
  trajectories and independent kinetic-parameter provenance remain unverified.
  Source hash: not established. It is not admitted as a qualifying campaign.

**Status: open.** Required data sheet: [ ] at least three run identifiers and
concentration-versus-time data at two or more temperatures; [ ] initial amounts
including water and catalyst loading; [ ] liquid volume and temperature versus
time; [ ] sampling and evaporation corrections; [ ] measurement uncertainties;
[ ] forward and reverse laws, units and independently established parameters;
[ ] explicit provenance showing which whole runs were used for fitting.

## Stationary column

* [Werle et al.](https://doi.org/10.1016/j.cep.2009.10.010): the abstract identifies
  a 13-tray ethanol/water apparatus. Complete stationary boundary conditions,
  stage efficiencies, heat losses, uncertainty and three suitable measured
  operating points have not been established from the accessible material.
  Full-text/hash: not established.
* [Continuous distillation dataset](https://zenodo.org/records/17628963): the
  repository describes water, butanol/water and reactive scenarios, and separate
  concentration measurements for butanol/water. The inspected normal-operation
  metadata for runs 006 and 007 did not provide numerical feed compositions;
  independent stage efficiency, heat loss and a compatible single-column boundary
  have not been established. The related [preprint](https://doi.org/10.31224/5631)
  was not accessible through the retrieval tool. `ContinuousDistillationData.zip`
  SHA-256 is `8076e5cec0c1137ab0902c305a2603a4c4ec672308db563ff9749a9b9a7a34ec`.

**Status: open.** Required data sheet: [ ] three steady operating-point identifiers
and stable windows; [ ] all feed flow rates, compositions, temperatures and pressures;
[ ] product and reflux flow rates; [ ] physical stage count, feed and measurement
stage mapping, pressure profile; [ ] independently measured or specified stage
efficiencies and heat losses; [ ] product compositions, stage temperatures and
available duties with uncertainties; [ ] documented parameter provenance.

The existing numerical equation and boundary tests remain separate from this
experimental status. The screening criteria in `model-validation.md` are fixed
before any eventual comparison and are not measurement uncertainties.
