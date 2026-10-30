"""Regression of binary activity-coefficient parameters to measured phase-equilibrium data.

The pipeline reads vapor-liquid equilibrium data from the NIST ThermoML Archive
(:mod:`.thermoml` parses, :mod:`.archive` downloads), checks their thermodynamic consistency
(:mod:`.consistency`), fits NRTL or UNIQUAC parameters by least squares (:mod:`.regression`)
and stores them with their provenance in versioned JSON files (:mod:`.store`) that
:mod:`fluxsheet.thermo_package` reads. ``fluxsheet fit-binary`` runs it from the command line.
"""
