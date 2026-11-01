# Contributing

Thanks for your interest! Fluxsheet is young, so there is plenty to do: see the roadmap in
the README, or open an issue to discuss an idea first.

## Development setup

```bash
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[server,dev]"
pytest
ruff check src tests

cd web && npm install && npm run dev   # with `fluxsheet serve` running
```

## Guidelines

- **Every physical result needs a test**, ideally against a hand calculation or a published
  value, plus mass/energy balance checks.
- **Test the regimes, not just the example.** Where physics gives the answer for a whole range
  (phase from the critical point and vapor pressure, a bubble point below its dew point), check
  it over a grid of conditions and compounds, as `tests/test_phase_physics.py` does: light gases
  at high pressure, near-critical and dense states are where models break. Where a unit solves
  equations that can be integrated another way, compare with that independent solution, as
  `tests/test_validation.py` does for a kinetic reactor at high pressure. Kinetics and other
  data in tests must be openly published or made up for the test, and say which.
- **Don't trust a library result unchecked.** A converged answer from `thermo` can still be the
  wrong one (a dew point returned as a bubble point). Verify the property that defines a result
  and raise a `SimulationError` when it fails, rather than pass on a plausible wrong number. Set
  library options that change results explicitly, with the reason, instead of relying on
  defaults.
- **Errors must explain themselves.** Raise `SimulationError(message, where=..., hint=...)`
  with a concrete suggestion.
- **SI inside, units at the edges.** Convert with `quantities.to_si`; never store mixed units.
- **Data needs a source.** Property data contributions must state origin and license; we only
  accept openly licensed data (see `docs/data-sources.md`).
- Keep pull requests focused; describe what changed and how you verified it.
