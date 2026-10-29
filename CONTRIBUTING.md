# Contributing

Thanks for your interest! OpenFlowsheet is young, so there is plenty to do: see the roadmap in
the README, or open an issue to discuss an idea first.

## Development setup

```bash
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[server,dev]"
pytest
ruff check src tests

cd web && npm install && npm run dev   # with `openflowsheet serve` running
```

## Guidelines

- **Every physical result needs a test**, ideally against a hand calculation or a published
  value, plus mass/energy balance checks.
- **Errors must explain themselves.** Raise `SimulationError(message, where=..., hint=...)`
  with a concrete suggestion.
- **SI inside, units at the edges.** Convert with `quantities.to_si`; never store mixed units.
- **Data needs a source.** Property data contributions must state origin and license; we only
  accept openly licensed data (see `docs/data-sources.md`).
- Keep pull requests focused; describe what changed and how you verified it.
