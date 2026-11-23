# Notes for coding agents

Fluxsheet is a chemical process simulator: a Python core (`src/fluxsheet`) and a React
editor (`web/`). Read these first:

- [docs/architecture.md](docs/architecture.md): how the code is organized, how to add a unit
  operation, plug-ins.
- [CONTRIBUTING.md](CONTRIBUTING.md): the engineering rules (tests against independent
  values, explanatory errors, SI inside).
- [docs/ROADMAP.md](docs/ROADMAP.md) (German): open work packages, who owns which files,
  and the rules for two agents working in parallel.
- [docs/model-validation.md](docs/model-validation.md): what is verified and what is
  validated against measurements. Do not describe a model as validated unless it is listed
  there.

## Commands

Python 3.14. The checks CI runs (`.github/workflows/ci.yml`):

```
pip install -e ".[server,dev]"
ruff check src tests          # line length 100
pytest -q                     # the full suite takes about 8 minutes

cd web
npm ci
npm run build                 # tsc --noEmit, then vite build into web/dist
npm test                      # vitest
```

- Run single test files while you work and the full suite once at the end, in the
  background.
- **Build the editor before running pytest, never at the same time.** `tests/test_server.py`
  serves `web/dist` (it skips if there is no build), and a build replacing it mid-run makes
  that test fail.
- `ruff check` also on `scripts/` and `examples/plugin/` when you change them; CI does not.

## Working in a git worktree

Agents work on their own branch, often in a worktree under `.claude/worktrees/`.

- The editable install points at the main checkout, so from a worktree run Python with the
  worktree's sources first: `PYTHONPATH=src python -m pytest ...` (PowerShell:
  `$env:PYTHONPATH = "src"`). Check with `python -c "import fluxsheet; print(fluxsheet.__file__)"`.
- A worktree has no `web/node_modules`. Link the main checkout's instead of installing:
  `cmd /c mklink /J web\node_modules <main checkout>\web\node_modules`. Before deleting the
  worktree, remove the link with `cmd /c rmdir web\node_modules` (this removes only the
  link). Then check that the main checkout's `web/node_modules` is still complete.

## Parallel work and publishing

- Never commit, stash, discard or overwrite changes you did not make, in any checkout. The
  git stash list is shared by all worktrees.
- Merges into `main` and pushes only after the maintainer approves them.
- **Do not publish anything without being asked:** no `v*` tag (it starts
  `.github/workflows/release.yml`, which creates a GitHub release and may upload to PyPI),
  no PyPI upload, and do not start the online demo.
- Shared files (`server.py`, `cli.py`, `units/__init__.py`, `web/src/types.ts`,
  `web/src/api.ts`, `README.md`) take additive changes only, so parallel branches merge
  cleanly.
- `CHANGELOG.md`: one entry under `[Unreleased]` when a package is finished.

## Conventions the tests enforce

- **Balance diagnostics.** `tests/conftest.py` records every successful solve. A balance
  `warning` or `not_checkable` fails the test unless the test declares it in the
  `expected_balance_exceptions` fixture, with the exact key
  `(scope, quantity, component, status)` and a reason. A declared exception that does not
  occur also fails.
- **Independent checks.** A new physical result is tested against a hand calculation, a
  published value or an independent solution of the same equations, not only against its
  own output.
- **Units register themselves.** Every `UnitOperation` subclass with a `type` in `units/` is
  picked up; the editor reads its schema from `/api/unit-types`.

## Pitfalls

- Do not name a script or scratch file after a standard-library module (`bisect.py`,
  `random.py`, ...) in the repository root or `tests/`; it shadows the real module for the
  code under test.
- Files use LF line endings (`.gitattributes`). On Windows, `Path.write_text` writes CRLF
  unless you pass `newline="\n"`.
- Tests import the validation scripts as `scripts.*` (`pythonpath = ["."]` in
  `pyproject.toml`).
- Keep scratch files out of the repository; use a temporary directory.
