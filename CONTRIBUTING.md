# Contributing

## Setup

```sh
pyenv install --skip-existing "$(cat .python-version)"
python -m venv venv && source venv/bin/activate
python -m pip install -r requirements.txt
pre-commit install
```

## Checks

Code standards are enforced by pre-commit (black, ruff, mypy, whitespace and
file hooks) on every commit and in CI.

```sh
pre-commit run --all-files
pytest
mkdocs build --strict
```

GitHub Actions workflows run locally with [act](https://github.com/nektos/act)
and Docker. `.actrc` selects the runner image, amd64 architecture, a fixed
package version, and a PATH that keeps node available after setup-python.
Artifact upload steps fail under act because its artifact server does not
support the current upload protocol; every other step is verified.

```sh
act -j lint
act -j test --matrix os:ubuntu-latest --matrix python-version:3.14
act -j extras
act -W .github/workflows/docs.yml -j build
```

## Rules

- Tests first. Each test checks a closed-form answer, a synthetic layout with
  a known structure, or a documented contract. Tests live in `tests/`, not in
  `src/`.
- Every numeric threshold or model switch is a field on `ExtractionParams`,
  `Stack`/`Material`, or `SolveOptions`, with a docstring and a row on the
  [parameters page](https://cao-rqm-lab.github.io/connflex-core/parameters/).
- Core code imports only shapely, numpy, scipy, and networkx. Readers live in
  `adapters/` behind optional extras.
- The package ships `py.typed` and mypy runs in strict mode with type stubs
  for every dependency. A new dependency brings its stubs into the `dev`
  extra and the mypy hook in `.pre-commit-config.yaml`. Only gdsfactory and
  klayout readers are exempt, through mypy overrides in `pyproject.toml`,
  because the hook does not install those extras.
- SI units inside the library. Readers convert on the way in.
- Every pipeline stage runs inside `stage()` from `connflex_core.progress`,
  which logs start and finish through the `connflex_core` logger and reports
  to the optional `progress` callback. Loops over polygons, bodies, or
  candidates advance once per item. No `print`, verbose flags, or logging
  handlers in the library.
- No GUI, plotting, or notebook code.

## Pull requests

Branch from `main`. A pull request needs green CI, one review, and a line
under **Unreleased** in `CHANGELOG.md`.

## Releases

1. Open a pull request that moves **Unreleased** in `CHANGELOG.md` under the
   new version.
2. After merging, create an annotated tag `vX.Y.Z` and a GitHub release from
   it.
3. The release workflow builds the distribution, checks that the tag matches
   the built version, smoke-installs the wheel, and publishes to PyPI after a
   maintainer approves the `pypi` environment. Releases marked as pre-release
   go to TestPyPI.

hatch-vcs derives the version from the tag.
