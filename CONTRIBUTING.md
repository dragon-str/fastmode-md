# Contributing to fastmode-md

Thank you for your interest. This project welcomes bug reports, feature
requests, and pull requests.

## Reporting a problem

Open an issue at
<https://github.com/dragon-str/fastmode-md/issues>. Include:

- the `fastmode-md` version (`pip show fastmode-md` or the Git commit),
- the GROMACS version,
- the command you ran and the system size,
- the complete error output, if the tool stopped.

## Development setup

The unit tests do not need GROMACS. Only NumPy and pytest are required.

```bash
python -m venv .venv
. .venv/bin/activate
pip install numpy pytest
python -m pytest -q
```

The integration workflow (`fastmode audit`, `fastmode selfcheck`) does need a
working GROMACS install, because it calls `grompp`, `mdrun`, `energy` and
`rdf`. Point the tool at the engine with `--gmx` when it is not on `PATH`.

## What to change

- Keep the style of `fastmode.py`: plain Python 3, standard library plus NumPy.
- Add a unit test in `tests/test_fastmode.py` for every new pure function.
- Do not add a new runtime dependency without a clear reason.
- Do not lower a validation limit to make a candidate pass. A failed check is
  a result, not a problem to hide.
- Update `CHANGELOG.md` under an `Unreleased` heading.

## Pull requests

1. Fork the repository and branch from `master`.
2. Run `python -m pytest -q` and make sure it passes.
3. Describe the change and the evidence for it in the pull request.
4. Link the issue the pull request closes, if one exists.

## Governance and support

The project is maintained by Daniel Reda. There is no paid support. Issues and
pull requests are answered on a best-effort basis.

## License

By contributing, you agree that your contributions are licensed under the MIT
License.
