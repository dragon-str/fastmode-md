# AI usage disclosure

This statement records the use of generative AI in the fastmode-md software and
in its paper. It follows the disclosure requirements of the Journal of Open
Source Software.

## Tools and versions

- Anthropic Claude, used through the `opencode` command-line tool. The model was
  `deepseek-v4p1-flash` (Fireworks) for part of the work and a Claude model for
  the rest. The exact model version is recorded in the session transcript.

## Nature and scope of the assistance

AI assistance was used for the following work:

- **Code:** generation and refactoring of utility functions, and generation of
  the unit tests in `tests/`.
- **Documentation:** drafting and copy-editing of `README.md`,
  `CONTRIBUTING.md`, `CHANGELOG.md`, and the issue templates.
- **Paper text:** drafting and copy-editing of this paper and of the methods
  note in `paper/`.
- **Packaging:** preparation of `pyproject.toml`, `CITATION.cff`, the release
  archive, and the OSF deposit.

AI was **not** used for conversational interaction with JOSS editors or
reviewers.

## Core design decisions

The scientific problem, the four validation checks, the limits that define a
pass, the noise-floor method, and the choice of systems were specified by the
human author. The human author made every decision about what counts as a
correct simulation.

## Confirmation of review

The human author reviewed, edited, and validated all AI-assisted output. The
human author is fully responsible for the accuracy, originality, licensing, and
compliance of the software and the paper.
