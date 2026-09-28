# AI usage disclosure

This statement records how generative AI was used in fastmode-md and in its
papers.

## Tools

- Claude models (Anthropic), through Claude Code and the `opencode` command-line
  tool.
- DeepSeek V4.1 Flash (`deepseek-v4p1-flash`), served by Fireworks AI, through
  `opencode`.

## What AI did

- **Code.** AI coding agents wrote most of the code in this repository,
  including `fastmode.py`, the search harness, the validation and experiment
  scripts, the machine-learning probes and the unit tests. They worked from a
  written specification of the task, the four checks and their limits.
- **Experiments.** AI agents ran the simulation sweeps and wrote the result
  reports under `results/` and `ml/`.
- **Text.** AI drafted `README.md`, `CONCLUSIONS.md`, the methods note, the JOSS
  paper and their revisions.
- **Review and correction.** In September 2026 an AI reviewer audited version
  0.1.0 and found the errors listed in `CHANGELOG.md`. An AI agent checked each
  finding against the committed data, ran the re-validation and timing
  experiments, and wrote the corrections in version 0.2.0.

## What the author did

The author set the problem and the written specification, chose the systems,
directed the work, and reviewed the results. The author is responsible for the
accuracy, originality and licensing of the software and the papers.

Version 0.1.0 shipped a false central claim that passed through AI-written
analysis and the author's review. An adversarial review of the released version
caught it. We record this because it bears on how much review AI-assisted
scientific work needs.
