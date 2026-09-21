# CLAUDE.md

Guidance for Claude Code working in the `linetools` repository.

## What this is

`linetools` is a public, community Python package for the analysis of 1d
astronomical spectra -- especially quasar and galaxy spectra.  It covers
absorption and emission lines, spectral line lists, abundances, continuum
fitting, and a set of `lt_*` command-line scripts and Qt GUIs.

- Docs: https://linetools.readthedocs.org
- GitHub: https://github.com/linetools/linetools
- DOI: 10.5281/zenodo.168270

## Git

**I (Xavier) will perform git commands.**

You may run read-only git: `git status`, `git diff`, `git log`, `git show`,
`git branch`.

You must **never** run anything that mutates repository state: `git commit`,
`git push`, `git reset`, `git rebase`, `git checkout`, `git switch`,
`git merge`, `git add`, `git restore`, `git stash`, `git clean`,
`git cherry-pick`, `git revert`, `git tag`, `git remote`, or similar.  Edit
files; I will stage and commit.

This is enforced in `.claude/settings.json`'s deny list, not just stated here.

- The default branch for pull requests is `master`.
- Development currently happens on the `2026` branch.

## Python environment

Use the **`astro`** conda environment for anything that runs Python:

```
conda run -n astro python ...
conda run -n astro pytest linetools -q
```

This is an **astronomy** repository, *not* one of my Oceanography
repositories -- do **not** use `ocean14` here.

**Always go through `conda run -n astro`.**  Do not invoke a bare `python`,
`python3`, `pytest`, `tox` or `pip install`: those resolve to whatever is on
`PATH` and silently bypass the environment rule.  They are set to *ask* in
`.claude/settings.json`, so reaching for one will interrupt me for approval --
which is the point.  If you genuinely need a bare interpreter (e.g. a throwaway
script that has nothing to do with linetools), say why when you ask.

`linetools` is installed into `astro` in development mode, so `import
linetools` picks up this working tree directly.  There is no need to
reinstall after editing.

## Calculations

If you do any calculation, generate it as a Python script and write it to disk
so that I can add it to the Repository.  Do not leave results only in inline
snippets or in the transcript -- they need to be committable and re-runnable.

## This is a public package

`linetools` has outside contributors and downstream users (notably `PypeIt`
and other spectroscopy packages).  Therefore:

- **Preserve the public API.**  Do not rename or re-signature public functions,
  classes, or keyword arguments without a clear reason.
- **Deprecate rather than delete.**  Leave a shim with a `DeprecationWarning`
  when something must go away.
- **Note any behaviour change** explicitly in your reply and in the log entry,
  so that I can add it to `CHANGES.rst`.
- Do not restructure the package layout.  Prefer surgical changes over
  rewrites.

## Package map

```
linetools/
  abund/         # solar abundances, ions, relative abundances
  analysis/      # absorption/emission line analysis, continuum, interpolation
  guis/          # Qt GUI widgets (xspecgui, xabssysgui, continuum fitting, ...)
  isgm/          # AbsComponent / AbsSystem / AbsSightline (IGM/CGM objects)
  lists/         # LineList and the underlying line databases
  scripts/       # lt_absline, lt_line, lt_xspec, lt_continuumfit, lt_plot,
                 #   lt_radec, lt_solabnd, lt_xabssys, lt_get_COS_LP
  spectra/       # XSpectrum1D and spectral I/O
  data/          # packaged line lists and reference data
  spectralline.py, line_utils.py, io.py, utils.py
  tests/         # top-level tests
```

Most subpackages also carry their own `tests/` directory.  Documentation
sources live in `docs/` (Sphinx, built on ReadTheDocs).

## GUIs

Code under `linetools/guis` and `linetools/scripts` depends on Qt through
`qtpy`.  **Do not attempt to launch a GUI in a session** -- there is no
display, and it will hang or crash.  Reason about GUI behaviour from the
source instead.  If a GUI change needs interactive verification, say so and
I will run it.

## Prompts and logging

`claude_prompts/` is the source of instruction for this repository.  Read the
relevant prompt doc before acting, and do the numbered task you were pointed
at -- not the whole file.

After completing a task, append a dated entry under that doc's `## Logs`
section:

```
### YYYY-MM-DD (Short summary of the work)

<what was done, what you learned about the repo, decisions and why,
 anything flagged for me to confirm, and a note that no git commands were run>
```

The "what you learned" half is not filler -- it is how repo knowledge survives
between sessions.
