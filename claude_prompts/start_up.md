# Getting started

## Goals

This repository is `linetools`, a public, community Python package for the
analysis of 1d astronomical spectra -- especially quasar and galaxy spectra
(absorption lines, emission lines, spectral line lists, continuum fitting, and
a set of `lt_*` GUIs and scripts).  It already exists, is released, is
documented on ReadTheDocs (https://linetools.readthedocs.org) and has a Zenodo
DOI; the goal here is to bring it up to my Claude working standard, modernize
its packaging and CI, and then continue development and maintenance with Claude
as a collaborator.

Unlike my `Oceanography/python` repositories, this is a *mature, multi-author,
public* codebase.  Backwards compatibility and the existing API matter.  Prefer
surgical changes over rewrites, and do not restructure the package layout.

## Prompts

1. Read this file.  Execute the 1st task under "Claude/CLAUDE.md file"
2. Read this file.  Execute the 1st task under "Claude/Settings"
3. Read this file.  Execute the 1st task under "Claude/Skills"
4. Read this file.  Execute the 1st task under "Basic start up"
5. Read this file.  Execute the 1st task under "Tests"

6. I have responded to your Q&A.  Please read and react.  If you have any
additional questions, put them in Q&A.  Use Opus 5.

## Claude

### CLAUDE.md file

1. Please generate a basic CLAUDE.md file for this project.  Have it indicate:

    - I will perform git commands.  You may run read-only git (`status`,
      `diff`, `log`, `show`, `branch`) but never `commit`, `push`, `reset`,
      `rebase`, `checkout` or anything else that mutates repository state.
    - If you do any calculation, generate it as a python script and write it to
      disk so that I can add it to the Repository.
    - If you need to run Python, use the `astro` conda environment.  This is an
      astronomy repository, *not* one of my Oceanography repositories -- do not
      use `ocean14` here.  `linetools` is already installed into `astro` in
      development mode, so `import linetools` picks up this working tree.
    - This is a public package with outside contributors.  Preserve the public
      API; deprecate rather than delete; note any behaviour change so I can
      mention it in `CHANGES.rst`.
    - The default branch for pull requests is `master`.  Development currently
      happens on the `2026` branch.
    - A short map of the package: `linetools/{abund, analysis, guis, isgm,
      lists, scripts, spectra}`, plus `spectralline.py`, `line_utils.py`,
      `io.py`, `utils.py`.  Tests live in `linetools/tests` and in per-subpackage
      `tests/` directories; docs live in `docs/`.
    - GUI code under `linetools/guis` and `linetools/scripts` depends on Qt via
      `qtpy`.  Do not attempt to launch GUIs in a session; reason about them
      from the source instead.

### Settings

1. Copy over the `settings.json` file from the `IOPtics` repository
   (`/Users/xavier/Oceanography/python/IOPtics/.claude/settings.json`) into
   `.claude/settings.json` here.  Copy the *policy*, not the accumulated
   path-specific cruft -- most of the long tail of allow entries references
   IOPtics' own scratchpad paths and test files, and none of it applies here.

   Then adapt it for this repository:

    - Replace the `conda run -n ocean14:*` allow entry with `conda run -n astro:*`.
    - Keep the deny list as-is: `sudo`, `rm -rf /`, `rm -rf ~`, `git push`,
      `git commit`, `git reset`, `git rebase`.  Add `git checkout` and
      `git merge` to the deny list.
    - Keep `rm:*` under `ask`.
    - Replace the oceanography/optics publisher domains in the `WebFetch`
      allow list with astronomy ones: `adsabs.harvard.edu`, `ui.adsabs.harvard.edu`,
      `arxiv.org`, `iopscience.iop.org`, `academic.oup.com`, `aanda.org`,
      `docs.astropy.org`, `linetools.readthedocs.io`.  Keep `doi.org` and the
      Crossref lookups -- literature checking is a first-class activity here too.
    - Allow `pytest:*` and `tox:*`.

### Skills

1. Copy over the `skills/` files from the `IOPtics` repository
   (`/Users/xavier/Oceanography/python/IOPtics/.claude/skills/`) into
   `.claude/skills/` here.  That is `critical-partner` and `grill-me`; copy
   each `SKILL.md` verbatim.

## Basic start up

1. This repository is *not* a fresh one -- the basic Python-repository files
   already exist, but several are stale or duplicated.  Please explore and
   write a Report below.  **Do no editing yet.**

   Specifically, look at and report on:

    - The packaging situation.  There is a nearly empty `pyproject.toml`
      (build-system only), a populated `setup.cfg` carrying all the metadata, a
      thin `setup.py`, and leftover `old_setup.py` / `old_setup.cfg`.  What is
      actually used at build time?  What can be deleted, and what would it take
      to move the metadata into `pyproject.toml`?
    - CI.  `.travis.yml` is present and Travis is dead for this project; there
      is a `.github/` directory and a `tox.ini`.  What runs today, what is
      vestigial, and what would a minimal GitHub Actions matrix look like
      (python versions, Qt/headless concerns, doc build)?
    - The README badges -- Travis and Coveralls -- which are likely broken.
    - Dependency floors in `setup.cfg` (`numpy>=1.20`, `astropy>=5.2.1`, ...)
      versus what is actually installed in the `astro` environment, and whether
      anything in the code still assumes older behaviour (e.g. the `np.int`
      fixes in recent commits suggest more numpy-2 era cleanup may be pending).
    - Anything else that diverges from how I organize my other repositories
      (look at `~/Oceanography/python/{IOPtics, PAB, cugn}` for the house
      layout), and which of their `settings.json` conventions are worth pulling
      in beyond what prompt #2 already covers.

   Write findings into the `### Report` section below.  Do not touch code.

### Report

*Written 2026-09-21.  Exploration only -- no code was edited.*

#### 1. Packaging

**What is actually used at build time.**  The live chain is `pyproject.toml`
(build-system: `setuptools` + `setuptools_scm`, backend `setuptools.build_meta`)
-> `setup.py` -> `setup.cfg`.  `setup.py` is a thin shim that does exactly two
things: it intercepts the legacy `python setup.py test` invocation with a help
message, and it calls `setup(use_scm_version={'write_to':
'linetools/version.py', ...})`.  All real metadata -- name, deps, extras, entry
points, package data, coverage config and `[tool:pytest]` -- lives in
`setup.cfg`.

**Dead weight, in order of confidence.**

| File | Status | Notes |
|---|---|---|
| `old_setup.py` | **Dead.** Safe to delete. | Pure `astropy_helpers` era: imports `ah_bootstrap` and `astropy_helpers.setup_helpers`, sets `builtins._ASTROPY_SETUP_`, branches on `sys.version_info[0] >= 3` for `ConfigParser`. None of those modules exist in the repo or in `astro`. It cannot run. |
| `old_setup.cfg` | **Dead.** Safe to delete. | `[ah_bootstrap]`, `[build_sphinx]`, `[upload_docs]`, plus a duplicate stale `[metadata]` block (says `license = BSD`, `edit_on_github = True`, `package_name =` -- all contradicted by the live `setup.cfg`). Nothing reads it. |
| `.travis.yml` | **Dead.** Safe to delete. | See CI below. |
| `MANIFEST.in` | **Stale but live.** Needs editing, not deleting. | References `README.rst` (the file is `README.md`), `ez_setup.py`, `ah_bootstrap.py`, `cextern/`, `astropy_helpers/` -- none of which exist. `recursive-include *.pyx *.c *.pxd` is malformed (missing the directory argument). `recursive-include scripts *` points at a non-existent top-level `scripts/`; the real one is `linetools/scripts`. |

**Moving metadata into `pyproject.toml`.**  Mechanically straightforward --
setuptools has read `[project]` from `pyproject.toml` since v61, and everything
in `setup.cfg`'s `[metadata]`/`[options]` has a direct PEP 621 equivalent:

- `[metadata]` + `[options]` -> `[project]` (`name`, `description`, `authors`,
  `license`, `requires-python`, `classifiers`, `dependencies`, `urls`)
- `[options.extras_require]` -> `[project.optional-dependencies]`
- `[options.entry_points] console_scripts` -> `[project.scripts]`
- `[options.package_data]` -> `[tool.setuptools.package-data]`
- `[coverage:report]` -> `[tool.coverage.report]`
- `[tool:pytest]` -> `[tool.pytest.ini_options]`

Three things make it more than a copy-paste:

1. **`setuptools_scm` version handling.**  `setup.py`'s `use_scm_version` with a
   custom `write_to_template` becomes `[tool.setuptools_scm]` with
   `version_file = "linetools/version.py"`.  The current template embeds a
   hard-coded fallback version; the modern `version_file` mechanism writes a
   plain literal instead.  Both work, but they differ for sdists built outside a
   git checkout, so this needs a deliberate choice.
2. **`github_project` / `edit_on_github`** in `[metadata]` are astropy-helpers
   leftovers.  Setuptools ignores them; PEP 621 has no home for them.  Drop.
3. **`use_2to3 = False`** in `[options]` is a removed setuptools option; recent
   setuptools errors on its mere presence in some configurations.  It should go
   regardless of whether we migrate.

If we migrate, `setup.py` can be deleted entirely (the legacy-`test` help
message is not worth a file) and `setup.cfg` can go with it.

**One packaging observation worth flagging:** `setuptools_scm` is **not
installed** in `astro`, so `linetools/version.py`'s `try:` block falls through
to its hard-coded fallback `'0.3.3.dev23+g99e6dbe3c'`.  That happens to match
what the installed distribution reports, but only because the file was generated
at the current commit.  It will silently go stale.

#### 2. CI

**`.travis.yml` is vestigial.**  Travis stopped serving open-source projects on
`travis-ci.org` years ago; the file is `language: c` with a conda bootstrap and
an astropy-helpers-era matrix.  Nothing runs it.  Delete.

**`.github/workflows/ci_tests.yml` is the live CI, and it has a real bug:**

```yaml
on:
  push:
    branches:
    - main
```

**The default branch of this repository is `master`, not `main`.**  The push
trigger has therefore never fired; only the `pull_request` trigger works.  This
is a one-word fix and probably the single highest-value item in this report.

Other issues in the workflow:

- `actions/checkout@v2` and `actions/setup-python@v2` are both deprecated and
  run on a Node version GitHub has sunset.  They warn today and will eventually
  hard-fail.  Bump to `@v4`/`@v5`.
- `SETUP_XVFB: True` is an astropy-helpers-era variable.  It does nothing for
  tox.  The modern equivalent is `MPLBACKEND=agg` (which `tox.ini` already
  sets), plus an explicit `xvfb-run` only if Qt tests are ever enabled.
- The "Install linetools requirements" step pip-installs deps into the *runner*
  environment, but tox then builds its own isolated env.  That step is wasted
  work.

**`tox.ini` is partly stale:**

- `envlist` is `py{38,39,310}` while the Actions matrix runs 3.11/3.12/3.13.
  These do not overlap.  It works only because the workflow calls bare
  `tox -e test`, which ignores `envlist` and uses the ambient Python.
- `extras = test, alldeps` -- **there is no `alldeps` extra** in `setup.cfg`
  (the extras are `pyside2`, `pyqt5`, `test`, `docs`, `dev`).  pip warns and
  continues, which means `test` and `test-alldeps` in the CI matrix run the
  *identical* environment.  One third of the matrix is redundant.
- `deps` defines pins for `numpy123`/`numpy124`/`astropylts`/`numpydev`/
  `astropydev`, but `envlist` references `numpy{119,120}` -- factors with no
  matching pin, so they silently install whatever pip resolves.
- `[testenv:conda]` references `{toxinidir}/environment.yml`, **which does not
  exist**, and requires `tox-conda`, which is unmaintained and incompatible with
  tox 4 (4.63.0 is what `astro` has).
- The `NIGHTLY` indexserver points at `pypi.anaconda.org/scipy-wheels-nightly`,
  retired in favour of `scientific-python-nightly-wheels`.  `numpydev` is
  therefore broken.

**A minimal, honest GitHub Actions matrix** would be:

- Trigger on `push` to `master` **and** `pull_request`.
- `ubuntu-latest` plus `macos-latest` (you develop on macOS; nothing currently
  tests it).
- Python 3.11, 3.12, 3.13 -- and a decision about 3.14, which is what `astro`
  runs locally, so local and CI currently disagree.
- Jobs: `test` (core deps), `test-devdeps` (astropy/numpy dev, allowed to
  fail), `test-oldestdeps` (pinned to the `setup.cfg` floors -- the only way
  those floors ever get validated), and a docs build.
- Drop `test-alldeps` until an `alldeps` extra actually exists, or define one.
- Qt/headless: the GUI tests are **unconditionally skipped** (see §5), so no Qt
  is needed in CI today.  `MPLBACKEND=agg` covers matplotlib.  If the GUI tests
  are ever turned on, that job needs `xvfb` *and* a real binding
  (`pyqt5`/`pyside6`) -- `qtpy` alone is not enough.
- A docs job (`sphinx-build -W -b html docs docs/_build/html`) would catch
  ReadTheDocs breakage before it ships.  There is no CI on docs at all today.

#### 3. README badges

Both status badges are broken:

- **Travis** -- `travis-ci.org/linetools/linetools.svg?branch=master`.  The
  whole `travis-ci.org` host is gone.  Dead image.
- **Coveralls** -- nothing in the current CI uploads coverage (`tox.ini` has no
  `cov` commands wired to an upload step, and the workflow never calls
  codecov/coveralls), so even when the badge renders it reports data from years
  ago.

Suggested replacements: a GitHub Actions status badge
(`github.com/linetools/linetools/actions/workflows/ci_tests.yml/badge.svg`), and
either wire up Codecov properly or drop the coverage badge.  The AstroPy and
Zenodo DOI badges are fine.

Minor but visible: the README is `README.md` but its body is written in
**reStructuredText heading style** (`linetools` / `=========`,
`Development status` / `------------------`).  GitHub renders that as plain text
with stray `===` lines rather than as headings.  It should be converted to real
Markdown or renamed to `README.rst`.  Note `MANIFEST.in` and
`licenses/README.rst` both still assume `.rst`.

#### 4. Dependency floors vs. reality

Installed in `astro` (Python **3.14.6**):

| Package | Floor in `setup.cfg` | Installed in `astro` |
|---|---|---|
| numpy | `>=1.20` | **2.5.2** |
| astropy | `>=5.2.1` | **8.0.1** |
| scipy | `>=1.6` | **1.18.0** |
| h5py | `>=3.7` | 3.16.0 |
| matplotlib | `>=3.3` | 3.11.1 |
| PyYAML | `>=5.1` | 6.0.3 |
| IPython | `>=7.10.0` | 9.17.1 |
| qtpy | `>=1.9` | 2.4.3 |
| importlib_resources | `>=5.7` | 7.1.0 |

The floors sit far below what is installed and are almost certainly untested --
nothing in CI pins to them.  `numpy>=1.20` in particular is not credible
alongside numpy 2.x support; numpy 2 removed enough public API that code working
on 2.5 has no reason to work on 1.20.  Either raise the floors to something we
will actually test (`numpy>=1.24`, `astropy>=6.0`) or add an `oldestdeps` CI job
that defends them.

Also stale: `python_requires = >=3.8`, and the classifier list stops at Python
3.10 while CI runs 3.11-3.13 and local runs 3.14.  Python 3.8 and 3.9 are both
end-of-life.

**numpy-2 cleanup status.**  Good news: grepping for the removed aliases
(`np.int`, `np.float`, `np.bool`, `np.str`, `np.object`, `np.complex`,
`np.long`, `np.unicode`) and the other numpy-2 removals (`np.alltrue`,
`np.sometrue`, `np.product`, `np.cumproduct`, `np.NaN`, `np.Inf`, `np.float_`,
`np.in1d`, `np.row_stack`) returns **zero hits** across `linetools/`.  The
`np.int` -> `int` commit (`f43aa9a`) appears to have finished that job, and the
suite passes cleanly on numpy 2.5.2.

**What is still pending is Python-3.16 and astropy deprecation cleanup, not
numpy:**

1. **`linetools/utils.py:68` -- `if verbose & (~sub_test):`.**  `sub_test` is a
   bool; `~False` is `-1` and `~True` is `-2`, both truthy.  This produces the
   `DeprecationWarning: Bitwise inversion '~' on bool ... removed in Python
   3.16`.  It happens to behave correctly for `verbose=True` (because
   `True & -1 == 1` and `True & -2 == 0`) but is wrong for any other truthy
   `verbose`.  Should be `if verbose and not sub_test:`.

2. **`linetools/tests/test_utils.py:149` and `:151` -- `assert ~f`.**  This one
   is worse: it is a **dead assertion**.  I confirmed that
   `ltu.overlapping_chunks(np.array([5,7]), np.array([1,3]))` returns the Python
   bool `False`, that `~False == -1` (truthy), and that `~True == -2` is *also*
   truthy.  **Both assertions pass regardless of what the function returns.**
   They should be `assert not f`.  Fixing them is the only way to learn whether
   `overlapping_chunks` is actually correct on that path.

   **Follow-up 2026-09-21 (prompt #6):** I ran the body of the permanently
   skipped `test_morton03` directly under Python 3.14.  `parse_morton03(orig=True)`
   passes both its assertions, but `mktab_morton03(do_this=True, fits=False,
   outfil=...)` raises `TypeError: a bytes-like object is required, not 'str'`
   at `lists/parse.py:768` -- a real Python-2 leftover (text-mode `open` feeding
   a binary `gzip` handle).  The identical bug sits at `lists/parse.py:493-495`
   in `parse_verner96`.  So the `skipif` was not just stale; it was concealing a
   live bug in code that `lists/linelist.py:161,163` depends on.

3. **`astropy.utils.isiterable` is deprecated** (6 `AstropyDeprecationWarning`s
   per run) at `abund/solar.py:95`, `analysis/absline.py:255` and `:337`,
   `guis/utils.py:223`, `isgm/abscomponent.py:294` and `:386`; also imported but
   unused in `analysis/emline.py`.  The astropy-recommended replacement is
   `np.iterable`.  Low-risk and mechanical.

#### 5. Divergences from the house layout, and other findings

Compared against `~/Oceanography/python/{IOPtics, PAB, cugn}`:

| House convention | linetools |
|---|---|
| `requirements.txt` | **Absent.**  Deps live only in `setup.cfg` extras. |
| `setup.py` | Present (shim). |
| `CLAUDE.md` | **Was absent** -- created by prompt #1 this session. |
| `claude_prompts/` | **Was absent** -- created for this start-up. |
| `.claude/{settings.json,skills}` | **Was absent** -- created by prompts #2/#3. |
| `LICENSE` at repo root | **Absent.**  It lives at `licenses/LICENSE.rst`.  Every sibling has a root `LICENSE`, and GitHub will not detect a license from `licenses/`. |
| `pytest.ini` (PAB) | Absent; `[tool:pytest]` in `setup.cfg` covers it.  Fine. |
| `docs/` | Present, and far richer than any sibling (real Sphinx + RTD). |

Two further findings outside the prompt's list that matter:

1. **The test suite writes into the source tree.**  After a run, `git status`
   shows four new untracked files: `J081227.432-122555.56_z2.929.json` and
   `tmp2.json` at the repo root, and `tmp2.fits` and `tmp2.hdf5` inside
   `linetools/spectra/tests/files/`.  `.gitignore` covers `tmp.*` but not
   `tmp2.*`, and nothing covers the `J08...json`.  Tests should write to
   `tmp_path`; failing that, `.gitignore` needs updating.  This is exactly how
   stray files end up committed.

   **Revised 2026-09-21 (prompt #6):** those four files are only what `git
   status` reveals -- `.gitignore`'s `tmp.*` rule was hiding the rest.  The real
   count is roughly twenty write sites across six test files, and several of
   them write *inside the package directory* via `data_path()` rather than into
   the repo root.  The full inventory is in
   `claude_prompts/2026_update_prompt.md` under "Test artifacts".

2. **`linetools/__init__.py` is completely empty**, so `linetools.__version__`
   raises `AttributeError`.  The version is reachable only via
   `importlib.metadata.version('linetools')` or `linetools.version.version`.
   Nearly every astronomy package exposes `__version__`, and downstream code
   will expect it.

**On `settings.json` conventions worth pulling in beyond prompt #2:** the
IOPtics allow-list has accreted roughly twenty hyper-specific entries pinned to
IOPtics scratchpad paths and its conda prefix; I dropped all of them, as the
prompt asked.  Two parts of the IOPtics *policy* I did **not** change and that
you may want to reconsider are raised in Q&A #9.

## Tests

1. Run the test suite in the `astro` environment and report the result:

   ```
   conda run -n astro pytest linetools -q
   ```

   Report how long it takes, which tests are skipped and why, whether anything
   fails, and whether any test requires network access, a display, or data
   files that are not in the repository.  Do not fix anything yet -- just
   characterize the baseline so we know what "green" looks like before making
   changes.  Write the result into the `## Logs` entry for this task.

## Q&A

Questions from the 2026-09-21 start-up session, in rough priority order.

1. **The `main` vs `master` CI bug -- fix it now?**  One word in
   `.github/workflows/ci_tests.yml`: the push trigger watches `main` but the
   default branch is `master`, so push CI has never run.  Highest value, lowest
   risk item in the report.  Unless you are mid-migration to a `main` default
   branch, in which case renaming the branch is the fix instead.

   > *Answer:*  We are going to do fixes like this in a separate prompt doc named `2026_update_prompt.md`.  You can start constructing that.

2. **Packaging: consolidate into `pyproject.toml`, or just delete the dead
   files?**  Two very different scopes.  (a) Minimal: delete `old_setup.py`,
   `old_setup.cfg`, `.travis.yml`; repair `MANIFEST.in`.  (b) Full: also migrate
   `setup.cfg` metadata into `pyproject.toml` and drop `setup.py` and
   `setup.cfg`.  (b) is a large diff on a public package and changes how sdists
   get their version.  I lean (a) now, (b) as its own staged prompt doc.

   > *Answer:* We will do (b).  Add to the new prompt doc.

3. **The two dead test assertions (`assert ~f`) -- fix them and see what
   breaks?**  Changing them to `assert not f` may reveal that
   `overlapping_chunks` has been wrong on that path for years, or may be a
   no-op.  I would rather find out deliberately than leave it.

   > *Answer:*  Same as above.

4. **The three permanently-skipped tests.**  `test_morton03` is disabled on all
   Python 3 (`skipif("sys.version_info >= (3,0)")`) and the two GUI tests are
   `skipif(True)`.  Delete them, fix them, or convert the GUI ones to a real
   conditional skip (`skipif(no Qt bindings)`) and install a Qt binding in
   `astro`?

   > *Answer:*  Delete them; in the prompt doc

5. **Python and dependency floors.**  `python_requires = >=3.8`, classifiers
   stopping at 3.10, CI on 3.11-3.13, and you on 3.14 locally.  What should the
   real supported range be?  And should I raise the numpy/astropy floors to
   something CI actually tests, or add an `oldestdeps` job to defend the current
   ones?

   > *Answer:*  Let's require 3.11 or higher.  And add tests for 3.14 and 3.15

6. **`README.md` is written in reStructuredText** and renders with stray `===`
   lines on GitHub.  Convert the body to Markdown, or rename to `README.rst`?
   (`MANIFEST.in` and `licenses/README.rst` both assume `.rst`.)

   > *Answer:* Fix; new prompt doc

7. **Should `linetools/__init__.py` export `__version__`?**  It is empty today,
   so `linetools.__version__` raises `AttributeError`.  This is a public-API
   *addition*, so it is safe -- but it is a deliberate choice.

   > *Answer:*  This is ok.

8. **Root `LICENSE` file.**  The license sits at `licenses/LICENSE.rst`, so
   GitHub does not detect it.  Add a root `LICENSE`?

   > *Answer:*  Yes, add a root LICENSE

9. **Should the permission deny-list also block `git add`, `git stash` and
   `git clean`?**  You asked me to add only `checkout` and `merge`, and I did
   exactly that.  But `CLAUDE.md` forbids `git add` in prose while
   `settings.json` does not enforce it -- and `Bash(python:*)` / `Bash(python3:*)`
   are allowed, which lets a session bypass the `astro`-only rule with the system
   interpreter.  Tighten both, or leave the prose to do the work?

   > *Answer:*  Tighten both

10. **Test artifacts in the working tree.**  Tests write `tmp2.*` and
    `J081227.432-122555.56_z2.929.json` into the source tree.  Fix properly (move
    those tests to `tmp_path`) or just extend `.gitignore`?  The former is right
    but touches more test files.

    > *Answer:*  Yes, move to tmp_path

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs

### 2026-09-21 (Execute start-up prompts 1-5: CLAUDE.md, settings, skills, packaging/CI report, test baseline)

Ran all five numbered prompts from this file in one session.

**Prompt 1 -- `CLAUDE.md`.**  Created `CLAUDE.md` at the repo root with the
directives listed under "Claude/CLAUDE.md file": git is Xavier's (read-only git
permitted, everything mutating forbidden), calculations become scripts on disk,
Python runs in the **`astro`** conda environment (explicitly *not* `ocean14` --
this is the astronomy repo), public-API preservation rules for a
multi-contributor package, the `master`/`2026` branch facts, a package map, a
do-not-launch-GUIs rule, and the `claude_prompts/` + logging convention.

**Prompt 2 -- `.claude/settings.json`.**  Copied the *policy* from
`IOPtics/.claude/settings.json` and dropped the roughly twenty accreted allow
entries pinned to IOPtics scratchpad paths, its conda prefix and its own test
files -- none of it transfers.  Adapted as instructed: `conda run -n astro:*`
replaces `conda run -n ocean14:*`; the oceanography/optics WebFetch domains are
replaced with ADS, arXiv, IOP, OUP, A&A, astropy docs and the linetools RTD
site; `doi.org` and generic Crossref lookups kept; `pytest:*` and `tox:*`
allowed; `git checkout` and `git merge` added to deny.  Validated as JSON.

**Prompt 3 -- skills.**  Copied `critical-partner/SKILL.md` and
`grill-me/SKILL.md` from IOPtics into `.claude/skills/`; verified byte identity
with `diff -rq`.

**Prompt 4 -- packaging/CI exploration.**  Written up in full under `### Report`
above.  No code was edited.  Headline findings: the GitHub Actions workflow
triggers on push to `main` while the default branch is `master`, so push CI has
never run; `old_setup.py`, `old_setup.cfg` and `.travis.yml` are provably dead;
`MANIFEST.in` is stale in four separate ways; `tox.ini` references a
non-existent `alldeps` extra (making a third of the CI matrix a duplicate) and a
non-existent `environment.yml`; the dependency floors are far below anything
tested; `README.md` is written in reStructuredText and renders wrong on GitHub;
there is no root `LICENSE`; and `linetools/__init__.py` is empty, so
`linetools.__version__` does not exist.

**Prompt 5 -- test baseline.**  `conda run -n astro pytest linetools -q`:

```
231 passed, 3 skipped, 453 warnings in 40.54s
```

Exit code 0.  **This is the green baseline.**

The three skips:

1. `linetools/guis/tests/test_guis.py:106` -- "test requires dev suite"
2. `linetools/guis/tests/test_guis.py:112` -- "test requires dev suite"
3. `linetools/lists/tests/test_parse_lists.py:15` -- `sys.version_info >= (3,0)`

Skips 1 and 2 come from `gui_test = pytest.mark.skipif(True, reason=...)` --
**hardcoded `True`**, so they are unconditionally dead rather than conditional
on the environment.  Skip 3 is `@pytest.mark.skipif("sys.version_info >= (3,0)")`,
i.e. disabled on every Python 3 -- `test_morton03` has not run in years.  All
three are dead code that looks like coverage.

**Network / display / data requirements: none.**  No test touches the network.
The only `urlopen` calls in the package are in `lists/parse.py` (lines 822-837),
inside `grab_galaxy_linelists()`, which returns early unless called with
`do_this=True`; `test_galaxy_lines` calls it with the default, so no request is
made.  (**Corrected 2026-09-21:** an earlier draft of this entry attributed
those calls to `parse_morton03(orig=True)` and the skipped `test_morton03`.
That was wrong -- `parse_morton03` reads a packaged data file and never touches
the network.  The conclusion, that the suite needs no network, is unchanged.)  No test needs a display: `qtpy` is
imported at module scope only in `guis/` and `scripts/` modules, and the single
import inside `spectra/xspectrum1d.py` is function-local (line 674, in a
plotting path).  Every data file the tests need is in the repository under the
various `tests/files/` directories.

Two caveats on that baseline worth recording:

- **There are no Qt bindings installed in `astro`** -- no PyQt5, PyQt6, PySide2
  or PySide6, so `import qtpy` raises `QtBindingsNotFoundError`.  The suite is
  green anyway, because the GUI tests that would need bindings are exactly the
  hardcoded-skip ones.  "231 passed" therefore says nothing about GUI code.
- **453 warnings**, by category: 49 `UserWarning`, 9 `DeprecationWarning`, 8
  `RuntimeWarning`, 6 `AstropyDeprecationWarning` (distinct source locations;
  453 is the raw occurrence count).  The actionable ones are detailed in §4 of
  the Report -- notably `assert ~f` at `test_utils.py:149` and `:151`, which I
  confirmed is a **dead assertion that passes regardless of the function's
  return value**, since both `~True` and `~False` are truthy.

**The suite leaves the working tree dirty.**  After a run, `git status` shows
four new untracked files: `J081227.432-122555.56_z2.929.json` and `tmp2.json` at
the repo root, plus `tmp2.fits` and `tmp2.hdf5` in
`linetools/spectra/tests/files/`.  `.gitignore` covers `tmp.*` but not `tmp2.*`.
I left them in place rather than deleting them, since deletion is a mutation you
did not ask for -- but they should not be committed.

**No git commands were run** beyond read-only `git status`, `git ls-files` and
`git log`.

**Files created this session:** `CLAUDE.md`, `.claude/settings.json`,
`.claude/skills/critical-partner/SKILL.md`, `.claude/skills/grill-me/SKILL.md`.
No pre-existing file was modified except this prompt doc.

### 2026-09-21 (Prompt 6: react to Q&A answers; create 2026_update_prompt.md)

Read the ten Q&A answers and acted on them.  Eight of the ten were deferred to a
new staged prompt doc, one was approved as-is, and one was applied immediately.

**Created `claude_prompts/2026_update_prompt.md`** -- eleven numbered prompts
covering Q&A #1-#8 and #10, in dependency order.  Prompts 2 -> 3 -> 4 -> 5
(packaging cleanup -> `pyproject.toml` migration -> `__version__` -> CI rebuild)
are sequential because they touch the same metadata; the rest are independent.
Each prompt names specific files and line numbers and states what the pass count
should be afterwards, so a fresh session can execute one without re-deriving the
analysis.  The doc carries its own `## Q&A` with six new questions.

**Applied immediately: Q&A #9 ("tighten both").**  This was the one answer not
about repository code, so deferring it to the staged doc made no sense -- it
governs how sessions operate.

- `.claude/settings.json` deny list grew from 9 to 18 entries.  Added `git add`,
  `git stash`, `git clean`, `git restore`, `git switch`, `git cherry-pick`,
  `git revert`, `git tag`, `git remote` alongside the existing
  `commit`/`push`/`reset`/`rebase`/`checkout`/`merge`.  `CLAUDE.md` had forbidden
  several of these in prose while `settings.json` permitted them.
- Moved `python:*`, `python3:*`, `pytest:*`, `tox:*` and `pip install:*` out of
  `allow` and into `ask`.  I chose *ask* rather than *deny* deliberately: a hard
  deny would block legitimate one-off scripting that has nothing to do with
  linetools, whereas *ask* makes every bypass of the `astro` rule visible and
  deliberate.  `conda run -n astro:*` stays in `allow`, so the correct path is
  the frictionless one.  Say the word if you want these denied outright.
- `CLAUDE.md` updated to match, and to state that the git rule is enforced by
  the permission system rather than by prose alone.

**Q&A #7** ("This is ok") is treated as approval to add `linetools.__version__`;
it became prompt #4 of the new doc rather than an immediate edit, since it
interacts with the `setuptools_scm` change in prompt #3.

**Pushback recorded, rather than silently accepted:**

1. **Q&A #4 is partly wrong, and I revised it in the new doc.**  It said to
   delete all three permanently-skipped tests.  I ran `test_morton03`'s body
   under Python 3.14 first.  `parse_morton03(orig=True)` passes; but
   `mktab_morton03(do_this=True, fits=False, ...)` raises `TypeError: a
   bytes-like object is required, not 'str'` at `lists/parse.py:768`, because
   a text-mode `open()` is being fed to a binary `gzip` handle.  **The same bug
   is at `lists/parse.py:493-495` in `parse_verner96`.**  So the py2-era
   `skipif` was hiding a live Python-3 bug, in a module that
   `lists/linelist.py:161,163` calls for the `ism` and `hi` line lists.  The new
   doc therefore says: fix both gzip sites, un-skip `test_morton03`, and delete
   only the two hardcoded-`skipif(True)` GUI tests.  Deleting the morton03 test
   would have buried the bug.

2. **Q&A #5's "add tests for 3.15" is not yet actionable.**  Python 3.15 is not
   final until around October 2026, so that job needs `allow-prereleases: true`
   and will almost certainly fail installing numpy/astropy/h5py wheels that do
   not exist for 3.15 yet -- a dependency-availability failure, not a linetools
   failure.  The new doc adds it as `continue-on-error: true` and raises the
   choice as Q&A #1 there.

3. **My own §5 finding understated the test-artifact problem.**  I reported four
   stray files, which is what `git status` showed -- but `.gitignore`'s `tmp.*`
   rule was masking the rest.  The real figure is about twenty write sites
   across six test files, and several write *into the package directory* via
   `data_path()`, not just into the repo root.  Corrected in §5 above and
   inventoried in full in the new doc.

**Correction to the previous log entry.**  I had written that the package's only
`urlopen` calls were reachable via `parse_morton03(orig=True)` behind the skipped
`test_morton03`.  That is wrong: they are in `grab_galaxy_linelists()`, guarded
by `do_this=False`, and `parse_morton03` never touches the network.  The
conclusion -- no test needs network access -- is unchanged.  Fixed in place above.

**No repository code was edited.**  The only files changed are `CLAUDE.md`,
`.claude/settings.json`, this prompt doc, and the new
`claude_prompts/2026_update_prompt.md`.  The probe script used to test
`test_morton03` was written to the session scratchpad, not the repo, since it
was a throwaway diagnostic rather than a calculation worth committing.

**No git commands were run** beyond read-only `git status` and `git ls-files`.
