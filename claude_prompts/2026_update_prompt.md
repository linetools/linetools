# 2026 Update

## Goals

This document carries out the fixes agreed in the Q&A of
`claude_prompts/start_up.md` (2026-09-21).  It modernizes `linetools`'
packaging, CI and Python support, removes accumulated dead weight, and repairs
several latent Python-3 bugs that were papered over rather than fixed.

`linetools` is a public, multi-author package with downstream users.  Preserve
the public API, deprecate rather than delete, and note any behaviour change so
it can go into `CHANGES.rst`.  Read `CLAUDE.md` first.  Run Python only via
`conda run -n astro`.  **Xavier performs all git commands.**

Each numbered prompt below is meant to be executed in its own session.  Do only
the task you are pointed at, then append a dated entry under `## Logs`.

## Prompts

1. Read this file.  Execute the 1st task under "CI trigger"
2. Read this file.  Execute the 1st task under "Remove dead packaging files"
3. Read this file.  Execute the 1st task under "Migrate to pyproject.toml"
4. Read this file.  Execute the 1st task under "Expose __version__"
5. Read this file.  Execute the 1st task under "Rebuild CI"
6. Read this file.  Execute the 1st task under "Root LICENSE"
7. Read this file.  Execute the 1st task under "README"
8. Read this file.  Execute the 1st task under "Morton03 and the gzip bug"
9. Read this file.  Execute the 1st task under "Dead assertions"
10. Read this file.  Execute the 1st task under "Test artifacts"
11. Read this file.  Execute the 1st task under "Astropy deprecations"

Order matters for 2 -> 3 -> 4 -> 5 (they touch the same metadata and depend on
each other).  6 through 11 are independent of one another and of 1-5.

Before starting any prompt, establish the baseline:

```
conda run -n astro pytest linetools -q
```

The green baseline as of 2026-09-21 is **231 passed, 3 skipped, 40.5s**.  After
your change, the count must not drop.  Say explicitly in the log what the new
count is and account for any difference.

---

## CI trigger

1. In `.github/workflows/ci_tests.yml`, the push trigger watches branch `main`:

   ```yaml
   on:
     push:
       branches:
       - main
   ```

   The default branch of this repository is **`master`**, so the push trigger
   has never fired -- only `pull_request` CI has ever run.

   Change `main` to `master`, and add the `2026` development branch:

   ```yaml
   on:
     push:
       branches:
       - master
       - 2026
     pull_request:
   ```

   Change nothing else in the file; prompt #5 rebuilds it properly.  This is
   deliberately a one-line fix so it can land immediately.

---

## Remove dead packaging files

1. Delete the following, all of which are provably unreachable (see the Report
   in `start_up.md` §1):

    - `old_setup.py` -- astropy-helpers era; imports `ah_bootstrap` and
      `astropy_helpers.setup_helpers`, neither of which exists anywhere on this
      machine.  It cannot run.
    - `old_setup.cfg` -- `[ah_bootstrap]`, `[build_sphinx]`, `[upload_docs]`,
      and a stale duplicate `[metadata]` block that contradicts the live
      `setup.cfg`.  Nothing reads it.
    - `.travis.yml` -- `travis-ci.org` no longer exists.

2. Repair `MANIFEST.in`, which is stale in four separate ways:

    - `include README.rst` -- the file is `README.md`.  (If prompt #7 renames
      it to `README.rst`, this line becomes correct instead; sequence
      accordingly and say in the log which way you resolved it.)
    - `include ez_setup.py` and `include ah_bootstrap.py` -- neither exists.
    - `recursive-include *.pyx *.c *.pxd` is malformed: `recursive-include`
      takes a directory as its first argument.  There are no Cython or C
      sources in this package, so delete the line rather than repairing it.
    - `recursive-include cextern *` and the four `astropy_helpers` lines --
      none of those directories exist.
    - `recursive-include scripts *` points at a top-level `scripts/` that does
      not exist; the real one is `linetools/scripts`, which is already covered
      as package data.  Delete the line.

   Keep: `include CHANGES.rst`, `include setup.cfg` (until prompt #3 removes
   it), `recursive-include docs *`, `recursive-include licenses *`, the three
   `prune` lines, and `global-exclude *.pyc *.o`.

3. Verify the sdist still contains what it should:

   ```
   conda run -n astro python -m build --sdist --outdir /tmp/lt_sdist
   ```

   If `build` is not installed in `astro`, say so and stop rather than
   installing it without asking.

---

## Migrate to pyproject.toml

1. Move all packaging metadata out of `setup.cfg` into `pyproject.toml`
   (PEP 621), then delete `setup.py` and `setup.cfg`.  This is Q&A answer #2,
   option (b).

   The translation table:

   | From `setup.cfg` | To `pyproject.toml` |
   |---|---|
   | `[metadata]` + `[options]` | `[project]` |
   | `[options.extras_require]` | `[project.optional-dependencies]` |
   | `[options.entry_points] console_scripts` | `[project.scripts]` |
   | `[options.package_data]` | `[tool.setuptools.package-data]` |
   | `[coverage:report]` | `[tool.coverage.report]` |
   | `[tool:pytest]` | `[tool.pytest.ini_options]` |

   Specific decisions to apply:

    - **`requires-python = ">=3.11"`** (Q&A #5).  Python 3.8/3.9/3.10 are
      dropped.  This is a **breaking change for downstream installers** and
      must go in `CHANGES.rst`.
    - **Classifiers**: replace the `3.8`/`3.9`/`3.10` entries with `3.11`,
      `3.12`, `3.13`, `3.14`.  Do not claim `3.15` until it is released.
    - **Drop `use_2to3 = False`** -- a removed setuptools option that newer
      setuptools versions reject outright.
    - **Drop `github_project` and `edit_on_github`** -- astropy-helpers
      leftovers with no PEP 621 equivalent.  Setuptools already ignores them.
    - Keep the dependency floors as they are for now; prompt #5 decides whether
      to raise them, and doing both at once makes the diff unreadable.

2. **`setuptools_scm`.**  `setup.py` currently calls
   `use_scm_version={'write_to': 'linetools/version.py', 'write_to_template': ...}`
   with a hard-coded fallback version baked into the template.  Replace with:

   ```toml
   [tool.setuptools_scm]
   version_file = "linetools/version.py"
   ```

   Note two things before you do:

    - `version_file` writes a **plain literal**, not the current
      try/except-around-`get_version` template.  The behaviour differs for
      sdists built outside a git checkout.  Confirm the sdist still reports a
      sane version.
    - **`setuptools_scm` is not installed in the `astro` environment.**  The
      existing `linetools/version.py` is therefore falling through to its
      hard-coded fallback (`0.3.3.dev23+g99e6dbe3c`), which currently happens to
      be right only because the file was generated at this commit.  Flag this;
      do not `pip install` into `astro` without asking.

3. `linetools/version.py` is generated and correctly gitignored (`*/version.py`).
   Leave it that way.

4. Verify: `conda run -n astro pip install -e .` still works, the `lt_*` console
   scripts are still registered, `conda run -n astro pytest linetools -q` still
   gives 231 passed, and package data (the line lists under
   `linetools/data/`) is still found.  That last one is the real risk in this
   migration -- the `[options.package_data]` glob is six levels deep and easy to
   get wrong.

---

## Expose __version__

1. `linetools/__init__.py` is completely empty, so `linetools.__version__`
   raises `AttributeError` (Q&A #7 approved adding it).

   Add a version export that does not depend on `setuptools_scm` being
   installed at runtime:

   ```python
   from importlib.metadata import PackageNotFoundError, version as _version

   try:
       __version__ = _version("linetools")
   except PackageNotFoundError:  # not installed, e.g. running from a source tree
       __version__ = "unknown"
   ```

   Prefer this over importing `linetools.version`, because `linetools/version.py`
   is a generated file that may be absent in a fresh checkout.

2. This is a public-API **addition**, so it is backwards compatible.  Add a line
   to `CHANGES.rst` anyway.

3. Verify: `conda run -n astro python -c "import linetools; print(linetools.__version__)"`.

---

## Rebuild CI

1. Rewrite `.github/workflows/ci_tests.yml` and prune `tox.ini`.  Q&A #5 sets
   the support range: **Python 3.11 and higher, with tests added for 3.14 and
   3.15.**

   Workflow changes:

    - Triggers: push to `master` and `2026`, plus `pull_request` (prompt #1
      already did this part).
    - Bump `actions/checkout@v2` -> `@v4` and `actions/setup-python@v2` -> `@v5`.
      Both v2 actions run on a Node runtime GitHub has already sunset.
    - Matrix: Python `3.11`, `3.12`, `3.13`, `3.14` on `ubuntu-latest`.  Add
      `macos-latest` on one Python version only -- Xavier develops on macOS and
      nothing currently tests it.
    - **Python 3.15**: as of 2026-09 this is not yet a final release.  It needs
      `allow-prereleases: true` on `setup-python` and should be marked
      `continue-on-error: true`.  Expect it to fail on *dependency wheels*
      (numpy/astropy/h5py) long before it fails on linetools code; that is not a
      linetools bug.  See Q&A #1 in this document.
    - Drop the `toxenv: test-alldeps` matrix entry.  **There is no `alldeps`
      extra** in the package metadata, so `test` and `test-alldeps` currently
      build the identical environment -- a third of the matrix is a duplicate of
      another third.  Either delete the entry or define a real `alldeps` extra.
    - Delete `env: SETUP_XVFB: True`.  It is an astropy-helpers variable and
      does nothing for tox.  `MPLBACKEND=agg` (already in `tox.ini`) is what
      actually matters, and no test needs a display.
    - Delete the "Install linetools requirements" step.  It pip-installs into
      the runner environment, which tox then ignores in favour of its own
      isolated env.  Pure waste.
    - Add a docs job: `sphinx-build -W -b html docs docs/_build/html`.  There is
      currently **no CI on the documentation at all**, which is how RTD
      breakage ships unnoticed.

   `tox.ini` changes:

    - `envlist` says `py{38,39,310}` while CI runs 3.11+.  They do not overlap;
      it only works because the workflow calls bare `tox -e test`, which ignores
      `envlist`.  Update to `py{311,312,313,314}`.
    - Remove the `numpy{119,120}` factors from `envlist`: there are no matching
      pins in `deps`, so those environments silently install whatever pip
      resolves and test nothing in particular.
    - Delete `[testenv:conda]`.  It references `{toxinidir}/environment.yml`,
      **which does not exist**, and requires `tox-conda`, which is unmaintained
      and incompatible with tox 4 (`astro` has tox 4.63.0).
    - The `NIGHTLY` indexserver points at `pypi.anaconda.org/scipy-wheels-nightly`,
      retired in favour of `scientific-python-nightly-wheels`.  Either update the
      URL or drop the `numpydev` factor; as written it is broken.

2. **Dependency floors.**  The declared floors (`numpy>=1.20`, `astropy>=5.2.1`,
   `scipy>=1.6`, `matplotlib>=3.3`, `IPython>=7.10.0`, `qtpy>=1.9`) are far below
   anything ever tested, and `numpy>=1.20` is not credible for a package that
   also supports numpy 2.x.  Now that `requires-python` is `>=3.11`, several
   floors are unreachable anyway (numpy 1.20 has no 3.11 wheels).

   Do **one** of the following and say which in the log:

    - Raise the floors to the oldest versions with Python 3.11 wheels
      (roughly `numpy>=1.24`, `astropy>=5.3`, `scipy>=1.10`), **or**
    - Add a `test-oldestdeps` CI job that pins to the declared floors, which is
      the only thing that would ever validate them.

   Raising them is the honest option; an untested floor is a false promise.

---

## Root LICENSE

1. The license lives at `licenses/LICENSE.rst`, so GitHub does not detect it and
   the repository shows no license in its sidebar or API metadata.  Q&A #8
   approved adding one.

   Copy `licenses/LICENSE.rst` to a root `LICENSE` file.  Use a real copy, not a
   symlink -- symlinks do not survive sdist/wheel packaging reliably and GitHub's
   license detector does not follow them.

2. Confirm the text is the BSD-3-Clause the metadata claims (`setup.cfg` says
   `license = BSD-3`; `old_setup.cfg` said `BSD`).  If they disagree, stop and
   ask rather than guessing -- relicensing is not a judgement call to make on
   someone else's behalf.

3. Leave `licenses/` in place; other files reference it.

---

## README

1. `README.md` has a `.md` extension but its body is written in
   **reStructuredText** heading syntax:

   ```
   linetools
   =========
   ...
   Development status
   ------------------
   ```

   GitHub renders that as paragraph text with stray `===` and `---` lines rather
   than as headings.  Q&A #6 said fix it.

   Convert the body to real Markdown (`# linetools`, `## Development status`,
   `## DOI`).  Prefer this over renaming to `.rst`: the badge syntax already in
   the file is Markdown, and `setup.cfg`'s commented-out `long_description`
   lines assume `README.md` too.  If you rename instead, `MANIFEST.in` and
   `licenses/README.rst` both need updating -- say which route you took.

2. Fix the badges:

    - **Travis** (`travis-ci.org/linetools/linetools.svg`) -- the host is gone.
      Replace with the GitHub Actions badge:
      `https://github.com/linetools/linetools/actions/workflows/ci_tests.yml/badge.svg`
    - **Coveralls** -- nothing in CI has uploaded coverage in years, so the badge
      reports stale data even when it renders.  Either wire up Codecov properly
      in the workflow or remove the badge.  Do not leave a badge that lies.
    - The AstroPy and Zenodo DOI badges are fine; leave them.

---

## Morton03 and the gzip bug

**Read this section before acting -- it revises Q&A answer #4.**

Q&A #4 said to delete the three permanently-skipped tests.  Investigation on
2026-09-21 showed that is the wrong call for one of them, and that the skip was
hiding a real bug.

1. **`linetools/lists/tests/test_parse_lists.py::test_morton03`** is decorated
   `@pytest.mark.skipif("sys.version_info >= (3,0)")`, i.e. disabled on every
   Python 3.  I ran its body directly under Python 3.14 (`astro`).  Results:

    - `parse.parse_morton03(orig=True)` **works** -- 3295 rows, and both
      assertions (`wrest[5] == 930.7482`, unit is Angstrom) pass.
    - `parse.mktab_morton03(do_this=True, fits=False, outfil='tmp.vo')`
      **raises `TypeError: a bytes-like object is required, not 'str'`** at
      `linetools/lists/parse.py:768`.

   The bug is a genuine Python-2 leftover:

   ```python
   with open(outfil) as src:              # text mode
       with gzip.open(outfil+'.gz', 'wb') as dst:   # binary mode
           dst.writelines(src)            # str -> binary handle: TypeError
   ```

   **The same bug exists twice**: `parse.py:493-495` (in `parse_verner96`) and
   `parse.py:766-768` (in `mktab_morton03`).

   Fix both by opening the source in binary mode: `with open(outfil, 'rb') as src:`.

   Then **remove the `skipif` from `test_morton03`** and confirm it passes.  Do
   not delete this test: `parse_morton03` is live code, called from
   `linetools/lists/linelist.py:161` and `:163` for the `ism` and `hi` line
   lists.  Deleting its only direct test to hide a fixable bug is the wrong
   trade.

   Note this in `CHANGES.rst` as a bug fix -- `mktab_morton03` and the
   `parse_verner96(write=True)` path have been broken for every Python 3 user.

2. **`linetools/guis/tests/test_guis.py`** -- the two tests gated on
   `gui_test = pytest.mark.skipif(True, reason='test requires dev suite')`.  The
   `True` is hardcoded, so these are unconditionally dead rather than
   environment-dependent.  Per Q&A #4, **delete them**, along with the now-unused
   `gui_test` marker if nothing else references it.

   Do not attempt to revive them: they need a "dev suite" of data files that is
   not in this repository, and there are no Qt bindings installed in `astro`
   (`import qtpy` raises `QtBindingsNotFoundError`).

3. Expected result: the skip count drops from 3 to 0, the pass count goes from
   231 to 232 (`test_morton03` now runs; the two GUI tests are gone).  State the
   actual numbers in the log.

---

## Dead assertions

1. **`linetools/tests/test_utils.py:149` and `:151`** contain

   ```python
   f = ltu.overlapping_chunks(chunk2, chunk1)
   assert ~f
   ```

   `f` is a Python `bool`.  `~False` is `-1` and `~True` is `-2`, and **both are
   truthy** -- so these assertions pass no matter what the function returns.
   They are dead.  I verified this directly: the call returns `False`, and
   `bool(~False)` and `bool(~True)` are both `True`.

   Change both to `assert not f` and run the suite.  Per Q&A #3, the point is to
   find out deliberately whether `overlapping_chunks` is actually correct on
   that path.  **If the test now fails, do not "fix" it by loosening the
   assertion** -- report what the function returns and stop.

2. **`linetools/utils.py:68`** has the same `~`-on-bool pattern in *library*
   code:

   ```python
   if verbose & (~sub_test):
   ```

   This emits `DeprecationWarning: Bitwise inversion '~' on bool ... will be
   removed in Python 3.16`.  It happens to behave correctly when `verbose` is
   `True` (`True & -1 == 1`, `True & -2 == 0`) but is wrong for any other truthy
   `verbose` value -- e.g. `verbose=2` prints on every line regardless of
   `sub_test`.

   Change to `if verbose and not sub_test:`.

3. Both changes together should remove 3 of the 9 `DeprecationWarning`s from the
   run.  Confirm the warning count drops and note the new total.

---

## Test artifacts

1. The test suite writes files into the source tree.  Q&A #10 said to fix this
   properly with `tmp_path`.

   **This is larger than the Report in `start_up.md` implied.**  That report
   listed only the four files visible in `git status`; in fact `.gitignore`'s
   `tmp.*` rule was hiding most of them.  The full inventory:

   Writing into the repo root (CWD):

   | Site | File |
   |---|---|
   | `linetools/tests/test_utils.py:107` | `tmp.json` |
   | `linetools/tests/test_utils.py:112` | `tmp.json.gz` |
   | `linetools/tests/test_utils.py:117` | `tmp2.json` |
   | `linetools/tests/test_init_absline.py:62` | `tmp.json` |
   | `linetools/tests/test_init_emissline.py:36` | `tmp.json` |
   | `linetools/isgm/tests/test_use_abssys.py:167` | `tmp.json` |
   | `linetools/isgm/tests/test_init_abssys.py:38` | `J081227.432-122555.56_z2.929.json` (bare `write_json()`, no path) |
   | `linetools/isgm/tests/utils.py:107` | bare `write_json()`, no path |

   Writing into `linetools/spectra/tests/files/` and
   `linetools/isgm/tests/files/` (i.e. **inside the installed package**, via
   `data_path()`):

   | Site | File |
   |---|---|
   | `linetools/spectra/tests/test_xspec_io.py:39,93,106` | `tmp.fits` |
   | `linetools/spectra/tests/test_xspec_io.py:40,56,62` | `tmp.hdf5` |
   | `linetools/spectra/tests/test_xspec_io.py:51,52` | `tmp.fits` |
   | `linetools/spectra/tests/test_xspec_io.py:54` | `tmp.ascii` |
   | `linetools/spectra/tests/test_xspec_io.py:68` | `tmp2.hdf5` |
   | `linetools/spectra/tests/test_xspec_io.py:98` | `tmp2.fits` |
   | `linetools/isgm/tests/test_use_abscomp.py:61` | `tmp.json` |

   The `data_path()` ones are the worse category: they write into the package
   directory itself, which is read-only in a normal installation and shared
   between parallel test runs.

2. Convert every site to pytest's `tmp_path` fixture.  Nothing in this
   repository uses `tmp_path` or `tmpdir` today, so establish the pattern
   cleanly:

   ```python
   def test_something(tmp_path):
       outfil = str(tmp_path / 'tmp.json')
       ltu.savejson(outfil, d, overwrite=True)
   ```

   For the two bare `write_json()` calls, pass an explicit output path rather
   than relying on the auto-generated name from the object's coordinates.  Check
   the `write_json` signature first -- if it does not accept a path, say so and
   stop rather than changing library code as a side effect of a test cleanup.

3. Afterwards, `git status` must be clean after a full test run.  Verify that
   explicitly and say so in the log.

4. Leave `.gitignore` alone.  Extending it would hide the problem again, which is
   exactly how the `data_path()` writes went unnoticed.

---

## Astropy deprecations

*This section was not part of the Q&A; it is carried over from the Report.  Skip
it if it is not wanted.*

1. `astropy.utils.isiterable` is deprecated and emits 6
   `AstropyDeprecationWarning`s per test run.  Replace with `np.iterable` at:

    - `linetools/abund/solar.py:95`
    - `linetools/analysis/absline.py:255`, `:337`
    - `linetools/guis/utils.py:223`
    - `linetools/isgm/abscomponent.py:294`, `:386`

   Also remove the unused `from astropy.utils import isiterable` import in
   `linetools/analysis/emline.py:12`.

2. `np.iterable` is a drop-in replacement for the `isiterable` semantics used at
   all six sites (it returns `True` for anything `iter()` accepts).  Confirm the
   test count is unchanged and the 6 warnings are gone.

3. Minor, while in `linetools/lists/parse.py`: the docstring of
   `grab_galaxy_linelists()` (line 802) describes `mktab_morton00` -- a
   copy-paste error.  Fix it if convenient.

---

## Q&A

1. **Python 3.15 CI.**  You asked for tests on 3.14 and 3.15.  3.14 is fine --
   it is what `astro` runs.  But 3.15 is not final until roughly October 2026,
   so that job needs `allow-prereleases: true` and will most likely fail while
   installing numpy/astropy/h5py wheels that do not exist yet for 3.15, rather
   than on anything in linetools.  I plan to add it as `continue-on-error: true`
   so it reports without blocking merges.  Good, or would you rather wait until
   3.15 ships?

   > *Answer:*

2. **`requires-python = ">=3.11"` is a breaking change** for anyone on 3.8-3.10
   who currently does `pip install linetools`.  They will silently get an older
   release instead of an error.  Do you want a final 3.8-compatible release
   tagged before this lands, or is that not worth the effort for this package?

   > *Answer:*

3. **`mktab_morton03` / `parse_verner96(write=True)` have been broken on Python
   3 for years** (the gzip `TypeError` in prompt #8).  These are "builder-only"
   functions -- they regenerate the packaged line-list data files.  Fixing them
   is easy, but it raises the question of whether the packaged
   `morton03_table2.fits.gz` and `verner96_tab1.fits.gz` are still
   byte-reproducible from the current source data.  Do you want me to check that
   the regenerated files match what is committed, or just fix the crash and stop
   there?

   > *Answer:*

4. **Dependency floors** (prompt #5, item 2).  I lean toward raising them to the
   oldest versions with Python 3.11 wheels rather than adding an `oldestdeps` CI
   job, because the current floors are demonstrably untested and `numpy>=1.20`
   is not installable on 3.11 anyway.  Confirm, or would you rather defend the
   existing floors with a CI job?

   > *Answer:*

5. **Coveralls.**  Prompt #7 says either wire up coverage properly or drop the
   badge.  Which?  Wiring it up means adding a `cov` factor to `tox.ini` and a
   Codecov upload step, and deciding whether a coverage drop should fail CI.

   > *Answer:*

6. **`setuptools_scm` is not installed in `astro`**, so the generated
   `linetools/version.py` is currently serving a hard-coded fallback version.
   May I `pip install setuptools_scm` into `astro`?  Prompt #3's verification
   step is not meaningful without it.

   > *Answer:*

## Logging

The "Logs" section will record Claude's work.  Please use the following format:

### <Date> (Short summary of the work)

<Detailed description of the work and what you learned>

...

## Logs
