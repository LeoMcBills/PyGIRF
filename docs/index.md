---
description: PyGIRF — installation, quick start, and links to the algorithm guide.
---

# PyGIRF

Python package to calculate and use the **gradient impulse response function (GIRF)** of an MRI gradient system. The installable module name is **`pygirf`** (`src/pygirf/`).

## Installation

```bash
python3 -m pip install -e .
```

From the repository root (editable install includes runtime dependencies: NumPy, SciPy, Matplotlib).

## Quick start

Run a demo:

```bash
python3 demo_pygirf.py
python3 demo2.py
```

Or with [uv](https://github.com/astral-sh/uv):

```bash
uv run demo_pygirf.py
```

Use in Python:

```python
from pygirf import GirfProvider, GirfApplier
```

## Documentation site (Zensical)

This site is built with [Zensical](https://zensical.org/docs/create-your-site/). From the repository root:

```bash
pip install zensical
```

Or install the package with the docs extra (declares Zensical in `pyproject.toml`):

```bash
pip install -e ".[docs]"
```

Then:

```bash
zensical serve
```

Open **http://localhost:8000** — the server rebuilds when you edit files under `docs/`.

Static output:

```bash
zensical build
```

HTML is written to `site/` (ignored by git by default).

## Concepts and algorithm

The full technical guide lives here:

- **[Concepts and algorithm](guide.md)** — LTI model, `GirfProvider` / `GirfApplier`, array shapes, estimation cases, prediction, checklist.

## Citations

If you use this package for research, cite the relevant paper(s):

**Gradient system characterization**  
Vannesjo, S.J., Haeberlin, M., Kasper, L., Pavan, M., Wilm, B.J., Barmet, C., Pruessmann, K.P., 2013. Gradient system characterization by impulse response measurements with a dynamic field camera. *Magn Reson Med* 69, 583–593. <https://doi.org/10.1002/mrm.24263>

**Image reconstruction based on the GIRF characterization**  
Vannesjo, S.J., Graedel, N.N., Kasper, L., Gross, S., Busch, J., Haeberlin, M., Barmet, C., Pruessmann, K.P., 2016. Image reconstruction using a gradient impulse response model for trajectory prediction. *Magn Reson Med* 76, 45–58. <https://doi.org/10.1002/mrm.25841>

**GIRF-based spiral fMRI**  
Graedel, N.N., Kasper, L., Engel, M., Nussbaum, J., Wilm, B.J., Pruessmann, K.P., Vannesjo, S.J., 2019. Feasibility of spiral fMRI based on an LTI gradient model. *bioRxiv* 805580. <https://doi.org/10.1101/805580>

**GIRF-based pre-emphasis of gradient or shim channels**  
Vannesjo, S.J., Duerst, Y., Vionnet, L., Dietrich, B.E., Pavan, M., Gross, S., Barmet, C., Pruessmann, K.P., 2017. Gradient and shim pre-emphasis by inversion of a linear time-invariant system model. *Magn Reson Med* 78, 1607–1622. <https://doi.org/10.1002/mrm.26531>

The project README in the repository root mirrors installation and citations for GitHub browsing.
