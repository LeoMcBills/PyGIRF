# PyGIRF

Gradient impulse response function (GIRF) tools in Python for MRI gradient systems. This code was ported from MATLAB which was initially written by [jvannjo](https://github.com/jvannesjo). The original MATLAB project is [MRI-gradient/GIRF](https://github.com/MRI-gradient/GIRF).

## Installation

```bash
python3 -m pip install -e .
```

Installs the **`pygirf`** package from `src/pygirf/`.

## Quick start

Run the demo workflow:

```bash
uv run demo_pygirf.py
```

Or:

```bash
python3 demo_pygirf.py
```

Use in Python:

```python
from pygirf import GirfProvider, GirfApplier
```

## Documentation

- **[Concepts and algorithm](docs/guide.md)** — how the package models the gradient system, how `GirfProvider` and `GirfApplier` work, array shapes, and the end-to-end flow.

### Zensical site (browse locally)

Install [Zensical](https://zensical.org/docs/create-your-site/) and from the repository root run:

```bash
pip install zensical
zensical serve
```

Open **http://localhost:8000** for a live preview. Build static HTML with `zensical build` (output in `site/`).

### GitHub Pages (hosted documentation)

The workflow [`.github/workflows/deploy-docs.yml`](.github/workflows/deploy-docs.yml) builds the Zensical site and deploys it with the official [GitHub Pages Actions](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site#publishing-with-a-custom-github-actions-workflow).

**One-time setup**

1. On GitHub, open the repository → **Settings** → **Pages**.
2. Under **Build and deployment**, set **Source** to **GitHub Actions** (not “Deploy from a branch”).
3. Push to **`main`** or **`master`** (or run the workflow manually under **Actions** → **Deploy documentation** → **Run workflow**).

After a successful run, the site is available at:

`https://<your-github-username>.github.io/<repository-name>/`

For example, if the repo is `github.com/octocat/PyGIRF`, the docs URL is `https://octocat.github.io/PyGIRF/`.

The workflow sets `site_url` automatically for that URL so links, search, and sitemap stay correct. For local `zensical serve`, keep `site_url` in `zensical.toml` pointed at `http://127.0.0.1:8000` (or your dev URL).

*If you use this code package for your research, please cite at least one of the following papers (depending on the use case):*  
**Gradient system characterization**  
Vannesjo, S.J., Haeberlin, M., Kasper, L., Pavan, M., Wilm, B.J., Barmet, C., Pruessmann, K.P., 2013. Gradient system characterization by impulse response measurements with a dynamic field camera. Magn Reson Med 69, 583–593. https://doi.org/10.1002/mrm.24263  
**Image reconstruction based on the GIRF characterization**  
Vannesjo, S.J., Graedel, N.N., Kasper, L., Gross, S., Busch, J., Haeberlin, M., Barmet, C., Pruessmann, K.P., 2016. Image reconstruction using a gradient impulse response model for trajectory prediction. Magn Reson Med 76, 45–58. https://doi.org/10.1002/mrm.25841  
**GIRF-based spiral fMRI**  
Graedel, N.N., Kasper, L., Engel, M., Nussbaum, J., Wilm, B.J., Pruessmann, K.P., Vannesjo, S.J., 2019. Feasibility of spiral fMRI based on an LTI gradient model. bioRxiv 805580. https://doi.org/10.1101/805580  
**GIRF-based pre-emphasis of gradient or shim channels**  
Vannesjo, S.J., Duerst, Y., Vionnet, L., Dietrich, B.E., Pavan, M., Gross, S., Barmet, C., Pruessmann, K.P., 2017. Gradient and shim pre-emphasis by inversion of a linear time-invariant system model. Magn Reson Med 78, 1607–1622. https://doi.org/10.1002/mrm.26531  
