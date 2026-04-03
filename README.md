# PyGIRF

Gradient impulse response function (GIRF) tools in Python for MRI gradient systems. This code was ported from MATLAB which was initially written by [jvannjo](https://github.com/jvannesjo). The original MATLAB project is [MRI-gradient/GIRF](https://github.com/MRI-gradient/GIRF).

## Installation

```bash
python3 -m pip install -e .
```

Installs the **`pygirf`** package from `src/pygirf/`.

After [publishing to PyPI](#publish-to-pypi-production), install with:

```bash
pip install pygirf
```

## Publish to PyPI (production)

### Option A — GitHub Actions (recommended)

The workflow [`.github/workflows/publish-pypi.yml`](.github/workflows/publish-pypi.yml) builds and uploads to PyPI using **[trusted publishing](https://docs.pypi.org/trusted-publishers/)** (OIDC). You do **not** store a long-lived PyPI token in GitHub secrets.

**One-time setup**

1. On **PyPI**, open **[your project](https://pypi.org/manage/project/pygirf/settings/publishing/)** → **Publishing** → **Add a new pending publisher** (or manage existing).
2. Choose **GitHub** as the publisher type and enter:
   - **Owner:** your GitHub user or organization  
   - **Repository name:** `PyGIRF` (or whatever the repo is called)  
   - **Workflow name:** `publish-pypi.yml`  
   - **Environment name:** leave blank unless you add a matching [GitHub Environment](https://docs.github.com/en/actions/deployment/targeting-different-environments/using-environments-for-deployment) later.
3. Save. PyPI will trust that workflow to upload **`pygirf`**.

**Each release**

1. Bump **`version`** in `pyproject.toml` (PyPI rejects duplicate versions).
2. Commit, then create and push a **tag** whose name starts with `v` and matches the release, e.g.:

   ```bash
   git tag v0.2.0
   git push origin v0.2.0
   ```

3. Watch **Actions** → **Publish to PyPI** on GitHub. Confirm on **https://pypi.org/project/pygirf/**.

### Option B — Manual upload (token)

Use an API token from **https://pypi.org/manage/account/token/**. Username for uploads is always `__token__`; the password is the token string.

```bash
pip install build twine
rm -rf dist && python -m build
twine check dist/*
TWINE_USERNAME=__token__ TWINE_PASSWORD=pypi-your-token-here twine upload dist/*
```

## Quick start

Run the demo workflow:

```bash
uv run demo_pygirf.py
```

Or:

```bash
python3 demo_pygirf.py
```

That demo builds **synthetic** calibration data (band-limited “measured” outputs + noise), **fits** a GIRF, **predicts** the response to the same nominal gradients, and prints **relative RMSE** and **correlation** vs the synthetic ground truth (strong correlation on the dominant Z-like term shows the forward model matches the data you trained on).

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
