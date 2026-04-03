# PyGIRF

Gradient impulse response function (GIRF) tools in Python for MRI gradient systems. This code was ported from MATLAB which was initially written by [jvannjo](https://github.com/jvannesjo). The original MATLAB project is [MRI-gradient/GIRF](https://github.com/MRI-gradient/GIRF).

## Installation

```bash
pip install pygirf
```

Requires Python 3.10+; dependencies are installed automatically (NumPy, SciPy, Matplotlib).

## Quick start

Run the demo:

```bash
python3 demo_pygirf.py
```

(If you use [uv](https://github.com/astral-sh/uv): `uv run demo_pygirf.py`.)

The demo uses **synthetic** calibration-style data, fits a GIRF, runs a forward prediction, and prints **relative RMSE** and **correlation** so you can sanity-check the pipeline.

Use in your own code:

```python
from pygirf import GirfProvider, GirfApplier
```

## Documentation

- **[Concepts and algorithm](docs/guide.md)** — LTI model, `GirfProvider` / `GirfApplier`, array shapes, estimation, prediction, and a short checklist.

## Citations

If you use this package for your research, please cite at least one of the following papers (depending on the use case):

**Gradient system characterization**  
Vannesjo, S.J., Haeberlin, M., Kasper, L., Pavan, M., Wilm, B.J., Barmet, C., Pruessmann, K.P., 2013. Gradient system characterization by impulse response measurements with a dynamic field camera. Magn Reson Med 69, 583–593. https://doi.org/10.1002/mrm.24263  

**Image reconstruction based on the GIRF characterization**  
Vannesjo, S.J., Graedel, N.N., Kasper, L., Gross, S., Busch, J., Haeberlin, M., Barmet, C., Pruessmann, K.P., 2016. Image reconstruction using a gradient impulse response model for trajectory prediction. Magn Reson Med 76, 45–58. https://doi.org/10.1002/mrm.25841  

**GIRF-based spiral fMRI**  
Graedel, N.N., Kasper, L., Engel, M., Nussbaum, J., Wilm, B.J., Pruessmann, K.P., Vannesjo, S.J., 2019. Feasibility of spiral fMRI based on an LTI gradient model. bioRxiv 805580. https://doi.org/10.1101/805580  

**GIRF-based pre-emphasis of gradient or shim channels**  
Vannesjo, S.J., Duerst, Y., Vionnet, L., Dietrich, B.E., Pavan, M., Gross, S., Barmet, C., Pruessmann, K.P., 2017. Gradient and shim pre-emphasis by inversion of a linear time-invariant system model. Magn Reson Med 78, 1607–1622. https://doi.org/10.1002/mrm.26531  
