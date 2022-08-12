# *infresnel*

*infresnel* is a Python package which facilitates computation of the
[Fresnel number](https://en.wikipedia.org/wiki/Fresnel_number) for
[infrasound](https://en.wikipedia.org/wiki/Infrasound) applications.

The dimensionless Fresnel number $N$ for wavelength $\lambda$ is given by
$$N = \frac{(R_\mathrm{t} - R_\mathrm{d})}{\lambda / 2},$$ where $R_\mathrm{t}$ is the
length of the shortest diffracted path, and $R_\mathrm{d}$ the length of the direct
path, from source to receiver (Maher et al., 2021).

## Quickstart

1. Obtain
   ```
   git clone https://github.com/liamtoney/infresnel.git
   cd infresnel
   ```

2. Create environment, install, and activate
   ```
   conda env create
   conda activate infresnel
   ```

3. Run using the Python interpreter
   ```python
   python
   >>> from infresnel import calculate_paths
   ```

## Installation details

*infresnel*'s dependencies are [PyGMT](https://www.pygmt.org/latest/) (for simplified
SRTM data downloading and caching) and
[rioxarray](https://corteva.github.io/rioxarray/stable/) (for DEM file I/O,
reprojection, and elevation profile interpolation). These dependencies are listed in
[`environment.yml`](environment.yml), and they are installed in step 2 above. Note that
[Matplotlib](https://matplotlib.org/) is an **optional** dependency; it's used only for
plotting.

You might want to install *infresnel* into an existing
[conda](https://docs.conda.io/en/latest/) environment, instead of making a new one. In
this case, after step 1 above run
```
conda activate <existing_environment>
pip install --editable .
```
which uses `pip` to install *infresnel*'s dependencies, if you don't already have them
installed in your existing environment.

In either installation case, your installation will be "editable" — meaning, you can
modify the source code in your local `infresnel/` directory — or run a `git pull` to
update with any new remote changes — and the installed package will be updated.

## References

Maher, S. P., Matoza, R. S., de Groot-Hedlin, C., Kim, K., & Gee, K. L. (2021).
Evaluating the applicability of a screen diffraction approximation to local volcano
infrasound. *Volcanica*, *4*(1), 67–85. https://doi.org/10.30909/vol.04.01.6785
