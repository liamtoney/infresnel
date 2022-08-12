# *infresnel*

*infresnel* is a Python package which facilitates computation of the
[Fresnel number](https://en.wikipedia.org/wiki/Fresnel_number) for
[infrasound](https://en.wikipedia.org/wiki/Infrasound) applications — specifically,
propagation over topography at local to regional scales.

## Background and motivation

The dimensionless Fresnel number $N$ for an acoustic wave with wavelength $\lambda$ is
given by $$N = \frac{(R_\mathrm{f} - R_\mathrm{d})}{\lambda / 2},$$ where $R_\mathrm{f}$
is the length of the shortest diffracted path over topography, and $R_\mathrm{d}$ the
length of the direct path (i.e., line-of-sight slant distance), from source to receiver
(Maher et al., 2021). In words, the Fresnel number is the extra distance a wave must
travel due to topography, normalized by half a wavelength. The path length difference
$(R_\mathrm{f} - R_\mathrm{d})$ can additionaly be used to estimate travel time delays
due to diffraction over topography. The travel time delay $\delta_\mathrm{f}$ is given
by $$\delta_\mathrm{f} = \frac{(R_\mathrm{f} - R_\mathrm{d})}{c},$$ where $c$ is the
estimated acoustic wavespeed assuming a homogenous atmosphere.

These are simple equations, but the practical computation of the quantity
$(R_\mathrm{f} - R_\mathrm{d})$ is somewhat involved. The goal of *infresnel* is to make
this computation as quick and convenient as possible.

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
which uses [pip](https://pip.pypa.io/en/stable/) to install *infresnel*'s dependencies,
if you don't already have them installed in your existing environment.

In either case, your installation will be "editable." This means that you can modify the
source code in your local `infresnel/` directory — or run a `git pull` to update with
any new remote changes — and the installed package will be updated.

## References

Maher, S. P., Matoza, R. S., de Groot-Hedlin, C., Kim, K., & Gee, K. L. (2021).
Evaluating the applicability of a screen diffraction approximation to local volcano
infrasound. *Volcanica*, *4*(1), 67–85. https://doi.org/10.30909/vol.04.01.6785
