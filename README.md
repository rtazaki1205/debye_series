# Introduction

This package provides a light scattering solution of a homogeneous sphere by means of the Debye series expansion 
(see also Appendix A in [Tazaki et al. 2021](https://ui.adsabs.harvard.edu/abs/2021arXiv210107635T)).
The Debye series is a geometric expansion of Lorenz-Mie coefficients with respect to the reflection coefficients in a rigorous manner. 
Thus, each terms of the series represents a light scattering component, such as diffraction, surface reflection, transmitted light, and internally reflected lights.

If all terms of the Debye series are used, the solution exactly recovers the one obtained by the Lorenz-Mie theory.
Also, the short-wavelength limit of the Debye series is equivalent to the geometrical optics approximation for a sphere, and therefore, 
it can be considered as a generalized version of the geometrical optics approximation.

The codes adopt the algorithm developed by [Shen & Wang (2010)](https://ui.adsabs.harvard.edu/abs/2010ApOpt..49.2422S).

# Terms of use

The codes are distributed under the [MITlicense](https://opensource.org/licenses/MIT) and can be used, changed
and redistributed freely. If you use this package to publish papers, please cite the relevant papers.

# History

 - Initial release
