# Milne-Eddington inversion code

---

## Overview

- A pixel-by-pixel Stokes inversion code, written in C and paralleled with MPI, follows the analytic solutions described in J. C. del Toro Iniesta (2003) and E. Landi Degl’Innocenti & M. Landolfi (2004). The Levenberg–Marquardt method is used to solve the nonlinear least-squares problem.

### Geometry

- The direction of positive Stokes U corresponds to a counterclockwise rotation of 45$^\circ$ with respect to that of positive Stokes Q, while Stokes V represents the difference between the right- and left-handed circular polarization. The observer is suppose to face the radiation source.

### Parameters

- The inverted inclination angle ($\theta_B$) is defined as the angle between the magnetic field vector and the direction pointing to the observer, and the azimuth angle ($\phi_B$) is the counterclockwise from direction of positive Stokes Q. The code adopts the sum ($S$ = $S_0$ + $S_1$) and a ratio $\beta$ = $S_0$/($S_0$ + $S_1$) instead of the source function ($S_0$) and its gradient ($S_1$), since S corresponds to the continuum intensity, which is uncorrelated with other parameters. In addition, a folded boundary condition was applied to $\phi_B$ to prevent solutions from becoming trapped at the bounds.

---

## Requirements

- C/C++ compiler
- [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/) library (for FITS file input/output)
- MPI library (for parallel execution)
- Make / GNU Make

---

## Installation

- Clone the repository and build the program:

```
git clone https://github.com/HaoLi-70/Milne-Eddington.git
cd Milne-Eddington
./configure --hard
make
```

---

## Example

```
cd niris_example
mpirun -np 4 ../MEINV
```

---
## Input/Output

* Input: 
  * niris_example.fits: FITS file [nl,4,nx,ny] containing the Stokes profiles.
  * niriswavelength.fits: FITS file [nl] containing the wavelength.
  * input: control file.
   
* Output: 
  * res.fits: FITS file [10,nx,ny] containg the inverted model parameters ($B$, $\theta_B$, $\phi_B$, $v_{\rm los}$, $\rm \Delta\lambda_D$, $a$, $\eta$, $S$, $\beta$) and the $\chi^2$.

---

## License

- This project is licensed under the MIT License.

---