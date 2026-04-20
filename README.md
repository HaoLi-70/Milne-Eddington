# Milne-Eddington inversion code

---

## Overview

- This project is a pixel-by-pixel Stokes inversion code written in C and paralleled with MPI, follows the analytic solutions described in Introduction to Spectropolarimetry by J. C. del Toro Iniesta (2003), and Polarization in Spectral Lines by E. Landi Degl’Innocenti & M. Landolfi (2004). The Levenberg–Marquardt algorithm is employed to solve the nonlinear least-squares fitting problem. 
 


- This code has been applied to spectropolarimetric observations of the Fe I 15648.5 $\rm \AA$ obtained by Goode Solar Telescope/NIRIS. A brief description of the inversion code can be found in: 

  - *Qiao, Fangfang; Li, Hao; Wang, Jiasheng; Duan, Yadan; Sun, Zheng; Li, Leping (2026), Fine Structure and Formation Mechanism of a Sunspot Bipolar Light Bridge in NOAA AR 13663, The Astrophysical Journal, 1000, 2, 205* 
  
  We kindly ask you to cite this paper if this code is used.
  
---

### Geometry

- The direction of positive Stokes U corresponds to a counterclockwise rotation of 45$^\circ$ with respect to that of positive Stokes Q. Stokes V is defined as the difference between the right- and left-handed circular polarization. The observer is suppose to face the radiation source.
- Therefore, the sign and angle may need to be converted depending on the polarization defintion adopted in the observations.

---

### Parameters

- The inverted inclination angle ($\theta_B$) is defined as the angle between the magnetic field vector and the direction pointing to the observer, and the azimuth angle ($\phi_B$) is the counterclockwise from direction of positive Stokes Q. 
- Instead of using the source function ($S_0$) and its gradient ($S_1$), the code adopts: 
  
  - $S$ = $S_0$ + $S_1$
  -  $\beta$ = $S_0$/($S_0$ + $S_1$)
  
  This parameterization is chosen because S corresponds to the continuum intensity, which is less correlated with other inversion parameters.
  
- In addition, a folded boundary condition is applied to $\phi_B$ to prevent solutions from becoming trapped at the bounds.

---

## Requirements

Required:
- C compiler
- [CFITSIO](https://heasarc.gsfc.nasa.gov/docs/software/fitsio/) (for FITS file input/output)
- Make / GNU Make

Optional:
- MPI library — for parallel execution.  
- OpenBLAS — for SVD.
  (If not available, the SVD implementation from Numerical Recipes will be used.)
- GSL Scientific Library — for random number generation.
  (If not available, the Numerical Recipes RNG will be used instead; it is slower but robust.)


---

## Installation

- Clone the repository and build the program:
  ```
  git clone https://github.com/HaoLi-70/Milne-Eddington.git
  cd Milne-Eddington
  ./configure --hard
  make
  ```

- if MPI is not available:
  ```
  ./configure --hard --disable-mpi --ccompiler=gcc
  make
  ```


- if cfitsio is not found automatically, specify its paths manully in the configure file:
  ```
  lflags="-L/opt/homebrew/opt/cfitsio/lib -lcfitsio"
  iflags="-I/opt/homebrew/opt/cfitsio/include"
  ``` 


---

## Examples

- Single FITS File Inversion
  Usage: 
  ```
  mpirun -np <core_number> ../MEINV <control_file/optional>
  ```

- Example:
  ```
  cd niris_example
  mpirun -np 8 ../MEINV
  ```

- Batch FITS File Inversion: 
  Usage:
  ```
  Usage: python RunMEINV.py <data_path> np=<core_number>
  ```
- Example:
  ```
  cd Tools
  python RunMEINV.py ../niris_example/data   np=8
  ```

- If MPI is not available:
  ```
  cd niris_example
  ../MEINV
  ```

---

## Input/Output

* Input Files: 
  * niris_example.fits
    FITS file with dimensions [nl,4,nx,ny], storing Stokes profiles.
  * niriswavelength.fits: 
    FITS file with dimensions [nl], storing wavelength values.
  * input: 
    control file, storing the inversion configuration.
   
* Output Files: 
  * res.fits: 
    FITS file with dimensions [10,nx,ny], storing the inverted model parameters ($B$, $\theta_B$, $\phi_B$, $v_{\rm los}$, $\rm \Delta\lambda_D$, $a$, $\eta$, $S$, $\beta$) and the $\chi^2$.

---

## Keywords in control file

- num_run = [number]
  Number of inversions with random initializations. 
  Initial parameters are randomly generated according to the initial guess. 

- chisq_criteria = [float]
  Stopping criterion for repeated random runs. 
  Note: 
  Small value means stricter fitting quality requirements but increase computation time.

- lines =  [wavelength_of_line_center],  [effective_Lande_factor]
  Line center wavelength and effective Lande factor.
  NOTE: 
  Only normal Zeeman effect is taken into account.

- data_path = [Stokes_profile_file]
  Path to the FITS file [nl,4,nx,ny], storing the observed Stokes profiles.

- wavelength_path = [wavelength_file]
  Path to the FITS file [nl], storing the wavelength value in $\rm \AA$.

- weights = [weight_on_I], [weight_on_Q], [weight_on_U], [weight_on_V]
  Weights assigned to Stokes I, Q, U, and V.

- Bounds_[X] = [min], [max] 
  limits on the parameters [X = Bmod, Btheta, Bphi, Vlos, Dopp, Damp, Eta, Beta].

- LCoeffi = [float]
  Coefficient used to estimate the initial transverse magnetic field strength from linear polarization.

- VCoeffi = [float]
  Coefficient used to estimate the initial longitudinal magnetic field strength from circular polarization.


---

## License

- This project is licensed under the MIT License.

---
