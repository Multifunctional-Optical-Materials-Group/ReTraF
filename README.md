# ReTraF
Reflectance and Transmittance fitter for arbitrary layered optical systems


## **Usage**:
ReTraF fuction call takes the following form:
```
[models_out,N,D,Data_exp,Data_theor,xbest,foptions_out] = ReTraF(wl,theta,models,data_file,foptions)
```
## Input arguments:

- ```data_file``` full path to .mat file containing Reflectance and Transmitance data to fit. Variables inside this file must be:
  - ```wl_exp``` wavelength in nm.
  - ```RSample_S``` reflectance (in %) for S polarization.
  - ```RSample_P``` reflectance (in %) for P polarization.
  - ```TSample_S``` transmittance (in %) for S polarization.
  - ```TSample_S``` transmittance (in %) for P polarization.

Reflectance and transmittance variables are matrices of dimension ```(m x n)``` where ```m``` is the number of measured wavelengths and ```n``` is the number of measurements for different incident angles.
- ```wl``` wavelength (in nm) vector
- ```theta``` incident angles struct:
   - ```theta.values``` vector containing the incident angles (in degrees) of the measurements to fit.
   - ```theta.index``` vector containing the indices of the corresponding angles.
 - ```models``` cell of structs containing the information of each layer.
 - ```foptions``` options struct.

### Models struct:

