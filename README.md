# ReTraF
Reflectance and Transmittance fitter for arbitrary layered optical systems


## **Usage**:
ReTraF fuction call takes the following form:
```
[models_out,N,D,Data_exp,Data_theor,xbest,foptions_out] = ReTraF(wl,theta,models,data_file,foptions)
```
where the input arguments are:

- ```data_file``` full path to .mat file containing Reflectance and Transmitance data to fit. Variables inside this file must be:
  - aa
- ```wl``` wavelength (in nm) vector
- ```theta```
