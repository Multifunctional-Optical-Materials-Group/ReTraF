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
All the required information of each layer is stored inside a ```model``` struct. We make a distinction between known layers and unknown layers.
  - Known layers:
    - Forouhi-Bloomer model (max. 4 oscillators):
      - ```model.type = "Fh-1"``` or ```model.type = "Fh-1"``` or ```model.type = "Fh-2"``` or ```model.type = "Fh-3"``` or ```model.type = "Fh-4"```
      - ```model.Eg``` Bandgap in Ev
      - ```model.n0``` low frequency refractive index
      - ```model.fi``` fi parameter (length shoul be equal to the number of oscillators)
      - ```model.fi``` Ei parameter (length shoul be equal to the number of oscillators)
      - ```model.Gi``` Gi parameter (length shoul be equal to the number of oscillators)
      - ```model.D``` layer thickness in nm
    - Real Caucy model
      - ```model.type = "Ch-n"```
