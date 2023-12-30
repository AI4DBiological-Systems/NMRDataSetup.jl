# NMRDataSetup.jl
Currently, this package is designed to work given the 1D 1H NMR time series free-induction decay (FID) data from Bruker spectrometers. The main function is `setupBruker1Dspectrum`, which requires the user to specify the folder path where the binary *fid* file and text file *acqu* (alternatively, *acqus*) are.

The following parameters must be present in the *acqu* (or *acqus*) file:
```
$TD # total number of data samples acquired.
$SW # spectral window, in ppm.
$SFO1 # carrier freuqnecy. Used as spectrometer frequency.
SO1 # offset frequency. Used to estimate roughly where the 0 ppm is.
$SW_h # sampling frequency.

$BYTORDA # control variable on endianess.
$DTYPA # control variable on data type of the FID binary file.
```

Please refer to a Bruker manual for the meaning of these parameters.

# Usage
Please see `/examples/load_experiment.jl` for a walk-through. 

# Install
Add the custom registries for dependencies, and then add the package.
```
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/AI4DBiological-Systems/PublicJuliaRegistry"))

pkg"add NMRDataSetup"
```

You also need to add this custom public registry if you want to run the example in `/examples/`:
```
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/RoyCCWang/RWPublicJuliaRegistry"))
```


# Citation
Our work is undergoing peer review. Please cite our ChemRxiv version if you use this software.
[https://doi.org/10.26434/chemrxiv-2023-0s196](https://doi.org/10.26434/chemrxiv-2023-0s196)

