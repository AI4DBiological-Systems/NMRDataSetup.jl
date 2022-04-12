# NMRDataSetup
This package extracts from a Bruker 1D 1H NMR experiment the following:
- a Julia complex-valued time series, `s_t`
- the sweep width, in ppm, `SW
- the sampling frequency, in Hz, `fs`
- the 0 ppm frequency, in Hz, `ν_0ppm`.
- and more, such as a dictionary of metadata provided by the NMRGlue Python library.

See the list of exported function in `\src\NMRDataSetup.jl`, and `\examples\load_experiment.jl` for an example script. The guide on this script can be found [here](https://ai4dbiological-systems.github.io/NMRDataSetup.jl/).

# Non-Julia Dependencies
Make sure the [nmrglue](https://www.nmrglue.com/) Python 3 package is installed on your system. In linux, you can try the terminal command `pip install nmrglue --user` to install this package.
