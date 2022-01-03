# NMRDataSetup
This package extracts from a Bruker 1D 1H NMR experiment the following:
- a Julia complex-valued time series, `s_t`
- the sweep width, in ppm, `SW
- the sampling frequency, in Hz, `fs`
- the 0 ppm frequency, in Hz, `ν_0ppm`.

Given multiple experiments, this package also additively (in a weighted manner for each experiment) combines the experiments to get a single experiment. It is up to the user to ensure that the input experiments have the same `fs` and `SW`, otherwise the output won't be a realistic NMR spectrum. This combination is referred to as weighted mixing.

See the two scripts in `./examples` for a demo of these two tasks.