## Project overview
This project includes my scripts working with statistical emulators. I was experimenting with three different packages (hetGP, laGP, mgcv) to find one to replace the current emulator (mlegp) used in PEcAn (https://github.com/PecanProject/pecan).

### Master function
The main function to work with is `fit_model()` from 'model_metrics.R'. It is a function that fits you a chosen model in given data and returns metrics of a fit in addition with different diagnostics plots. Its construction is divided into smaller scripts e.g. `helpers.R` or `autobuild_mgcv.R`.

### Other scripts in the project
In addition to `model_metrics.R` the project includes the following R.scripts:
* `helpers.R` (helper functions for `fit_model()` to use)
* `autobuild_mgcv.R` (to automatise mgcv model building)
* `plotting_functions.R` (to visualize results from experiments)
* `runs.R` (to run experiments)
* `PIT_test.R` (to draw PIT histograms)
* `dharma.R` (to do residual diagnostics using DHARMa-package)

To get information about certain script and functions in it, open the script and read the inline comments.

### How to get started?
I've made a special script called `intro.R` to introduce the basics of using these scripts mentioned above. Go it through to get an idea of the workflow I have used for experimenting between packages.
