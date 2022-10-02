# ParallelProcessingHEP

All the packages required in `requirements.txt` (excepted from ROOT).

For `RDataFrame` use `runtime_measurement_rdf.py`. For `uproot` use `runtime_measurement_uproot.py`. 

Basic usage `python3 runtime_measurement_rdf.py` or `python3 runtime_measurement_uproot.py`:

- `--step` step size for the scan;
- `--max` the max value to scan;
- `--loops` how many times to run the measurements. The time mean will be taken;
- `--variable` to choose if to run for different file sizes (with `EnableImplicitMT()`and no threads selected) or a different number of threads (taking all the files);
- `--name` path tho the files;
- `--output` output folder for the results (a csv file with the measurements);
- `--n_files` if used with `threads` mode fix the number of files to be used;
- `--n_threads` if used with `size_mt` mode will allow to fix the number of threads.

A first warming up measurement is run but ignored. Then for `loops` times the load and filtering are repeated and timed. This will produce a `csv` file in the `output` folder with all the different measurments with the corresponding `variable` value.
