# CTC-sim
Open-source computational framework for mechanistic modeling and simulation of circulating tumor cell clusters during metastatic migration.

## Funding

CTC-sim is being developed with support from the [National Science Foundation CAREER Award No. 2543185](https://www.nsf.gov/awardsearch/show-award/?AWD_ID=2543185)

## Contact

Agnieszka Truszkowska  
Department of Chemical and Materials Engineering  
The University of Alabama in Huntsville  
301 Sparkman Drive, Engineering Building 117  
Huntsville, AL 35899, USA  

Email: [at0175@uah.edu](mailto:at0175@uah.edu)

## Initial configuration

After cloning the repository, run `config.py` to configure the software before compiling or running the simulations.

## Requirements and Compatibility

The software is currently developed and tested on Linux operating systems. The software is designed to be readily adaptable to macOS and Windows operating systems, and the PI will provide assistance with such adaptations to facilitate broader use of the software.

Compilation requires support for the C++17 standard. The Python scripts were developed and tested using Python 3.8.10; other Python 3.x versions may also be compatible, although they have not been explicitly tested.

The current implementation is serial and runs on a single CPU core.

## Structure

* `include` - header files
* `src` - source files
* `scripts` - utility scripts
* `simulations` - input files for all simulations
* `tests` - test suite for the code

## Documentation

After running `config.py`, source code documentation can be generated using [Doxygen](https://www.doxygen.nl/). Run Doxygen using the provided `Doxyfile`:

```text id="k5q2rs"
doxygen Doxyfile
```

The generated documentation will be placed in the `documentation` directory created by `config.py`.

Instructions for running individual simulations are provided in `README` files within the corresponding directories under `simulations`.

## Testing and validation

All software components used in the simulations were thoroughly tested. Other functionality may not have been tested to the same extent, and users are advised to use it with caution. The `tests` directory is currently partially out of sync with the codebase and will be made fully operational by the end of October 2026.

## Citation

If you find this software useful, please consider citing the following publications:

* Ezeobidi EI, Truszkowska A (2025) Modeling the dynamics of circulating tumor cell clusters inside a microfluidic channel. *Biomicrofluidics* 19(1):014103. https://doi.org/10.1063/5.0249165

## Questions and Feedback

Questions about using the software, suggestions for improvements, or reports of unexpected behavior are welcome. Feel free to open an issue or contact me directly.



