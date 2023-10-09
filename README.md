# Spillover

This repository contains simulation files used in the publication Trpevski, et al. "Glutamate spillover provides robust all-or-none behavior of plateau potentials in multicompartment models of striatal projection neurons", Authorea. March 29, 2023. DOI: 10.22541/au.168012861.16093883/v1

The simulations are implemented in Python 3, using the NEURON simulator through its Python 3 interface. This README file provides instructions for running the simulation scripts and reproducing the figures in the paper.

## Requirements

- NEURON simulator
- Python 3
- Required Python packages (list them here)

## Usage

The repository contains the following files:

- `fig_2B.py` and `fig_2B_hm.py` for reproducing plots in Fig. 2B.
- `iv.py` for generating Fig. 3.
- `single_run.py` for running a single simulation with clustered inputs.
- `nmda_plateaus.py` for obtaining results for variable cluster size (as in Fig. 2).
- `alpha_mpi.py` and `eta_mpi.py` for generating full simulated data, which were run on the [Fenix Infrastructure resorces](https://www.cscs.ch/) for 5 hours with the setup in the files `alpha.sh` and `eta.sh`, respectively.
- `analyze_alpha.py` and `analyze_eta.py` for analyzing the full simulated data.

The simulation code is organized around two classes:

- `d1msn.py` for the cell (inherits from the general class `neuron_cls.py`)
- `spillover_experiment.py` for the experiment performed with the cell (inherits from the general class `experiment.py`)

All parameters in the neuron model and the simulations are set in the `parameters.py` file.

## Running a Simulation

To run an elementary simulation, such as `single_run.py`, follow these steps:

1. Create a neuron using the `d1msn` class.
2. Create an experiment with that neuron using the `spillover_experiment` class.
3. Set up the input to the neuron using the `insert_synapses()` method. Important arguments for this method are:
   - `noise_SPN` - inserts the background noise
   - `my_spillover` - creates a cluster of synapses with extrasynaptic NMDARs
   - `no_spillover` - creates a cluster of synapses without extrasynaptic NMDARs
4. Set up the recording via the `set_up_recording()` method.
5. Run the simulation via the `simulate()` method.
6. Plot the results via the `plot_results()` method.

## Results

The results of the simulations are provided in the `results` folder. This includes the simulated data used in the `fig_2B.py` and `fig_2B_hm.py` scripts.

## DOI

The DOI for the publication is DOI: [10.3389/fncel.2023.1196182]

## License

The code is released under GNU General Public License v3.0.

## References

- [[NEURON Simulator documentation link]](https://neuron.yale.edu/neuron/)
