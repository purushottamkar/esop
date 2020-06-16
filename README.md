# ESOP: Epidemiologically and Socio-economically Optimal Policies via Bayesian Optimization

This repository presents an implementation of the ESOP algorithm and the VIPER epidemiological model. The accompanying paper describing the development of the ESOP algorithm as well as the VIPER (Virus-Individual-Policy-EnviRonment) epidemiological model can be accessed at [[1]][ref1]

## Setup

### Installing Python Library Dependencies

Implementations of ESOP and VIPER make use of various standard libraries. A minimal list is provided in the [requirements.txt](requirements.txt) file. If you are a pip user, please use the following command to install these libraries:

```sh
pip3 install -r requirements.txt
```

**Note about version dependency**: although the requirements file specifies version dependencies to be exact, this is to err on the side of caution. For most modules (e.g. `numpy` or `scikit-learn`), a more recent version should work too and we are not aware of any version specific dependency on the functioning of ESOP or VIPER. If you wish to avoid version-related complications, you may try using the following command instead. This does not impose a version requirement on any module and will either use the version you currently have, or else download the latest one from the pip repository in case you dont have a certain module.

```sh
pip3 install -r requirements_noversion.txt
```

### Generating population and patch files
If you wish to simply reproduce the results of the accompanying paper [[1]][ref1], skip this step. Otherwise, if you wish to repeat the experiments on freshly generated population data and demographic patches, with possibly different population sizes or parameters, please do the following
1. (Optional) Modify the file [setup.py](setup.py) to specify parameters as desired
1. Use the following command to generate population and India-specific demographic patch files
```sh
python3 setup.py
```

Doing the above should produce two pickled files, with default names as `pop_generic` and `patch_India_fit`.

## Executing ESOP with the VIPER model

The repository presents five Jupyter notebooks, each one containing a different experiment, as listed below. Please execute the various cells in the respective notebooks to execute the ESOP algorithm for these experiments.

1. [demo_optimal_initiation.ipynb](demo_optimal_initiation.ipynb): This file reproduces the experiment corresponding to Fig 2(d) in [[1]][ref1] and uses ESOP to find out the optimal initiation point of a level 5 lock-down that lasts 30 time steps, when the virus has an incubation period of 10 days. The file also reproduces the experiment corresponding to Fig 3(b) that analyzes the global convergence of ESOP.
1. [demo_constrained_lockdown.ipynb](demo_constrained_lockdown.ipynb):  This file reproduces the experiment corresponding to Fig 4(b) in [[1]][ref1] and uses ESOP to find out an optimal lockdown that starts no earlier than t = 12, lasts no longer than 40 time steps, and optimizes the sum of f_epi and f_eco (see the accompanying paper [[1]][ref1] for details).
1. [demo_multi-phase_lockdown.ipynb](demo_multi-phase_lockdown.ipynb):  This file reproduces the experiment corresponding to Fig 4(d) in [[1]][ref1] and uses ESOP to find out an optimal third lockdown given that two previous lock-downs have already taken place. The constraints here are that the third lock-down must start no earlier than 10 time steps after the second lock-down ended, and last no longer than 40 time steps, as well as optimize the sum of f_epi and f_eco (see the accompanying paper [[1]][ref1] for details).
1. [demo_quarantine_effects.ipynb](demo_quarantine_effects.ipynb): This file reproduces the experiment corresponding to Fig 5(d) in [[1]][ref1] and demonstrates how ESOP can offer better outcomes with milder lock-downs, if testing and quarantining is more aggressive (see the accompanying paper [[1]][ref1] for details).
1. [demo_India_demographics.ipynb](demo_India_demographics.ipynb): This file reproduces the experiment corresponding to Figs 6(a), 6(b) and 7(b) in [[1]][ref1] and demonstrates how demographic variations and changes in geographical distributions of people affect outcomes offered by ESOP (see the accompanying paper [[1]][ref1] for details).

## Expected Results
If executed using the population file `pop_generic` and the demographic patch file `patch_India_fit` supplied with this repository, running ESOP with VIPER using the Jupyter notebooks supplied with this repository should closely reproduce results reported in accompanying paper [[1]][ref1]. However, minor differences may arise due to varitions in the RNG engines in different machines (although ESOP does seed all its simulations).

In particular, the following results are expected to be closely reproduced. Please see the **Final Results** section in each notebook to compare these results.

1. [demo_optimal_initiation.ipynb](demo_optimal_initiation.ipynb):
	- 15533 Infected
	- 2224 Expired
	- Peak 6672
	- Recovery time 89.929371
	- Expiry time 28.792266
	- Quarantine time 18.143396
	- 14847 Quarantined
	- Objective acheived is 6672.000000
1. [demo_constrained_lockdown.ipynb](demo_constrained_lockdown.ipynb): 
	- 15630 Infected
	- 2266 Expired
	- Peak 7353
	- Recovery time 82.490871
	- Expiry time 21.882171
	- Quarantine time 11.201635
	- 14928 Quarantined
	- Objective acheived is 7353.000000 (fEpi) + 560.000000 (fEco) = 7913.000000 (overall)
1. [demo_multi-phase_lockdown.ipynb](demo_multi-phase_lockdown.ipynb): 
	- 13977 Infected
	- 1999 Expired
	- Peak 6052
	- Recovery time 82.791667
	- Expiry time 22.046523
	- Quarantine time 11.152221
	- 13303 Quarantined
	- Objective acheived is 6052.000000 (fEpi) + 350.000000 (fEco) = 6402.000000 (overall)
1. [demo_quarantine_effects.ipynb](demo_quarantine_effects.ipynb): 
	- 844 Infected
	- 119 Expired
	- Peak 706
	- Recovery time 82.900690
	- Expiry time 22.033613
	- Quarantine time 5.823460
	- 844 Quarantined
	- Objective acheived is 706.000000 (fEpi) + 192.000000 (fEco) = 898.000000 (overall)
1. [demo_India_demographics.ipynb](demo_India_demographics.ipynb): 
	- 10849 Infected
	- 0 Expired
	- Peak 1653
	- Recovery time 17.019725
	- Expiry time nan
	- Quarantine time 12.008547
	- 117 Quarantined
	- Objective acheived is 1653.000000 (fEpi) + 480.000000 (fEco) = 2133.000000 (overall)

## Contributing
This repository is released under the MIT license. If you would like to submit a bugfix or an enhancement to ESOP or VIPER, please open an issue on this GitHub repository. We welcome other suggestions and comments too (please mail the corresponding author at purushot@cse.iitk.ac.in)

## License
This repository is licensed under the MIT license - please see the [LICENSE](LICENSE) file for details.

## References
[1] Amit Chandak, Debojyoti Dey, Bhaskar Mukhoty, and Purushottam Kar. Epidemiologically and Socio-economically Optimal Policies via Bayesian Optimization. arXiv:2005.11257 [q-bio.PE], 2020 (available at [[this arXiv link]][ref1]).

[ref1]: https://arxiv.org/abs/2005.11257