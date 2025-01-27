# Flavio with FCC-ee Higgs and electroweak precision observables

This package extends flavio v2.6.1 by introducing the projected FCC-ee sensitivities for a number of Higgs and electroweak precision observables (see also the related [smelli_fcc repository](https://github.com/eetuloisa/smelli_fcc)).

Whilst we eventually hope to incorporate the features as part of the main program, currently the best way to access the FCC-ee sensitivity estimates is to 
1. Create a virtual environment with a functioning flavio v2.6.1 default installation.
2. Make a copy of the 'flavio' directory found in this repository and replace the default flavio installation with the FCC-ee version found here.
3. The package can now be used in the same way as the unmodified flavio package. 

For the full set of observables added, see the appendices in [arXiv:2501.08321](https://arxiv.org/abs/2501.08321). 
If new to flavio, please consult the main [flavio page](https://github.com/flav-io/flavio) and the references therein.

## Information on parametric theory uncertainties

In contrast to the main flavio program, the FCC-ee theory uncertainties are treated as 'extra sources of errors' on the measurements. 
This choice was made to minimise the number of changes to the core program.
As a consequence, when accessing the FCC-ee observables *directly with flavio*, the user should execute the line 
```
flavio.parameters.read_file_values('/path/to/flavio/data/fcc_projected_parameters.yml', flavio.default_parameters)
```
at the beginning of the script or notebook. 
*When using smelli* to construct a likelihood, this command is not required.
