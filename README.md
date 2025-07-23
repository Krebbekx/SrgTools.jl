# SrgTools.jl

Julia toolbox for computing and plotting the Scaled Relative Graph for dynamical systems.

In the `paper/` directory you will find scripts that reproduce all results in the journal paper:
> J.P.J. Krebbekx, R. Tóth, A. Das. "Graphical Analysis of Nonlinear Multivariable Feedback Systems," *arXiv:2507:16513*, 2025, Submitted to IEEE-TAC.

Please cite us when you generate results using our toolbox:
```
@article{krebbekxGraphicalAnalysisNonlinear2025,
  title={Graphical Analysis of Nonlinear Multivariable Feedback Systems},
  author={Krebbekx, Julius P. J. and T{\'o}th, Roland and Das, Amritam},
  note={Submitted to IEEE-TAC},
  year={2025},
  journal={arXiv:2507:16513}
}
```
The code is released under the BSD 3-Clause license.

Please feel free to ask questions, propose changes or new features by emailing to [j.p.j.krebbekx@tue.nl](mailto:j.p.j.krebbekx@tue.nl).

## Installation

To use this package, `cd` to the folder of this project, and run
1. `julia` to open the Julia REPL (read-eval-print-loop),
2. `]activate .` to activate the environment,
3. `]instantiate` to resolve the package dependencies of your environment

Now, you can include `using SrgTools` in any script as desired. 

In VSCode, select `Julia env: SrgTools.jl` in the bottom left corner. A VSCode workspace will remember this choice of Julia environment for the folder.

To install the latest version of Julia in a safe and easy way (Unix), run
```
curl -fsSL https://install.julialang.org | sh
```
For Windows, one achieves the same by running
```
winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore
```

For more information on installing Julia, see [the Julia website](https://julialang.org/install/).

## Getting Started

To get started with this toolbox, we suggest looking at the examples in `examples/` first. 

## Some Remarks on Performance

In practically all SRG calculations there is a "resolution" degree of freedom. 

In the case of `compute_lti_srg_boundary_upperhalfplane`, which underpins all LTI SRG computations, the computational burden depends on how "fine" the arrays `alphas`, `phis`, and `frequencies` are taken.

The same goes for computing interconnections of SRGs. 

## Notes

This package is still under continuous development. Therefore, things like tests and stuff are not implemented (yet). 
