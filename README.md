# MATLAB code for "Low-Complexity Blind Parameter Estimation in Wireless Systems with Noisy Sparse Signals"

If you are using the simulator (or parts of it) for a publication, then you must clearly mention it and cite our paper:

A. Gallyas-Sanhueza and C. Studer, "Low-Complexity Blind Parameter Estimation in Wireless Systems with Noisy Sparse Signals," in IEEE Transactions on Wireless Communications, doi: 10.1109/TWC.2023.3247887.

## Running simulations

- To generate figures 1, 2, 4, and 5 of the paper, simply run `figures_1_2_4_and_5.m`
  - Be patient, as the simulation takes long. (To check that the code will run properly, we advise to run fewer trials first, by setting `trials = 10` in line 11. This takes about 1 minute on a 1.8 GHz Dual-Core Intel Core i5. With 10000 trials it takes hours)

- To generate figure 3 of the paper, simply run `figure_3.m`

- Code for figures 6, 7 and 8 comming soon...

## Third-party software

- To compute the median with the quickselect algorithm, we use a modified version of the _qselect_ function from https://www.mathworks.com/matlabcentral/fileexchange/68947-qselect. We include the licence for that function inside the `qselect` folder.

- To generate our plots, we use the _shadedErrorBar_ function from https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar. We include the licence for that function inside the `shadedErrorBar-master` folder.

## Version history

Version 1: ag753@cornell.edu - initial version for GitHub release
