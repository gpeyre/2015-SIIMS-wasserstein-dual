This toolbox reproduces the numerical results of the paper:

A Smoothed Dual Approach for Variational Wasserstein Problems,
Marco Cuturi and Gabriel Peyré
Preprint arXiv:1503.02533
2015

These code reproduces the experiments on the TV regularization:
* test_shapes.m: barycenters of images with regularization.
* test_gradflow.m: same but for gradient flows instead of barycenters.
* test_meg.m: barycenters for MEG data, i.e. defined on a graph.

These codes are 
* compute_tv_barycenters.m: main function to compute regularized barycenter using FB or BFGS.
* compute_barycenter_grad.m: compute the gradient of the dual cost function for barycenters.
* compute_dual_wasserstein.m: compute the gradient of the dual Wasserstein distance. 