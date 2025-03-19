Inverse Medium Problem Solver

This solver uses the forward scattering volume solver of the paper:

Felipe Vico, Leslie Greengard, Miguel Ferrando, Fast convolution with free-space Green's functions, Journal of Computational Physics, Volume 323, 2016, Pages 191-203, ISSN 0021-9991, https://doi.org/10.1016/j.jcp.2016.07.028.

For the inverse solver, we recast the inverse problem as an optimization problem with multiple frequencies. Then we solve from the lowest to the highest frequency each single frequency inverse problem using steepest descent. 
The solution at frequency k is the initial guess at frequency k+dk. For details in continuation in frequency, see

Inverse scattering via Heisenberg's uncertainty principle, Yu Chen, Inverse Problems, Volume 13, Number 2

High resolution inverse scattering in two dimensions using recursive linearization, C Borges, A Gillman, L Greengard - SIAM Journal on Imaging Sciences, 2017.

The domain is represented by the product of sine functions. For details on the domain representation see:

High resolution inverse scattering in two dimensions using recursive linearization, C Borges, A Gillman, L Greengard - SIAM Journal on Imaging Sciences, 2017.
