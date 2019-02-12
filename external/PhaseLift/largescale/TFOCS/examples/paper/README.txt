These scripts generate the figures in the paper (with perhaps less formatting...)
TFOCS version 1.0a, code by Michael Grant (mcg@cvxr.com) and Stephen Becker (srbecker@caltech.edu)

Below is a list of figures and the corresponding code used to generated them
Before running any of the test codes, please run "setup_path_test.m" in this
directory, and it will add a necessary folder to the path.

Note: so that figure #'s are consistent, the figure numbers I refer to here are
from the "TFOCS_Oct2010_reference.pdf" copy included in this file (since, e.g.,
a journal version of the paper may be a bit different).

Figure 1: Trivial to generate, so not including code for this.

Figure 2: Exact penalty property for Dantzig Selector. see
"ExactPenaltyTest_Fig2" folder

Figure 3: Accelerated continuation, outer iterations only, on LASSO. see
"ContinuationOuter_Fig3" folder

Figure 4: Accelerated continuation, inner iterations, on Dantzig.  see
"ContinuationDantzig_Fig4" folder

Figure 5: Test on a strongly convex quadratic problem.  see
"StronglyConvexQuadratic_Fig5"

Figure 6: Comparison of gradient descent variants on the Dantzig Selector.  see
"CompareFirstOrderSolvers_Fig6" folder

Figure 7: Comparison of solver vs. SPGL1 on LASSO.  see "CompareWithSPGL1_Fig7" folder

Figure 8: Wavelet + TV denoising.  see "TV_and_Analysis_denoising_Fig8" folder

Figure 9: Noiseless matrix completion.  see "MatrixCompletion_Figs9-10" folder, and "noiselessMatrixCompletion.m"

Figure 10: Noisy matrix completion.  see "MatrixCompletion_Figs9-10" folder, and "nosisyMatrixCompletion.m"

Figure 11: Radar signal via RMPI.  see "RMPI_Fig11" folder
