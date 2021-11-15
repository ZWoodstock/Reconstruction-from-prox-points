===============DIRECTIONS FOR RECREATING RESULTS===================
Run ex1_script in Octave with all subfunctions and files in the
same directory.

===============DESCRIPTIONS OF EACH FILE IN THE DIRECTORY==========
ex1_data
- Octave file containing important variables. To load, run
  "load 'ex1_data'". File contains:
  	-x_true = Original solution \overline{x} used to create
	observations
	-E = matrix whose columns store e_{k,j} from the problem
	description. E was randomly generated, so it must be saved
	in order to reproduce the same results.
	-xnew_h = Approximate recovery without using affine
	extrapolation (determined by terminating after 1000
	iterations).
	-xnew_k = Approximate recovery with affine extrapolation
	(determined by terminating after 1000 iterations).
	-x_sol = Solution x_{\infty} to the problem defined in the
	numerics section (consistent with all observations derived
	from x_true). This solution is numerically approximated by
	running 'ex1_script' for 2500000 iterations, until both
	versions of the algorithm agreed.
	-err_h = norm(xnew_h - x_sol), taken every 10 iterations.
	-err_k = norm(xnew_k - x_sol), taken every 10 iterations.
	-err_iters = convenient vector for plotting errors, e.g.
	to get an error-vs-iteration plot begin by using
	semilogy(err_iters,err_h,err_iters,err_k).

ex1_script
- Main script for running both versions of the algorithm for
  this particular problem (with and without affine extrapolation).
  The following variables may be of interest.
    -block_size : dimension of \mathcal{G}_k (10 in experiment)
    -max_evals : integer used to determine the maximum iteration
    count (max_it : round(max_evals / block_size). 
    -xnew_h : current iterate without affine extrapolation
    -xnew_k : current iterate with affine extrapolation
    -lambda_list_h : extrapolation parameter for algorithm without
    affine extrapolation.
    -lambda_list_k : extrapolation parameter for algorithm with
    affine extrapolation.
    -err_h and err_k : same description as in 'ex1_data'.   
    -df: difference between the current iterates of both
    instantiations of the algorithm. Since both versions are
    solving the same problem, df should go to zero; this is
    periodically printed to the screen for diagnostic purposes.

bandlim
- Bandlimits a real signal to a prescribed Fourier frequency band.

Q_haug
- Computes the Haugazeau projection onto the intersection of
  half-spaces defined in the paper.

subgrad_proj
- Computes a subgradient projection onto a lower-level set of a
function. Requires input of a function, a selection of the
function's subgradient, and a level set parameter. 

tvl1_loss
- Computes the function evaluation and subgradient selection for the
1-D total variation (1-norm of the finite differences). To be used
with 'subgrad_proj'.

project_monotone
- Computes the projection onto the isotonic cone.
