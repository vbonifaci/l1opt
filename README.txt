=================================================================
MATLAB implementation of Laplacian-based gradient methods for
l1-Norm Minimization, as described in the paper: 

[B19] V. Bonifaci. A Laplacian Approach to l1-Norm Minimization. 
arXiv:1901.08836 [cs.DS] -- http://arxiv.org/abs/1901.08836

Author: Vincenzo Bonifaci, Università Roma Tre, Italy
<vincenzo.bonifaci@uniroma3.it>

These files are to be used with the *l1benchmark* MATLAB package: 
https://people.eecs.berkeley.edu/~yang/software/l1benchmark/

For more details about the l1benchmark MATLAB package, 
please refer to the paper:
[YGZ+10] A. Yang, A. Ganesh, Z. Zhou, S. Sastry, and Y. Ma. 
Fast L1-Minimization Algorithms for Robust Face Recognition. 
arXiv:1007.3753 [cs.CV] -- https://arxiv.org/abs/1007.3753
=================================================================

DESCRIPTION

The following MATLAB R2020b modules are included in this collection: 

* compare_noise_free.m
	-- a sample driver that tests the methods from [B19] against the other methods in the l1benchmark suite. This file should replace the file with the same name in the l1benchmark distribution. Please note that the revised version measures the difference in objective function value instead of the Euclidean distance in order to benchmark the algorithms, hence some of the original code had to be slightly adapted. 
* SolvePGS.m
	-- a MATLAB implementation of the Primal Gradient Scheme from [B19]
* SolveAGS.m
	-- a MATLAB implementation of the Accelerated Gradient Scheme from [B19]
* SolveAGS2.m
	-- a MATLAB implementation of the revised Accelerated Gradient Scheme from [B19]
* SolveHomotopy.m
	-- MATLAB implementation of the Homotopy method from [YGZ+10], adapted to the new measure
* SolvePALM.m
	-- MATLAB implementation of the Primal Augmented Lagrangian method from [YGZ+10], adapted to the new measure
* SolveDALM.m
	-- MATLAB implementation of the Dual Augmented Lagrangian method from [YGZ+10], adapted to the new measure
* SolvePDIPA.m
	-- MATLAB implementation of the Primal-Dual Interior Point method from [YGZ+10], adapted to the new measure
* SolveAMP.m
	-- MATLAB implementation of the Approximated Message Passing method from [YGZ+10], adapted to the new measure

