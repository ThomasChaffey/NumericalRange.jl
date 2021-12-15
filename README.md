# NumericalRange

Provides the Julia function

        nrange(A; nk = 1, thmax = 32, noplot = false)

which evaluates and plots the numerical range of the
`nk` largest leading principal submatrices of `A`, using `thmax`
equally spaced angles in the complex plane.
The defaults are `nk = 1` and `thmax = 32`
The eigenvalues of `A` are plotted as `x`.  The function returns `f` and `e`, 
where `f` is the numerical range and `e` is a vector of eigenvalues of `A`.
Setting `noplot = true` suppresses the plot.

This function is a direct translation of `fv.m` in Professor Nick Higham's [Matrix Computation Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/2360-the-matrix-computation-toolbox),
based on an original routine by A. Ruhe.

References:
* R. A. Horn and C. R. Johnson, Topics in Matrix Analysis, Cambridge
     University Press, 1991; sec. 1.5.
* A. S. Householder, The Theory of Matrices in Numerical Analysis,
     Blaisdell, New York, 1964; sec. 3.3.
* C. R. Johnson, Numerical determination of the field of values of a
     general complex matrix, SIAM J. Numer. Anal., 15 (1978),
     pp. 595-602.


[![Build Status](https://github.com/ThomasChaffey/NumericalRange.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ThomasChaffey/NumericalRange.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ThomasChaffey/NumericalRange.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ThomasChaffey/NumericalRange.jl)
