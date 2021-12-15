module NumericalRange

using LinearAlgebra
using Plots

"""
        nrange(A; nk = 1, thmax = 32, noplot = false)

Numerical range (or field of values) of the matrix `A`.
Evaluates and plots the numerical range of the
`nk` largest leading principal submatrices of `A`, using `thmax`
equally spaced angles in the complex plane.
The defaults are `nk = 1` and `thmax = 32`
The eigenvalues of `A` are plotted as `x`.  The function returns `f` and `e`, 
where `f` is the numerical range and `e` is a vector of eigenvalues of `A`.
Setting `noplot = true` suppresses the plot.

This function is a direct translation of `fv.m` in Professor Nick Higham's Matrix Computation Toolbox,
https://www.mathworks.com/matlabcentral/fileexchange/2360-the-matrix-computation-toolbox,
based on an original routine by A. Ruhe.

References:
R. A. Horn and C. R. Johnson, Topics in Matrix Analysis, Cambridge
     University Press, 1991; sec. 1.5.
A. S. Householder, The Theory of Matrices in Numerical Analysis,
     Blaisdell, New York, 1964; sec. 3.3.
C. R. Johnson, Numerical determination of the field of values of a
     general complex matrix, SIAM J. Numer. Anal., 15 (1978),
     pp. 595-602.


"""
function nrange(B; nk=1, thmax=32, noplot=false)
        function rq(A, x)
                return x'*A*x/(x'*x)
        end
        thmax -= 1 # the function uses thmax + 1 angles

        (n, p) = size(B)
        if n != p
                DimensionMismatch("Matrix must be square.")
        end

        z = Matrix{ComplexF64}(undef, 2*thmax + 1, 1)
        F = eigen(B)
        e = F.values
        f = []

        # filter out cases where B is Hermitian or skew-Hermitian
        if B == B'
                f = [minimum(e), maximum(e)]
        elseif B == -B'
                e = imag(e)
                f = [minimum(e), maximum(e)]
                e *= im
                f *= im
        else
                for m = 1:nk
                        ns = n + 1 - m
                        A = B[1:ns, 1:ns]
                        for i = 0:thmax
                                th = i/thmax*pi
                                Ath = exp(im*th)*A
                                H = 0.5*(Ath + Ath')
                                X = eigen(H, sortby = x -> real(x))
                                V = X.vectors
                                z[i + 1] = rq(A, V[:, 1])
                                z[i + 1 + thmax] = rq(A, V[:, ns])
                        end
                        f = [f; z]
                end
                # join up the boundary
                f = [f; f[1, :]]
        end
        if thmax == 0
                f = e
        end

        if !noplot
                p = plot(real(f), imag(f), framestyle=:origin, aspect_ratio = 1, legend=false)
                scatter!(p, real(e), imag(e), markershape=:x, legend=false)
                display(p)
        end

        return f, e
end

end
