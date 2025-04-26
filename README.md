# KdVSolitonGas.jl

This repository contains the code for the paper "Effective computation of soliton gas primitive potentials": [arXiv link here](arXiv link here).

To install as a Julia package, run the following lines:
```
using Pkg
Pkg.add(url="https://github.com/cade-b/KdVSolitonGas.jl")
```
The package can then be loaded by running
```
using KdVSolitonGas
```

Soliton gas solutions are generated with the function ```precompute```. This function requires at least 3 arguments:

```intervals``` designates the intervals $$\Sigma_+$$ on which the solitons have accumulated in the spectral plane. For instance, $$\Sigma_+=(0.5,1)\cup(1.5, 2.5)\cup(3,4.5)$$ is obtained by setting
```
intervals = [0.5 1.; 1.5 2.5; 3 4.5]
```
The manner of accumulation is designated by a function ```h``` and a vector ```typevec```, which together define the function 
```math
r(z)=2r_1(\mathrm{i} z) = h_j(z)(z-a_j)^{\alpha_j}(b_j-z)^{\beta_j},\quad z\in(a_j,b_j).
```
```h``` should take a positive integer as input and return a function of $$z$$. For instance, $$h_1(z)=1$$, $$h_2=\exp(z-2)$$, $$h_3=1+(x-4)^2$$ is obtained by defining
```
function h(j)
    if j == 1
        return z -> 1
    elseif j == 2
        return z -> exp(z-2)
    elseif j == 3
        return z -> 1+(z-4)^2
    end
end
```
The $$j$$ th entry of ```typevec``` designates the endpoint behavior $$\alpha_j,\beta_j$$. 1 corresponds to $$\alpha_j=\beta_j=-\frac{1}{2}$$ (Chebyshev 1st kind), 2 corresponds to $$\alpha_j=\beta_j=\frac{1}{2}$$ (Chebyshev 2nd kind), 3 corresponds to $$\alpha_j=\frac{1}{2}$$, $$\beta_j=-\frac{1}{2}$$ (Chebyshev 3rd kind), and 4 corresponds to $$\alpha_j=-\frac{1}{2}$$, $$\beta_j=\frac{1}{2}$$ (Chebyshev 4th kind). For example, one may define
```
typevec = [1; 2; 3]
```
With these defined, a soliton gas solution is obtained by running
```
u = precompute(intervals,h,typevec)
```
The output ```u``` denotes the solution $$u(x,t)$$ and may be evaluated at any point in the $$(x,t)-plane$$. For instance, ```u(-10,0.01)``` will return $$u(-10,0.01)$$. 

A solution plot at a given time can then be generated via the following code:
```
using Plots, LaTeXStrings
xv = -10:0.05:10
uv = u.(xv,0.) #solution at time zero
plot(xv,real.(uv), fillrange = -1, fillalpha = 0.35, linewidth = 2, label = false)
xlabel!(L"x")
ylabel!(L"u(x,0)")
```
