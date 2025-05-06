# KdVSolitonGas.jl

This repository contains the code for the paper "Effective computation of soliton gas primitive potentials": [https://arxiv.org/abs/2505.02029](https://arxiv.org/abs/2505.02029).

**Installation**

To install as a Julia package, run the following lines:
```
using Pkg
Pkg.add(url="https://github.com/cade-b/KdVSolitonGas.jl")
```
The package can then be loaded via the command
```
using KdVSolitonGas
```
**Basic usage**

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
**Soliton addition**

Additional solitons can be inserted into the soliton gas by specifying additional arguments when calling ```precompute```. For instance a single soliton with spectral parameter $$\kappa=5$$ and norming constant $$\chi=10^{-6}$$ can be inserted by running
```
κ, χ = 5., 1e-6
u = precompute(intervals,κ,χ,h,typevec)
```
Inserting multiple solitons requires inputting vectors of spectral parameters and norming constants. For instance, another soliton with spectral parameter $$\kappa=1.2$$ and norming constant $$\chi=100$$ may be inserted by running
```
κvec = [1.2; 5.]
χvec = [100; 1e-6]
u = precompute(intervals,κvec,χvec,h,typevec)
```
**Optional arguments**

The function ```precompute``` has two optional arguments that the user may modify:

```nmat``` designates the number of collocation points used on each circle and interval to solve the Riemann--Hilbert problem. By default, 120 points on each circle and 20 points on each interval are used. For solitons that have accumulated on $$n$$ pairs of intervals, ```nmat``` should be a $$n\times 2$$ matrix of integers (same dimensions as ```intervals```). ```nmat[j,1]``` designates the number of collocation points on the circle surrounding $$(a_j,b_j)$$, and ```nmat[j,2]``` designates the number of collocation points on the interval $$(a_j,b_j)$$.

```circ_size``` designates the size of the circular contours in the Riemann--Hilbert problem. The circle around $$(a_j,b_j)$$ will have diameter ```circ_size*(bⱼ-aⱼ)```. By default, ```circ_size=1.25```. One should ensure that this value is small enough so that no circles intersect each other or poles corresponding to inserted solitons.

When additional solitons are inserted, the soliton gas solution ```u(x,t)``` has several optional arguments:

```flip_tol``` denotes the threshold at which the triangularity of residue conditions is flipped. This corresponds to the constant $$c$$ in Section 3.1 and defaults to 10.

```pole_circ``` denotes the radius of the circle around each residue condition in which the solution to the pure gas problem is expanded in a Laurent series. It defaults to 0.001 and the user should ensure that it is small enough to avoid intersecting with other poles or the lenses of the pure gas problem.

```flip``` or ```flips``` allows the user to manually designate which residue conditions should be flipped. For a single soliton addition, ```flip=1``` will flip the residue conditions and ```flip=0``` will preseve them. For multiple solitons, flips should be a binary vector with entries corresponding to each soliton.

```max_deriv_terms``` designates the maximum number of terms used in the Laurent expansion of the solution to the pure gas problem around each residue. It defaults to 25. 

```verbose``` is a flag that allows the user to diable warnings. This is done by setting ```verbose=false```.
