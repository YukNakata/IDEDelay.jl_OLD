# IDEDelay

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://YukNakata.github.io/IDEDelay.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://YukNakata.github.io/IDEDelay.jl/dev)
[![Build Status](https://travis-ci.com/YukNakata/IDEDelay.jl.svg?branch=master)](https://travis-ci.com/YukNakata/IDEDelay.jl)

"IDEDelay" is a package to numerically compute solutions of Delay Differential Equations. 
 
# Requirement
 
* Julia (VERSION ≥ v"1.0")
 
# Installation
 
Install IDEDelay 
 
```bash
import Pkg
Pkg.add("https://github.com/YukNakata/IDEDelay.jl")
```
 
# Example
 
```bash
using IDEDelay
using Plots

tspan = [0 40];
idefun(t,y,int) = int;
K(t,s,y)        = -8.5*sin(y)*exp(-(t-s));
delays(t)       = t-1;
history1(t)     = 1.5*sin(t);
tolerance = 1e-2;

sol = ide_delay_rk(idefun,K,delays,history1,tspan,tolerance);

plot(sol)
plot(sol, linewidth=3, xlabel="t", ylabel="y", title="IDE Delay Runge-Kutta")

sol
```
 
# Note

* Extension to the system of RE and DDE.

 
# Author
 
* Alexander Lobaskin (Saint Petersburg State University)
* Alexey Eremin (Saint Petersburg State University)
* Yukihiko Nakata (Shimane University)
 
# License
 
"IDEDelay" is under [MIT license](https://en.wikipedia.org/wiki/MIT_License).
 