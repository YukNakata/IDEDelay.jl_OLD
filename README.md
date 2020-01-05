# IDEDelay

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://YukNakata.github.io/IDEDelay.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://YukNakata.github.io/IDEDelay.jl/dev)
[![Build Status](https://travis-ci.com/YukNakata/IDEDelay.jl.svg?branch=master)](https://travis-ci.com/YukNakata/IDEDelay.jl)

"IDEDelay" is a package to numerically compute solutions of Delay Differential Equations, in particular, including "Distributed Delay".
 
# Requirement
 
* Julia (VERSION ≥ v"1.0")
 
# Installation
 
Install IDEDelay 
 
 Entering Pkg mode ()
```bash
add https://github.com/YukNakata/IDEDelay.jl
```
 
# Example

## scalar Distributed DDE

### 1 
```bash
using IDEDelay
using Plots

tspan = [0 20];
idefun(t,y,int) = -2.5*int;
K(t,s,y)        = sin(y); 
delays(t)       = t-1;
history(t)      = 1.5;
stepsize		= 1e-3;

sol = ide_delay_rk(idefun,K,delays,history,tspan,stepsize);

# Output vector of times t
# Output vector of solutions y
#print(sol[1])
#print(sol[2])

# create a line plot
plot(sol, linewidth=2, xlabel="t", ylabel="y", title="IDE Delay Runge-Kutta")
sol
```

### 2

```bash
using IDEDelay
using Plots

tspan = [0 25];
idefun(t,y,int) = -15*int;
K(t,s,y)        = sin(y); 
delays(t)       = t-1;
history(t)      = 1.5;
stepsize		= 1e-3;

sol = ide_delay_rk(idefun,K,delays,history,tspan,stepsize);

# Output vector of times t
# Output vector of solutions y
#print(sol[1])
#print(sol[2])

# create a line plot
plot(sol, linewidth=2, xlabel="t", ylabel="y", title="IDE Delay Runge-Kutta")
```

## Distributed DDE system

### 1 

```bash
using IDEDelay
using Plots

tspan = [0 25];
idefun(t,y,int) = [ int[1];
                    int[2];
                    int[3] ]
K(t,s,y)        = [ -2.5*sin(y[1])-sin(y[2]);
                    -5.5*sin(y[2])-sin(y[1]);
                    -15*sin(y[3]) ];
delays(t)       = t-1;
history(t)      = [ 1.5;
                    -1.0;
                    1.5 ];
stepsize		= 1e-3;

sol = ide_delay_rk(idefun,K,delays,history,tspan,stepsize);

# Output vector of times t
# Output vector of solutions y
#print(sol[1])
#print(sol[2])

# create a line plot

plot(sol, linewidth=2, xlabel="t", ylabel="y", title="IDE Delay Runge-Kutta")
```

### 2

```bash
using IDEDelay
using Plots

tspan = [0 50];
idefun(t,y,int) = [ int[1];
                    int[2]]
K(t,s,y)        = [ -6*sin(y[1])-9*sin(y[2]);
                    -9*sin(y[1])-7*sin(y[2])];
delays(t)       = t-1;
history(t)      = [5.5*sin(t);
                    -5.0*sin(t)];
stepsize		= 1e-3;

sol = ide_delay_rk(idefun,K,delays,history,tspan,stepsize);

# Output vector of times t
# Output vector of solutions y
#print(sol[1])
#print(sol[2])

# create a line plot
plot(sol, linewidth=2, xlabel="t", ylabel="y", title="IDE Delay Runge-Kutta")
``` 
# Note

* Extension to a system of RE and DDE is planned.

 
# Author
 
* Alexander Lobaskin (Saint Petersburg State University)
* Alexey Eremin (Saint Petersburg State University)
* Yukihiko Nakata (Shimane University)
 
# License
 
"IDEDelay" is under [MIT license](https://en.wikipedia.org/wiki/MIT_License).
 