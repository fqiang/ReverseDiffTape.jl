# ReverseDiffTape 

[![Build Status](https://travis-ci.org/fqiang/ReverseDiffTape.jl.svg?branch=ad2016)](https://travis-ci.org/fqiang/ReverseDiffTape.jl)

Welcome to the ReverseDiffTape.jl, a Julia package for reverse mode differentiation on a tape.

#Overview
This julia package implements reverse mode automatic/algorithmic differentiation algorithm for computing gradient and Hessian of a scalar function. The Hessian evaluation and pattern finding use a state-of-the-art reverse Hessian algorithm, named <strong>Edge_Pushing<strong>.


#Interface with [JuMP](https://github.com/JuliaOpt/JuMP.jl/) throught [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl)
* Delegate the AD evaluation routines to use ReverseDiffTape implementation
```julia
include("interface_tape.jl")               
using TapeInterface                            
m = Model(solver=TapeSolver(IpoptSolver())) #creating TapeSolver instance
```
* Rest proceeds as the same of a normal JuMP model

```julia
@variable(m, -100<=x[1:10]<=100)
@NLobjective(m, Min, sum{ cos(sum{ x[i+j], j=1:3}) + sum{ (x[i] + x[j])^2, j=1:3} , i=1:7} )
status = solve(m)
```

# Benchmarking results

The bencharmking results are presented in The 7th International Conference on Algorithmic Differentation, titled "On efficient Hessian computation using the edge pushing algorithm in Julia". 

An improved performance results are submitted in the post-conference publication in Optimization Methods and Software. The results can be reproduced using scripts at [benchmark](https://github.com/fqiang/ReverseDiffTape.jl/tree/master/benchmark) directory of this repository. 


#References: 
R.M. Gower and M.P. Mello. "A new framework for the computation of Hessians", Optimization Methods and Software 27-2, pp. 251–273, 2012. [paper](http://www.ime.unicamp.br/rel_pesq/2010/rp16-10.html)
