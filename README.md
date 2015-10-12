# ReverseDiffTape 

[![Build Status](https://travis-ci.org/fqiang/ReverseDiffTape.jl.svg?branch=master)](https://travis-ci.org/fqiang/ReverseDiffTape.jl)

Welcome to the ReverseDiffTape.jl, a Julia package for reverse mode differentiation on a tape.

#Overview
This julia package implements reverse mode automatic/algorithmic differentiation algorithm for computing gradient and Hessian of a scalar function. The package uses an operator overloading interface for user to build the computation graph on an Tape--an array of integers. Then the interface methods allow the function values, gradient vector, Hessian matrix and also the sparsity pattern of the gradient and Hessian to be evaluated on the indices array. 

#Highlights
- Operator overloading interface for building computational graph.
- Everything is implemented on a plain array. 
- Reverse mode algorithm for both gradient and Hessian evaluation.
- Implements a state-of-the-art reverse Hessian algorithm, named <strong>Edge_Pushing<strong>

Jump to see [Examples](https://github.com/fqiang/ReverseDiffTape.jl/blob/master/test/runtests.jl) in the tests script of using this package. 

#Package Exports
##Types
* TT_TYPE

    The array of indicies to represent the computational graph.

* TV_TYPE

    The array of real numbers for storing parameter and independent variable values.

* AD_V

    The type used to create an independent variable.

* AD_P

    The type used to create a parameter. 

* Edge

    The edge type is return for Hessian to represent a nonlinear relationships between two independent variables. 

##Interface Functions
* Function evaluation
    - feval(tt::TT_TYPE, vvals::TV_TYPE, pvals::TV_TYPE)

* Gradient evaluation
    - grad_reverse(tt::TT_TYPE,vvals::TV_TYPE,pvals::TV_TYPE)
    - nzg(tt::TT_TYPE)

* Hessian evaluation
    - reverse_hess_ep(tt::TT_TYPE,vvals::TV_TYPE,pvals::TV_TYPE)
    - nzh(tt::TT_TYPE)

#Using this Package
This package is designed to be very user friendly for using it.
Just type 
```julia 
using ReverseDiffTape 
``` 
in julia console, then you are already to try the examples below.

##Example
- To evaluate function `sin(x1)+cos(x2^2)*1-x3 `, given `x1=1.1, x2=2.2, x3=3.3`
```julia
function test1() 
    tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt, vvals, 1.1)
	x2 = AD_V(tt, vvals, 2.2)
	x3 = AD_V(tt, vvals, 3.3)
	p1 = AD_P(tt, pvals, 1)
	p2 = AD_P(tt, pvals, 2)
	c = sin(x1)+cos(x2^p2) * p1 - x3*p2
	val = feval(tt,vvals,pvals)
    return val
end
```
- Reverse gradient evaluation for function `sin(x1)+cos(x2)`, given `x1=1.1, x2=2.2`
```julia
function test2()
    tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	a = AD_V(tt,vvals,1.1)
	b = AD_V(tt,vvals,2.2)
	c=sin(a)*cos(b)
	grad = grad_reverse(tt,vvals,pvals)
    return grad
end
```
- Reverse Hessian evaluation for function `cos(x1*x2)`, given `x1=1.1, x2=2.2`
```julia
function test3()
    tt = TT_TYPE()
	pvals = TV_TYPE()
	vvals = TV_TYPE()
	x1 = AD_V(tt,vvals,1.1)
	x2 = AD_V(tt,vvals,2.2)
	c = cos(x1*x2)
	eset = reverse_hess_ep(tt,vvals,pvals)
    return eset
end
```

#Future work
- At this time, only supports unary and binary operator types, <em>i.e.</em>, +, -, *, /, ^, sin, cos. Later the packages will be extended to support more operation types. 
- Julia's multiple dispatch idea can be applied to make the code more structured and extensible. <em>i.e.</em>, by dispatching the operator symbol as types. 

#References: 
R.M. Gower and M.P. Mello. "A new framework for the computation of Hessians", Optimization Methods and Software 27-2, pp. 251–273, 2012. [paper](http://www.ime.unicamp.br/rel_pesq/2010/rp16-10.html)
