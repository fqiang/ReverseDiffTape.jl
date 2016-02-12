# ReverseDiffTape 

[![Build Status](https://travis-ci.org/fqiang/ReverseDiffTape.jl.svg?branch=master)](https://travis-ci.org/fqiang/ReverseDiffTape.jl)

Welcome to the ReverseDiffTape.jl, a Julia package for reverse mode differentiation on a tape.

#Overview
This julia package implements reverse mode automatic/algorithmic differentiation algorithm for computing gradient and Hessian of a scalar function. The Hessian evaluation and pattern finding use a state-of-the-art reverse Hessian algorithm, named <strong>Edge_Pushing<strong>.


#Interface
- Building tape object from the operator overloading interface. 
- Building tape object from the Julia expression object.

#Computational Graph
Internally, the ReverseDiffTape represents each computational graph using three arrays. 
- An index array is use to record using Postfix notation.
- A value array to record the parameter values.
- A value array to supply the independent variable values.


Jump to see [Examples](https://github.com/fqiang/ReverseDiffTape.jl#example-of-using-this-package) in the tests script of using this package. 

#Operator overloading interface

* AD(data)
    
    creating an AD type object from data array. The data array is already in postfix notation. 

* AD_P( p, val )
    
    creating an AD type object that represents a fixed parameter with value equals to val. val is appended to p and its index on p is recorded in this AD object.  

* AD_V( x, val )

    creating an AD type object that repsents an independent variable. val is appended to x array and its index on x is recorded in this AD object. 

AD type is a wrapper of a data array which repsents the computation graph in postfix notation. This allow the postfix array to be build easily. For example,

```julia
using ReverseDiffTape
p = Vector{Float64}()
x = Vector{Float64}()
x1 = AD_V(x, 1.1)
x2 = AD_V(x, 2.2)
x3 = AD_V(x, 3.3)
p1 = AD_P(p, 1.0)
p2 = AD_P(p, 2.0)
c = sin(x1)+cos(x2^p2) * p1 - x3*p2
@show c.data
```
Then, c.data will give us the correponding postfix notation for the function expression. 

* tapeBuilder(data)
    
    returns a tape object from the postfix data array.

#Julia expression interface

* tapeBuilder(expr, p)
    
    returns a tape object from Julia expression object. The fixed parameter values are pushed into the parameter value array p. It is assume that independent variable in the Julia expression is represented by the ref symbol.


#Tape operations

Once we have a tape object, the function values and gradient and Hessian matrix can be evaluated by the following interface operations. 

* Function evaluation
    - feval(tape, x, p)

* Gradient evaluation
    - grad_reverse(tape, x, p)

        evaluating the gradient vector. The corresponding nonzero values are recorded in tape.g. 

    - grad_structure(tape)

        evaluating the gradient structure. Nonzero indicies are recorded in tape.g_I. 

* Hessian evaluation
    - hess_structure2(tape)

        evaluating the sparsity pattern of the Hessian matrix. The nonzero row and column indicies are recorded in tape.h_I and tape.h_J correspondingly.

    - hess_reverse2(tape,x,p,factor=1.0)

        evaluating the nonzero values of the Hessian matrix. The nonzero values are recorded in tape.hess. 


#Example of using this package

Just type 
```julia 
using ReverseDiffTape 
``` 
in julia console, then you are ready to try the examples below.

##Example
- To evaluate function `sin(x1)+cos(x2^2)*1-x3 `, given `x1=1.1, x2=2.2, x3=3.3`

```julia
p = Vector{Float64}()  #a empty vector of parameters
x = Vector{Float64}()  #a empty vector of independent variables
x1 = AD_V(x, 1.1)      #creating x1
x2 = AD_V(x, 2.2)      #creating x2
x3 = AD_V(x, 3.3)      #creating x3
p1 = AD_P(p, 1.0)      #creating a parameter 1.0
p2 = AD_P(p, 2.0)      #creating a parameter 2.0
c = sin(x1)+cos(x2^p2) * p1 - x3*p2   #make a function expression sin(x1)+cos(x2^2)*1.0 - x3*2.0
tt = tapeBuilder(c.data)              #building the tape object 

val = feval(tt,x,p)    #compute the function value
grad_structure(tt)     #compute the nonzero indicies in gradient vector
grad_reverse(tt,x,p)   #compute the nonzero values in gradient vector
hess_structure2(tt)    #compute the nonzero indicies in Hessian matrix
hess_reverse2(tt,x,p)  #compute the nonzero values in the Hessian matrix

@show val
@show tt.g_I, tt.g
@show tt.h_I, tt.h_J, tt.hess
```

#References: 
R.M. Gower and M.P. Mello. "A new framework for the computation of Hessians", Optimization Methods and Software 27-2, pp. 251–273, 2012. [paper](http://www.ime.unicamp.br/rel_pesq/2010/rp16-10.html)
