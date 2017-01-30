Repository containing the benchmark tests for "On efficient Hessian computation using the edge pushing algorithm in Julia" by C.G. Petra, F. Qiang, M. Lubin and J. Huchette.

# Pre-requsite to run the test
Install JuMP and Ipopt
```julia
Pkg.add("JuMP");
Pkg.add("Ipopt");
```
# Running individual model in command line
Set *opt* argument to 1 to enable JuMP's coloring AD implementation, and to 2 to enable *edge_pushing* AD implementation.
* Arrow-head model

```bash
julia arrowhead.jl [n] [k] [opt]     
```
* Random sparsity model

```bash
julia -e 'include("data_gen.jl"); genRandData([n],[k]);'  #generate random data
julia random_sparsity [n] [k] [opt]
```
* Logistic regression model
```bash
julia -e 'include("data_gen.jl"); genLogData([m],[n]);'   #generate logistic regression data
julia logmod.jl [m] [n] [opt]
```

* Optimal power flow model

These acopf models are from [JuMPSupplement](https://github.com/mlubin/JuMPSupplement/tree/master/acpower)
```bash
julia opf_me.jl [nbus] [opt]  #avaliable nbus are, 662, 6620 and 66200
```

# Runing benchmark test sets from script

* Arrow-head models
```bash
run_arrowhead.jl.sh [opt] | tee arrowhead.jl.out
```
* Random sparsity models
```bash
julia -e 'include("data_gen.jl"); n=4000; for i=1:6; genRandData(n,2^i); end'           #generate random data for fixed n=4000
julia -e 'include("data_gen.jl"); k=32; n=1000; for i=1:5; genRandData(n,k);n=n*2; end' #generate random data for fixed k=32
run_random_sparsity.jl.sh [opt] | tee random_sparsity.jl.out
```
* Logistic regression models
```bash
julia -e 'include("data_gen.jl"); m=2; n=2000; for i=1:4; genLogData(m,n); n=n+2000; end'  #generate logistic model data
run_logmod.jl.sh [opt] | tee logmodel.jl.out
```
* Optimal power flow models
```bash
run_opf_me.jl.sh [opt] | tee opf_me.jl.out
```
Once the benchmark test scripts finished running, the timing stats are recorded in the *.out files. 
The Hessian evaluation time is reported in "Lagrangian Hessian" and total AD time is reported in "Total CPU secs in NLP function evaluations". Some sample output is:

    Total CPU secs in NLP function evaluations           =      0.689
    Lagrangian Hessian.................:      0.209 (sys:      0.003 wall:      0.213)

# Q&A 
* **How to obtain the statistics for the AD's timing?**
The timing statistics will be printed on the screen once the JuMP model is solved by Ipopt. Please make sure the ipopt.opt has ```print_level 5```. 
* **How to get a single evaluation timing results?**
This can be obtained from the timing statistics reported by Ipopt as well by setting ```max_iter 1``` in ipopt.opt. 
