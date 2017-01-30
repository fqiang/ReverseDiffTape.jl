option presolve 0;
param n;
param m;
set N:={1..n};     
set M:={1..m};
var theta{N} := -0.0001;
param lambda;
param y{M};
param x{M,N};

#var sum1{M} :=1.0;
#subject to con1 {i in M}:
#    y[i]*sum{j in N}theta[j]*x[i,j] == sum1[i];

# testing expression
#minimize obj: sum{i in M} log(1+1/exp(sum{j in N}theta[j])); 

#original lgorithmic regression (no parameter)
minimize obj: sum{i in N}(theta[i]^2) + sum{i in M}log(1+1/exp(sum{j in N}theta[j])); 

#original lgorithmic regression
#minimize obj: lambda*(sum{i in N} theta[i]^2) + sum{i in M}log(1+1/exp(y[i]*(sum{j in N}theta[j]*x[i,j]) )); 


#logx->x^2
#minimize obj: lambda*(sum{i in N} theta[i]^2) + sum{i in M}(1+1/exp(y[i]*(sum{j in N}theta[j]*x[i,j]) ))^2 ; 

#exp x - > x^2
#minimize obj: lambda*(sum{i in N} theta[i]^2) + sum{i in M}log(1+1/(y[i]*(sum{j in N}theta[j]*x[i,j]) )^2 ); 

# 1/exp x -> x^2
#minimize obj: lambda*(sum{i in N} theta[i]^2) + sum{i in M}log(1+(y[i]*(sum{j in N}theta[j]*x[i,j]) )^2 ); 

