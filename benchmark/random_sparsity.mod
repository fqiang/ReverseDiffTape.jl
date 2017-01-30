#random_sparsity.mod

option presolve 0;
param n;
param k;
set N = {1..n};
set K = {1..k};

var x{1..n}, <=100 , >=-100;

param rand{b0 in N, b in K};
minimize obj: sum{i in N} ( (x[i] -1)^2 + prod{j in K} x[rand[i,j]] );
