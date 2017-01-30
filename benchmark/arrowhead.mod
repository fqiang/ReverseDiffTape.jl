#arrowhead.mod

option presolve 0;
param n;
param k;
set N := {1..n};
set K := {1..k};
set J := {1..n+k};
var x{1..n+k}, <=100 , >=-100;

minimize obj: sum{i in J} (x[i]^2) + sum{i in N} ( cos( sum{j in K} x[i+j] ) + sum{j in K}( x[i] + x[j] )^2  );