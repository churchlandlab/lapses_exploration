function [c,ceq] = optimality_constraint(x)
ceq= (x(6))^(-2)-(x(2))^(-2)-(x(10))^(-2) ;
c = [];