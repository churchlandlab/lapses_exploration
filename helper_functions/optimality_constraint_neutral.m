function [c,ceq] = optimality_constraint_neutral(x)
ceq(1) = (x(6))^(-2)-(x(2))^(-2)-(x(14))^(-2) ;
ceq(2) = x(3)-x(11);
c = [];