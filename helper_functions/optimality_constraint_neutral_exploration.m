function [c,ceq] = optimality_constraint_neutral_exploration(x)
ceq(1) = (x(6))^(-2)-(x(2))^(-2)-(x(14))^(-2) ;
ceq(2)= (1-1./(1+exp(-x(3)*x(4))) + 1./(1+exp(x(3)))) - (1-1./(1+exp(-x(11)*x(12))) + 1./(1+exp(x(11))));
c = [];