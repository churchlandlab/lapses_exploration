function [c,ceq] = inactivation_optimality_constraint(x)
% ceq(1)= (x(6)).^(-2)-(x(2)).^(-2)-(x(10)).^(-2) ;
% ceq(2)= (x(18)).^(-2)-(x(14)).^(-2)-(x(22)).^(-2) ;
ceq(1)= (x(7)).^(-2)-(x(2)).^(-2)-(x(12)).^(-2) ;
ceq(2)= (x(22)).^(-2)-(x(17)).^(-2)-(x(27)).^(-2) ;
%Right inactivation
% ceq(3) = x(15)/x(3) - x(19)/x(7);
% ceq(4) = x(23)/x(11)-x(19)/x(7);
%Left inactivation
% ceq(3) = x(16)/x(4) - x(20)/x(8);
% ceq(4) = x(24)/x(12)-x(20)/x(8);
c = [];