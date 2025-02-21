function dyu = Dy(U)
% y-axis forward difference
dyu = U([2:end,1],:)-U;  