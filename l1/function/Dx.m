function dxu = Dx(U)
% x-axis forward difference
dxu = U(:,[2:end,1])-U;
