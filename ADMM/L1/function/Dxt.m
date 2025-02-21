function dxtu = Dxt(U)
% -(x-axis backward difference)
dxtu = U(:,[end,1:end-1])-U;