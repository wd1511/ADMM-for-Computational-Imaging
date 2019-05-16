function dytu = Dyt(U)
% -(y-axis backward difference)
dytu = U([end,1:end-1],:)-U;