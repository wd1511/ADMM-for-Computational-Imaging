  function [out] = JR_ADMM(Qr,Uzr,Ur,Vr,gamma,sgn)
% this function just calculates the criterion to be minimized w.r.t. the current values of components
% -> Qr - observations, 
% -> Uzr - estimate of the sensor wave field denoted in the AL algorithm as ur, t+1/2 
% -> Vr - estimation (preliminary) of the sensor wave field denoted in the AL algorithm as ur, t+1 
% -> Lr - Lagrange multipliers
% -> n_u0est - a norm of the estimate of the object wave field
% -> factor = [gamma, mu] - only two parameters !!! why so see below
% -> sgn - a matrix of sign of the denuminator
% Vr must be real and positive
cc1 = 4/gamma*sgn.*(abs(Ur).^2-Qr)+1;
Ur = (Uzr - Vr)./cc1;
out = abs(Qr-abs(Ur).^2).^2 + gamma/2 * abs(Ur-Uzr+Vr).^2-gamma/2 * abs(Vr).^2;
