function [F,B,strs,t]=TransformMatrix(x,z1,dz,lambda,deltaA,deltaB,N,M,K,gamma,mu,av)
% calculation of the F-DDT or ASD transfer functions for the forward (F)
% and backward (B) diffraction propagation

F = zeros(N,M,K);  
if x,  QQ = zeros(N,M); strs = 'F-DDT '; 
else   QQ = 0;          strs = 'ASD ';
end
t = 0;    
for index = 1:K
    if x       
        [F(:,:,index),tmp] = TransferFunctionFDDT(z1+(index-1)*dz,lambda,deltaA,deltaB,N,M,av);         %function_TransferFunction_FDDT(z1+(index-1)*dz,lambda,deltaA,N,M,av);
    else
        [F(:,:,index),tmp] = TransferFunctionASD(z1+(index-1)*dz,lambda,deltaA,deltaB,N,M);
    end; t = t+tmp;
end
tic
for index = 1:K
    if x, QQ = QQ + abs(F(:,:,index)).^2 * gamma(index);
    else  QQ = QQ + gamma(index);
    end
end
QQ = QQ + mu;
B = zeros(N,M,K);
for index = 1:K
    B(:,:,index) = conj(F(:,:,index));
    %B(:,:,index) = conj(F(:,:,index)) ./ QQ;
end
t = t + toc;