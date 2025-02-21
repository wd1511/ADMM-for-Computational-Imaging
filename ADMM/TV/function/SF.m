function [Y] = SF(X,V,gamma)
[Nv,Mv,Kv] = size(V); [Nx,Mx,Kx] = size(X); 

if Nx == Nv && Mx == Mv && Kx == Kv
    Y = zeros(Nx,Mx);
    for index=1:Kx
        Y = Y + fft2(X(:,:,index)).*V(:,:,index)*gamma(index);
    end
end

if Nx < Nv && Mx < Mv && Kx == 1
    Y = zeros(Nv,Mv); XX = zeros(Nv,Mv);
    bourdery = Nv/2-Nx/2+1:Nv/2+Nx/2; bourderx = Mv/2-Mx/2+1:Mv/2+Mx/2;
    XX(bourdery,bourderx) = X;
    
    for index=1:Kv
        Y = Y + fft2(XX).*V(:,:,index)*gamma(index);
    end
end