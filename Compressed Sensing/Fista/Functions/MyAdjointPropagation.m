function eta=MyAdjointPropagation(S,E,Nx,Ny,Nz,phase)

cEsp=ifftshift(conj(ifft2(conj(real(S)))));

cEs=conj(phase).*repmat(cEsp,[1 1 Nz]);

eta=zeros(Nx,Ny,Nz);
for i=1:Nz
    eta(:,:,i)=conj(fft2(conj(ifftshift(cEs(:,:,i)))));
end

eta=conj(E).*eta;
%x=S;
%a=fftshift(fft2(x)).*phase;

%function y = f1(x) 
%    y=ifftshift(ifft2(x));
%    TV=MyTV3D_conv(((y)));
%    y=(sum(abs(TV),4));
%end
%y1 = f1(a);

%function y = f2(x) 
    %y=fftshift(fft2(x));
    %TV=MyTV3D_conv(((y)));
    %y=(sum(abs(TV),4));
%end
%y2 = f2(x);

%eta=y1.*phase;
%eta=eta.*y2;
%eta=conj(E).*eta;

%end
