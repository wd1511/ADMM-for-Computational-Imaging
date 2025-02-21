function y=MyV2C2(x)

[nx,ny]=size(x);
y=zeros(nx*ny/2,1);
y=x(1:nx*ny/2).*exp(x(nx*ny/2+1:nx*ny)*1i);