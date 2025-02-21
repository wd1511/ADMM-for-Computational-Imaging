function y=MyTVphi(x,Nvx,Nvy,Nvz)

%%xx = phase_unwrap(angle(x));

X=reshape(x,Nvx,Nvy,Nvz);

v1=var(x);

[y2,dif2]=MyTVnorm1((X),v1);

y = y2;


