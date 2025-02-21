function y=MyTVphi(x,Nvx,Nvy,Nvz)

%%xx = phase_unwrap(angle(x));

X=reshape(x,Nvx,Nvy,Nvz);

v1=var(x);

%[y1,dif1]=MyTVnorm((X));
[y2,dif2]=MyTVnorm1((X),v1);

y = y2;
%v1 = var(x);
%a = 0;
%m=Nvx*Nvy*Nvz;
%for i = 1 : m
%    for j = 1 : m
%        a = a +abs(x(i)-x(j))*exp(-(x(i)-x(j))*(x(i)-x(j))/v1/v1);
%    end;
%end;
%y=a;

