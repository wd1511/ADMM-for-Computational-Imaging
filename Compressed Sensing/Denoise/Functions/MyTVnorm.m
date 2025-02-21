function [y,dif]=MyTVnorm(x)

%TV=MyTV3D_conv(((x)));

%dif=(sum(abs(TV),4));

TV = tvlpls(((x)));
dif=(sum(abs(TV),3));
y=sum(dif(:));
