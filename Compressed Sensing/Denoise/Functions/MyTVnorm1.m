function [y,dif]=MyTVnorm1(x,v1)

TV=MyTV3D_conv1(((x)),v1);

%dif=(sum(abs(TV),4));

%TV = tvlpls(((x)));
dif=(sum(abs(TV),3));
y=sum(dif(:));