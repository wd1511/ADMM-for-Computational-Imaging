function TV=tvlpls(x)

[nx,ny]=size(x);
TV=zeros(nx,ny,1,4);

eo=circshift(x,[-1 0])-x;
eo(nx,:)=0.0;

wo=circshift(x,[1 0])-x;
wo(1,:)=0.0;

so=circshift(x,[0 -1])-x;
so(:,ny)=0.0;
so(:,ny/2)=0.0;

no=circshift(x,[0 1])-x;
no(:,1)=0.0;
no(:,ny/2+1)=0.0;

nese=circshift(no,[-1 0])-circshift(so,[-1 0]);
nese(nx,:,:)=0.0;

nwsw=circshift(no,[1 0])-circshift(so,[1 0]);
nwsw(1,:,:)=0.0;

nenw=circshift(eo,[0 1])-circshift(wo,[0 1]);
nenw(:,1)=0.0;
nenw(:,ny/2+1)=0.0;

sesw=circshift(eo,[0 -1])-circshift(wo,[0 -1]);
sesw(:,ny)=0.0;
sesw(:,ny/2)=0.0;

v1=var(x(:));

%TV(:,:,1)=sqrt(eo.*eo + (nese+no-so).*(nese+no-so)/16).*exp(-eo.*eo/v1/v1);
%TV(:,:,2)=sqrt(wo.*wo + (nwsw+no-so).*(nwsw+no-so)/16).*exp(-wo.*wo/v1/v1);
%TV(:,:,3)=sqrt(no.*no + (nenw+eo-wo).*(nenw+eo-wo)/16).*exp(-no.*no/v1/v1);
%TV(:,:,4)=sqrt(so.*so + (sesw+eo-wo).*(sesw+eo-wo)/16).*exp(-so.*so/v1/v1);

TV(:,:,1,1)=eo.*exp(-eo.*eo/v1/v1) + (nese+no-so).*exp(-(nese+no-so).*(nese+no-so)/v1/v1)/4;
TV(:,:,1,2)=wo.*exp(-wo.*wo/v1/v1) + (nwsw+no-so).*exp(-(nwsw+no-so).*(nwsw+no-so)/v1/v1)/4;
TV(:,:,1,3)=no.*exp(-no.*no/v1/v1) + (nenw+eo-wo).*exp(-(nenw+eo-wo).*(nenw+eo-wo)/v1/v1)/4;
TV(:,:,1,4)=so.*exp(-so.*so/v1/v1) + (sesw+eo-wo).*exp(-(sesw+eo-wo).*(sesw+eo-wo)/v1/v1)/4;

