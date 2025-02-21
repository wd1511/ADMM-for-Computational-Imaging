function TV=MyTV3D_conv1(x,v1)

%[Nvx,Nvy,Nvz]=size(x);
%x=x(:);
%x=MyV2C(x);
%x=reshape(x,Nvx,Nvy,Nvz);

[nx,ny,nz]=size(x);
TV=zeros(nx,ny,nz,3);
TV(:,:,:,1)=(circshift(x,[-1 0 0])-x).*exp(-(circshift(x,[-1 0 0])-x).*(circshift(x,[-1 0 0])-x)/v1/v1);
TV(nx,:,:,1)=0.0;

TV(:,:,:,2)=(circshift(x,[0 -1 0])-x).*exp(-(circshift(x,[0 -1 0])-x).*(circshift(x,[0 -1 0])-x)/v1/v1);
TV(:,ny,:,2)=0.0;
TV(:,ny/2,:,2)=0.0;

TV(:,:,:,3)=(circshift(x,[0 0 -1])-x).*exp(-(circshift(x,[0 0 -1])-x).*(circshift(x,[0 0 -1])-x)/v1/v1);
TV(:,:,nz,3)=0.0;
TV(:,:,:,3)=TV(:,:,:,3).*(1.0);


