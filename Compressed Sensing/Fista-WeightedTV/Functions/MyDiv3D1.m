function y=MyDiv3D1(TV,v1)

n=size(TV);

x_shift=circshift(TV(:,:,:,1),[1 0 0 0]);
yx=TV(:,:,:,1)-x_shift;
yx(1,:,:)=TV(1,:,:,1);
yx(n(1),:,:)=-x_shift(n(1),:,:);
yx1=yx.*exp(-yx.*yx)./max(abs(yx),0.01);

y_shift=circshift(TV(:,:,:,2),[0 1 0 0]);
yy=TV(:,:,:,2)-y_shift;
yy(:,1,:)=TV(:,1,:,2);
yy(:,n(2),:)=-y_shift(:,n(2),:);
yy1=yy.*exp(-yy.*yy)./max(abs(yy),0.01);

z_shift=circshift(TV(:,:,:,3),[0 0 1 0]);
yz=TV(:,:,:,3)-z_shift;
yz(:,:,1)=TV(:,:,1,3);
yz(:,:,n(3))=-z_shift(:,:,n(3));
yz1=yz.*exp(-yz.*yz)./max(abs(yz),0.01);

y=yx1+yy1+yz1;
