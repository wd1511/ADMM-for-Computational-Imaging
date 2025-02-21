function  y= TV(a,op)
%1范数的TV
%op=1,水平方向梯度
%op=2,垂直方向梯度
%op为其他，取0

[nx,ny,nz] = size(a);
if op==2
	y=circshift(a,[-1 0 0])-a;
    y(nx,:,:)=0.0;
elseif op==1
    y=circshift(a,[0 -1 0])-a;
    y(:,ny,:)=0.0;
else
    y=zeros(nx,ny,nz);
end

