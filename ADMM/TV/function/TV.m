function  y= TV(a,op)
%1������TV
%op=1,ˮƽ�����ݶ�
%op=2,��ֱ�����ݶ�
%opΪ������ȡ0

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

