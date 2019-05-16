function y= Dh(a,gamma0,gamma1,gamma2)
    [nx,ny]=size(a);
    
    a=a(:);
    
    %水平梯度
    a1=eye(nx*ny);
    a1=circshift(a1,[-ny,0,0]);
    a1=a1-eye(nx*ny);
    a1((nx-1)*ny+1:nx*ny,:)=0;
    
    %垂直梯度
    a2=eye(nx*ny);
    a2=circshift(a2,[-1,0,0]);
    a2=a2-eye(nx*ny);
    for index = 1:nx*ny
        if mod(index,nx)==0
            a2(index,:)=0;
        end
    end
   
    A = gamma0*eye(nx*ny)+gamma1*a1'*a1+gamma2*a2'*a2;
    y = inv(A)*a;
end

