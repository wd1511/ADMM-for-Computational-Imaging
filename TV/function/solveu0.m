function y = solveu0(a,gamma0,gamma1,gamma2)
    [nx,ny]=size(a);
    a1=eye(nx);
    a1=circshift(a1,[-1,0]);
    a1=a1-eye(nx);
    a1(nx,:)=0;
    
    a2=eye(ny);
    a2=circshift(a2,[0,-1]);
    a2=a2-eye(ny);
    a2(:,ny)=0;
    
    A=gamma0*eye(nx)+gamma1*a1'*a1+gamma2*a2'*a2;
    
   	y =inv(A)*a;
end

