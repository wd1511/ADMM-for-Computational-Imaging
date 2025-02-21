function cc = cubic_equation_ADMM(Qr,Uzr,Vr,gamma)

% A.*X.^3 + B.*X.^2 +C.*X + D = 0 - general cubic equation, change the variables: X -> Y=X+B./3/A  (Tschirnhaus transformation), then
% P = (3*A.*C-B.B) ./ (3*A.A), Q = (2*B.^3 - 9*A*B*C + 27*A*A*D) ./ (27*A.^3) and we are forced to conclude that Y.^3 + P.*Y + Q = 0, 
% and solve it with the method due to Scipione del Ferro and Tartaglia, published by Gerolamo Cardano in 1545.

P = gamma/4-Qr;  
Q = -gamma/4.*abs(Uzr-Vr);  % # 1 cc>0

i = sqrt(-1);
% DD is the discriminant:
DD = (P/3).^3+(Q/2).^2;
% DD>0 -> 1Re, 2complex-conjugate (X1 or X3 only !!!)
% DD<0 -> 3Re
% DD=0 -> 2Re same (3Re same if P=Q=0)

% 1. solving the uncertainty problem: 0/0, x/0, 0/x
[NN,MM] = size(DD);
if 1    % FAST:  Cardano solution via assumtion of real u, v
    tmp = -Q/2+sqrt(DD);  u = zeros(NN,MM);  
    Y1 = imag(tmp)==0;
    u(Y1) = sign(tmp(Y1)).*(abs(tmp(Y1)).^(1/3));
    u(~Y1) = tmp(~Y1).^(1/3);

    msk = (real(u) > 1e-6 & imag(u) > 1e-6);
    X1 = zeros(NN,MM);X2 = X1;X3 = X1;
    
    if any(~msk(:))
        nmsk = ~msk;
        % v must be real because of Cardano assumption
        tmp = -Q(nmsk)/2-sqrt(DD(nmsk));
        Y2 = imag(tmp)==0;
        v = tmp*0;
        v(Y2) = sign(tmp(Y2)).*(abs(tmp(Y2)).^(1/3));
        v(~Y2) = tmp(~Y2).^(1/3);

        X1(nmsk)=(u(nmsk)+v);
        w = i*sqrt(3)/2;
        X2(nmsk)=-X1(nmsk)/2+w*(u(nmsk)-v);
        X3(nmsk)=-X1(nmsk)/2-w*(u(nmsk)-v);    
    end
        
    z = P(msk)/3./u(msk);
    w = sqrt(3)/2*i*(u(msk)+z);
    X1(msk)=(u(msk)-z);
    X2(msk)=-X1(msk)/2-w; % complex conjugate to X1 or X3  
    X3(msk)=-X1(msk)/2+w;
    
    
else % SLOW, BUT MORE ACCURATE for problem points??? 
    for index = 1:NN
        for jndex = 1:MM
            rr = roots([1 0 P(index,jndex) Q(index,jndex)]);
            X1(index,jndex) = rr(1);
            if length(rr)>=2
                X2(index,jndex) = rr(2);
            end
            if length(rr)==3
                X3(index,jndex) = rr(3);    
            end    
        end
    end
    
end
    
msk = (DD>0);   % 1Re

rX1 = 10000*ones(NN,MM);Y1 = (msk & imag(X2)==imag(-X3)) | (~msk & imag(X1)<1e-10);rX1(Y1) = abs(real(X1(Y1)));
rX2 = 10000*ones(NN,MM);Y2 = (msk & imag(X1)==imag(-X3)) | (~msk & imag(X2)<1e-10);rX2(Y2) = abs(real(X2(Y2)));
rX3 = 10000*ones(NN,MM);Y3 = (msk & imag(X1)==imag(-X2)) | (~msk & imag(X3)<1e-10);rX3(Y3) = abs(real(X3(Y3)));

clear P Q DD


JJtmp = 10000*ones(NN,MM);JJtmp(Y1) = JR_ADMM(Qr(Y1),Uzr(Y1),rX1(Y1),Vr(Y1),gamma,sign(X1(Y1)));
JJ(:,:,1)=JJtmp;JJtmp = 10000*ones(NN,MM);JJtmp(Y2) = JR_ADMM(Qr(Y2),Uzr(Y2),rX2(Y2),Vr(Y2),gamma,sign(X2(Y2)));
JJ(:,:,2)=JJtmp;JJtmp = 10000*ones(NN,MM);JJtmp(Y3) = JR_ADMM(Qr(Y3),Uzr(Y3),rX3(Y3),Vr(Y3),gamma,sign(X3(Y3)));
JJ(:,:,3)=JJtmp;minJJ3 = min(JJ,[],3);rrY = double( minJJ3 == JJ(:,:,1)  ).*rX1+double(  minJJ3 == JJ(:,:,2)  ).*rX2+double(  minJJ3 == JJ(:,:,3)   ).*rX3;
cc = 4/gamma *(abs(rrY).^2-Qr)+1;

