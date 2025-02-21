function [out,tm,RMSE]=ADMM(o,u0,op,z1,dz,lambda,deltaA,deltaB,gamma,mu,a0,ITER,alpha)
% default parameters (if not specified)
dgamma = 0.1;
dmu = 5e-3;
da0 = 1;
dITER = 100;
beta = 0.1;
q = [1,1];
dalpha = 1;
method = 0;
av = true;

i = sqrt(-1);RMSE = []; tm = 0;
if ~isreal(op) || ~isnumeric(op) || ceil(op)-floor(op)==1 
    disp('Error: Wrong setting of the operation, please choose op={0,1,2,3}.'); out = []; return 
end

[N,M,K] = size(o);                % the size of the input set of the (noisy) intensity observations
if mod(N,2)~=0 || mod(M,2)~=0
   disp('Error: The size of the intensity observations must be even integers'); out = []; return 
end
% if K<2
%     disp('Error: It must be at least 2 intensity observations, otherwise phase retrieval doesnt make sense'); out = []; return 
% end
if (nargin<9 || isempty(gamma)),        gamma=dgamma;  end
if (nargin<14 || isempty(alpha)),       alpha=dalpha; end
if (nargin<11 || isempty(a0)), a0=da0; end
if any(~isreal(gamma)) || any(~isnumeric(gamma)) || any(gamma<=0)
    disp('Error: "gamma" must be a vector of K positive real values'); out = []; return 
end
if (nargin<10 || isempty(mu)),    mu=dmu; end
if ~isreal(mu) || ~isnumeric(mu) || mu<0
    disp('Error: "mu" must be a non-negative real value'); out = []; return 
end
if (nargin<12 || isempty(ITER)),  ITER=dITER; end
if q(1) < 1 || q(2) < 1, disp('Error: q must be not smaller than 1'); out = []; return ; end
Nya = ceil(q(1)*N/2)*2;  Nxa = ceil(q(2)*M/2)*2;     % the support of the computed Fourier images
cRMSE = false;                   % control parameter which control whether RMSE should be calculated
                                 % it the true u0 is given - RMSE is computed and the output variable 'out' is the calculated object 
                                 % otherwise - RMSE=[] is empty variable and 'out' is the calculated volumetric sensor plane wave field
if ~isempty(u0)                  % the input object is used to estimate the reconstruction accuracy in RMSE
    [Ny0,Nx0,~] = size(u0);  % the size of the input object
    if Ny0>Nya || Nx0>Nxa
        disp('Error: The size of the object "u0" should not be larger than the size of the intensity observations!');out = []; return;
    end
    if mod(Ny0,2)~=0 || mod(Nx0,2)~=0
        disp('Error: The size of the object "u0" must be even integers'); out = []; return 
    end
    cRMSE = true; RMSE = zeros(2,ITER); 
    
end 
% IF mod(op,2)=0, i.e. op IS EVEN (op={2,4,6,8...}) THEN WE USE F-DDT, OTHERWISE - ASD
% if method < 2
    Nya = Nya+Ny0;Nxa = Nxa+Nx0;    
% end
coordy = Nya/2-N/2+1:Nya/2+N/2; coordx = Nxa/2-M/2+1:Nxa/2+M/2; % center part of the sensor plane wave field
bourdery = Nya/2-Ny0/2+1:Nya/2+Ny0/2; bourderx = Nxa/2-Nx0/2+1:Nxa/2+Nx0/2; % center part of the object
gammar = (1-gamma)/K *ones(K,1);
gamma_scaler = 2;
% initialization£¨preallocation)
Ar_u0 = zeros(Nya, Nxa, K);
% INITIAL GUESS
if mod(method,1)
    Vr = 0;
    str = 'simple,';
else
    Vr = zeros(Nya, Nxa, K);
    str = 'opt';
end

% initial guess for the object: u0[k] = 1/2
u = zeros(Nya,Nxa); v = zeros(Nya,Nxa);
u0_iter = zeros(Nya, Nxa); a_=1/2;
u0_iter(bourdery, bourderx)=a_; u(:,:)=a_;
% figure(3),close(3),
% V=function_AL_DataSynthesis(u0,lambda,delta,z1,dz,K,[],1,1);
[V,B,strs,tm]=TransformMatrix(method<1,z1,dz,lambda,deltaA,deltaB,Nya,Nxa,K,gammar,mu,av);

% 0. calculation of the wave fields at the sensor planes
for index = 1:K
    Ar_u0(:,:,index) = ifft2(fft2(u0_iter).*V(:,:,index));
    Ar_u0(coordy,coordx,index)=sqrt(o(:,:,index)).*exp(i*angle(Ar_u0(coordy,coordx,index)));
end
Ur = Ar_u0;

% ----------------------------- ADMM ----------------------------------------
for jndex = 1:ITER
    tic;
    for index = 1:K
        if jndex>1
            Ar_u0(:,:,index)=ifft2(fft2(u0_iter).*V(:,:,index));
            Ur(:,:,index) = Ar_u0(:,:,index);
        end
        % minimization on ur
        if ~mod(method, 2) % fitting on the observation / chosen w.r.t. the min value of the ADMM criterion
            [cc] = cubic_equation_ADMM(o(:,:,index), Ar_u0(coordy,coordx,index),Vr(coordy,coordx,index), gamma_scaler * gammar(index));
            Ur(coordy,coordx,index) = (Ar_u0(coordy,coordx,index) - Vr(coordy,coordx,index)) ./ cc;
        else % keep the amplitude, change the phase as in the conventional Gerchberg-Saxton-Fienup algorithm
            Ur(coordy,coordx,index) = Ar_u0(coordy,coordx,index) ./ abs(Ar_u0(coordy,coordx,index)) .* sqrt(o(:,:,index));
        end
    end
    % minimization on u0
    switch floor(op)
        case 0 % complex-valued
            if jndex == 1, opp='complex'; end
            PP = SF(Ur+Vr,B,gammar);
            u0_iter = ifft2(PP) + gamma * (u+v);
            u(bourdery,bourderx) = soft(u0_iter(bourdery,bourderx)-v(bourdery,bourderx),mu,gamma_scaler * gamma, 0);
            out = (u0_iter(bourdery,bourderx));
            u_out = u(bourdery,bourderx);
        case 1 % AM
            if jndex == 1, opp='AM'; end
            PP = SF(Ur+Vr,B,gammar);
            u0_iter = ifft2(PP) + gamma * (u+v);
            u(bourdery,bourderx) = soft(u0_iter(bourdery,bourderx)-v(bourdery,bourderx),mu,gamma_sclaer * gamma, 0);
            out = real(u0_iter(bourdery,bourderx));
            u_out = u(bourdery,bourderx);
            out = out.*(out>0)+(out<0).*0.0001;% positive projection
        case 2 % PM with known amplitude a0, gradient descent AL
            if jindex == 1, opp=' OM,known a0'; end
            ph = angle(u0_iter(bourdery,bourderx));
            PP = SF(Ur+Vr,conj(V), gammar);
            for kndex=1:JTER
                XX = SF(a0*exp(i*ph),conj(V).*V, gamma);XX=ifft2(PP-XX);
                YY = 2*imag(conj(conj(a0*exp(i*ph)).*XX(coordy,coordx)));
                ph = ph - beta*YY;
            end
            out = a0*exp(i*ph);
        case 3     % PM with unknown amplitude a_ , gradient descent AL
             if jndex == 1, opp = ' PM, unknown a0'; end
             ph = angle(u0_iter(bourdery,bourderx)); 
             PP = SF(Ur+Vr,conj(V),gamma); 
             for kndex=1:JTER
                 XX = SF(a_(end)*exp(i*ph),conj(V).*V,gamma);XX = ifft2(PP - XX);
                 YY = 2*imag(conj( conj(a_(end)*exp(i*ph)) .*XX(coordy,coordx) ));
                 ph = ph - beta*YY;
             end

             % recalculation of the constant amplitude
             S = 0; SS=0; ph = ph-mean(mean(ph)); u0est = exp(i*angle(u0_iter)); u0_iter(bourdery,bourderx) = exp(i*ph);
             PP = Izk_iter+LAMBDA;
             for index=1:K                 
                 %tmp_ = ifft2(RightHandSideTerm(u0est,1,conj(V(:,:,index))));
                 XX = ifft2(fft2(u0_iter).*V(:,:,index));                 
                 SS = SS + real( sum(sum( conj(XX(coordy,coordx)) .* PP(coordy,coordx,index) )) ).*gammar(index);
                 S = S + sum(sum( conj(XX(coordy,coordx)) .* XX(coordy,coordx) )).*gammar(index);
             end
             a_=[a_ SS/(S+mu*N*M)];
             out = a_(end)*exp(i*ph);
        otherwise
            disp('Error: Wrong setting of the operation, please choose op={0,1,2,3}.'); out=[];t=0; return     
    end
    u0_iter = zeros(Nya,Nxa); u = zeros(Nya, Nxa);
    u0_iter(bourdery,bourderx) = out;
    u(bourdery, bourderx) = u_out;
    tm = tm + toc;
    
    if ~mod(method,2)
        for kndex=1:K
            Vr(:,:,kndex) = Vr(:,:,kndex) - alpha * (ifft2(fft2(u0_iter).*V(:,:,kndex)) - Ur(:,:,kndex));
        end
        v = v - alpha * (u0_iter - u);
        tm = tm + toc;
    end
    if cRMSE
        figure(5)
         v1 = abs(out) - abs(u0);  RMSE(1,jndex)=sqrt(mean(mean( ( v1 ).^2 )));
         x1 = angle(out) - angle(u0);   txmp = mean(x1(:)); x1 = x1-txmp;        RMSE(2,jndex)=sqrt(mean(mean( ( x1 ).^2 )));   
         subplot(2,3,1),    imshow(abs(out),[]), title('amplitude, ADMM'),xlabel(['RMSE=' num2str(RMSE(1,jndex))])
         subplot(2,3,2),    plot(RMSE(1,1:jndex)),  xlabel(['\gamma=' num2str(mean(gamma)) ',\alpha_r=' num2str(alpha) ',\mu=' num2str(mu)]) 
         if length(a_)>1
             subplot(2,3,2),hold on,plot(a_,'r--'),hold off
             legend('error','est. ampl.','Location','Best')
             axis([1,max(jndex,2),min(min(RMSE(1,1:jndex)),min(a_)),max(max(RMSE(1,1:jndex)),max(a_))]),title('convergence: amplitude'),ylabel('RMSE, a_0^{est}'),grid on
         else
             axis([1,max(jndex,2),min(min(RMSE(1,1:jndex)),0),max(max(RMSE(1,1:jndex)),eps)]),title('convergence: amplitude'),ylabel('RMSE'),grid on
         end
         subplot(2,3,4),imshow((angle(out)-txmp),[]), title('phase, ADMM'), xlabel(['RMSE=' num2str(RMSE(2,jndex))])
         subplot(2,3,5),    plot(RMSE(2,1:jndex)), xlabel(['K=' num2str(K) ',iter=' num2str(jndex) '/' num2str(ITER)]),title('convergence: phase'),ylabel('RMSE'), grid on
         if jndex>1
             axis([1,max(jndex,2),min(min(RMSE(2,1:jndex)),0),max(max(RMSE(2,1:jndex)),eps)])
         end
         subplot(2,3,3),    plot(1:Nx0,abs(out(Ny0/2+1,:)),'b',1:Nx0,abs(u0(Ny0/2+1,:)),'r')
         if length(a_)>1
             xlabel(['a_0^{true}=' num2str(a0) ' , a_0^{est}=' num2str(a_(end))])
         end
         axis([1,M,0,max( max(abs(out(Ny0/2+1,:))), max(abs(u0(Ny0/2+1,:))) )]),title(['amplitude: cross-section along the ' num2str(Ny0/2+1) '-th line'])
         subplot(2,3,6),    plot(1:Nx0,angle(out(Ny0/2+1,:))-txmp,'b',1:Nx0,angle(u0(Ny0/2+1,:)),'r'),title(['phase: cross-section along the ' num2str(Ny0/2+1) '-th line'])
         axis([1,M,min(min((angle(out(Ny0/2+1,:))-txmp)),min(angle(u0(Ny0/2+1,:)))),max([max((angle(out(Ny0/2+1,:))-txmp)), max(angle(u0(Ny0/2+1,:))) ,eps])]),xlabel([strs str opp])
         drawnow 
    end
end
fprintf('elapsed time = %8.4f seconds\n',tm);

function [F,B,strs,t]=TransformMatrix(x,z1,dz,lambda,deltaA,deltaB,N,M,K,gamma,mu,av)
% calculation of the F-DDT or ASD transfer functions for the forward (F)
% and backward (B) diffraction propagation

F = zeros(N,M,K);  
if x,  QQ = zeros(N,M); strs = 'F-DDT '; 
else   QQ = 0;          strs = 'ASD ';
end
t = 0;    
for index = 1:K
    if x       
        [F(:,:,index),tmp] = TransferFunctionFDDT(z1+(index-1)*dz,lambda,deltaA,deltaB,N,M,av);         %function_TransferFunction_FDDT(z1+(index-1)*dz,lambda,deltaA,N,M,av);
    else
        [F(:,:,index),tmp] = TransferFunctionASD(z1+(index-1)*dz,lambda,deltaA,deltaB,N,M);
    end; t = t+tmp;
end
tic
for index = 1:K
    if x, QQ = QQ + abs(F(:,:,index)).^2 * gamma(index);
    else  QQ = QQ + gamma(index);
    end
end
QQ = QQ + mu;
B = zeros(N,M,K);
for index = 1:K
    B(:,:,index) = conj(F(:,:,index));
    %B(:,:,index) = conj(F(:,:,index)) ./ QQ;
end
t = t + toc;


function [Y] = SF(X,V,gamma)
[Nv,Mv,Kv] = size(V); [Nx,Mx,Kx] = size(X); 

if Nx == Nv && Mx == Mv && Kx == Kv
    Y = zeros(Nx,Mx);
    for index=1:Kx
        Y = Y + fft2(X(:,:,index)).*V(:,:,index)*gamma(index);
    end
end

if Nx < Nv && Mx < Mv && Kx == 1
    Y = zeros(Nv,Mv); XX = zeros(Nv,Mv);
    bourdery = Nv/2-Nx/2+1:Nv/2+Nx/2; bourderx = Mv/2-Mx/2+1:Mv/2+Mx/2;
    XX(bourdery,bourderx) = X;
    
    for index=1:Kv
        Y = Y + fft2(XX).*V(:,:,index)*gamma(index);
    end
end


  
function [out] = JR_ADMM(Qr,Uzr,Ur,Vr,gamma,sgn)
% this function just calculates the criterion to be minimized w.r.t. the current values of components
% -> Qr - observations, 
% -> Uzr - estimate of the sensor wave field denoted in the AL algorithm as ur, t+1/2 
% -> Vr - estimation (preliminary) of the sensor wave field denoted in the AL algorithm as ur, t+1 
% -> Lr - Lagrange multipliers
% -> n_u0est - a norm of the estimate of the object wave field
% -> factor = [gamma, mu] - only two parameters !!! why so see below
% -> sgn - a matrix of sign of the denuminator
% Vr must be real and positive
cc1 = 4/gamma*sgn.*(abs(Ur).^2-Qr)+1;
Ur = (Uzr - Vr)./cc1;
out = abs(Qr-abs(Ur).^2).^2 + gamma/2 * abs(Ur-Uzr+Vr).^2-gamma/2 * abs(Vr).^2;


function u = soft(s, mu, gamma,sep)
sa = real(s); sb=imag(s);
if sep
    ua = sign(sa).*max(abs(sa)-mu/gamma,0);
    ub = sign(sb).*max(abs(sb)-mu/gamma,0);
else
    coef = sqrt(sa.^2./(sa.^2+sb.^2));
    ua = coef.*sign(sa).* max(sqrt(sa.^2+sb.^2)-mu/gamma,0);
    ub = sb ./ sa .* ua;
end
u = ua + ub * sqrt(-1);



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


