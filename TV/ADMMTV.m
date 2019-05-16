function [out,tm,RMSE]=ADMMTV(o,u0,op,z1,dz,lambda,deltaA,deltaB,gamma,mu,a0,ITER,alpha)
% default parameters (if not specified)
dgamma = 0.1;
dmu = 5e-3;
da0 = 1;
dITER = 20;
beta = 0.1;
q = [1,1];
dalpha = 1;
method = 0;
av = true;

%gamma0是||Aru0-o||的系数，gamma1和gamma2是水平和竖直方向的TV的系数


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
if method < 2
    Nya = Nya+Ny0;Nxa = Nxa+Nx0;    
end
coordy = Nya/2-N/2+1:Nya/2+N/2; coordx = Nxa/2-M/2+1:Nxa/2+M/2; % center part of the sensor plane wave field
bourdery = Nya/2-Ny0/2+1:Nya/2+Ny0/2; bourderx = Nxa/2-Nx0/2+1:Nxa/2+Nx0/2; % center part of the object

gamma0=1.8/K;
gamma1=0.02;
gamma2=0.02;
gamma3=0.9;

D1 = abs(psf2otf([1,-1],[Nya,Nxa])).^2;
D2 = abs(psf2otf([1,-1],[Nya,Nxa])).^2;

gammar = gamma3/K*ones(K,1);
gamma_scaler = 2;

% initialization（preallocation)
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
% 1为水平方向，2为垂直方向
u1 = zeros(Nya,Nxa); 
u2 = zeros(Nya,Nxa); 
v1 = zeros(Nya,Nxa);
v2 = zeros(Nya,Nxa);
u0_iter = zeros(Nya, Nxa); a_=1/2;
u0_iter(bourdery, bourderx)=a_; 
u1(:,:)=a_;u2(:,:)=a_;
% figure(3),close(3),
% V=function_AL_DataSynthesis(u0,lambda,delta,z1,dz,K,[],1,1);
[V,B,strs,tm]=TransformMatrix(method<1,z1,dz,lambda,deltaA,deltaB,Nya,Nxa,K,gammar,mu,av);
% 0. calculation of the wave fields at the sensor planes
%mse(u0_iter,u0_iter.*V(:,:,2).*B(:,:,2))
for index = 1:K
    Ar_u0(:,:,index) = ifft2(fft2(u0_iter).*V(:,:,index));
    Ar_u0(coordy,coordx,index)=sqrt(o(:,:,index)).*exp(i*angle(Ar_u0(coordy,coordx,index)));
end
Ur = Ar_u0;

% ----------------------------- ADMM ----------------------------------------
for jndex = 1:ITER
    tic;
    %ur
    for index = 1:K
        if jndex>1
            Ar_u0(:,:,index)=ifft2(fft2(u0_iter).*V(:,:,index));
            Ur(:,:,index) = Ar_u0(:,:,index);
        end
        % minimization on ur
        %lllll = (Ar_u0(:,:,index)-u0_iter).^2;
        %lll=sum(sum(abs(lllll),1),2)/Nya/Nxa;
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
            %u0
            PP = SF(Ur+Vr,B,gammar*gamma_scaler); %PP=gammar*B*(Ur+Vr)
            PP1 = ifft2(PP + gamma1 * TV(u1+v1,1)+ gamma2 * TV(u2+v2,2)); %u0_iter = PP + gamma*(u+v)
            
            %PP2 = 0.5*gamma0+0.5*gamma0 * sum(abs(V).^2,3)+gamma1 * D1 +gamma2 * D2;
            %PP2 = 0.5*gamma0 * sum(abs(V).^2,3)+gamma1 * D1 + gamma2 * D2;
            PP2 = gamma0 * sum(abs(V).^2,3) + gamma1 * D1 +gamma2 * D2;
            u0_iter = ifft2(1./PP2.*fft2(PP1));
            %u0_iter = solveu0(u0_iter,gamma0,gamma1,gamma2);
            %u1，u2
            u1(bourdery,bourderx) = soft(u0_iter(bourdery,bourderx)-v1(bourdery,bourderx),mu,gamma1, 0);
            u2(bourdery,bourderx) = soft(u0_iter(bourdery,bourderx)-v2(bourdery,bourderx),mu,gamma2, 0);
            
            out = (u0_iter(bourdery,bourderx));
            u1_out = u1(bourdery,bourderx);
            u2_out = u2(bourdery,bourderx);
        otherwise
            disp('Error: Wrong setting of the operation, please choose op={0,1,2,3}.'); out=[];t=0; return     
    end
    u0_iter = zeros(Nya,Nxa); 
    u1 = zeros(Nya, Nxa);
    u2 = zeros(Nya, Nxa);
    u0_iter(bourdery,bourderx) = out;
    u1(bourdery, bourderx) = u1_out;
    u2(bourdery, bourderx) = u2_out;
    tm = tm + toc;
    
    %更新vr和v
    if ~mod(method,2)
        for kndex=1:K
            Vr(:,:,kndex) = Vr(:,:,kndex) - alpha * (ifft2(fft2(u0_iter).*V(:,:,kndex)) - Ur(:,:,kndex));
        end
        v1 = v1 - alpha * (TV(u0_iter,1) - u1);
        v2 = v2 - alpha * (TV(u0_iter,2) - u2);
        tm = tm + toc;
    end
    if cRMSE
        figure(5)
         v3 = abs(out) - abs(u0);  RMSE(1,jndex)=sqrt(mean(mean( ( v3 ).^2 )));
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
    gamma1 = gamma1 / 1.04;
    gamma2 = gamma2 / 1.04;
end
fprintf('elapsed time = %8.4f seconds\n',tm);
