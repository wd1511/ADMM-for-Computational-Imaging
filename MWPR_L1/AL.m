function [out,tm,RMSE]=AL(o,u0,op,z1,dz,lambda,deltaA,deltaB,gamma,mufactor,a0,ITER,alpha)
% AL() realizes the parallel augmented Lagrangian (AL) phase-retrieval algorithm: 
% the reconstruction of a complex-valued object from (noisy) intensity observations
% as it is described in [1]. The diffraction wave field propagation is performed by 
% default with the double size F-DDT.
%
% Notation:
% [out,t,RMSE] = AL(o,u0,op,z1,dz,lambda,deltaA,deltaB,gamma,mufactor,a0,ITER,q,alpha)
% ------------------------------- input -----------------------------------
% o                     3D matrix of the (noisy) intensity observations of the size NxM 
%                       for K sensor planes, i.e. the size of the matrix must be NxMxK,
%                       where N, M and K are positive even integers, K>1
% u0                    true object of the size NxM or smaller used to calculate the 
%                       reconstruction accuracy (and algorithm convergence) in RMSE
% op                    this parameter determines what method of AL to use:
%                       If op = 0 - the original AL in case the object structure is unknown
%                       If op = 1 - the AL-A to reconstruct the amplitude-only object (AM)
%                       If op = 2 - the gradient descent AL-Ph to reconstruct the phase-only object (PM) with the known constant amplitude
%                       If op = 3 - the gradient descent AL-Ph to reconstruct the phase-only object with the unknown constant amplitude
% z1                    distance between the object and sensor plane [m]
% dz                    distance between the observation planes [m]
% lambda                wavelength [m]
% deltaA                pitch size w.r.t. rows of the discrete object/sensor plane wave field [m]
% deltaB                pitch size w.r.t. columns of the discrete object/sensor plane wave field [m]
% gamma                 regularization parameter for the AL penalties 
% mufactor              regularization parameter for the diffraction propagation operators
% a0                    the true constant amplitude in case of PM (by default a0=1)
% ITER                  number of iterations of the SBMIR algorithm (optional, by default ITER=100)
% alpha                 updating step of the Lagrange multipliers (optional, alpha is from [0,1], by default alpha=1)
% ------------------------------- output ----------------------------------
% out                   If u0=[] is not defined, the output is the NxMxK complex-valued matrix of the 
%                       reconstructed volumetric wave field in K sesor planes
%                       Otherwise 'out' is the reconstructed object of the size of the input 'u0' 
% t                     the computational time in seconds (optional)
% RMSE                  the reconstruction accuracy (via RMSE) at all ITER iteration of the AL
%                       phase retrieval algorithm w.r.t. the object amplitude (RMSE(1,:)) and object phase (RMSE(2,:))    
%
% Example:
% Reconstruction of a complex-valued object from K=5 noiseless observations
% by AL (no prior information on the object type).
%
% u0 = im2double(imread('Lena256.png')); u0 = u0+0.1;    u0 = u0/max(u0(:));
% [N,M]=size(u0);                             % object size
% Nz = N; Mz = M;                             % sesnor size (be default, the same as the object size)
% lambda = 532e-9;                            % wavelength [m]
% delta = 6.7e-6;                             % square pixel size with 100% fill factor [m]
% ITER = 100;                                 % the number of iterations
% z1 =   0.0432;                              % distance to the first sensor plane
% dz = 2e-3;                                  % distance between measurement planes [m]
% K = 5;                                      % number of measurements
% sigma = 0*ones(K,1);                        % noise level of the intnesity observations
% gamma = 10*ones(K,1);                       % penalty coefficient of the AL base algorithms
% mufactor = 0.0005;                          % regularization parameter
% uz = zeros(Nz,Mz,K); o = zeros(Nz,Mz,K);noise = randn(Nz,Mz,K);
% for index = 1:K
%     uz(:,:,index) = FDDT(u0,0,z1+(index-1)*dz,lambda,delta,delta,[1,1]);
%     o(:,:,index)= abs(uz(:,:,index)).^2 + sigma(index)*noise(:,:,index);
%     o(:,:,index)=o(:,:,index).*(o(:,:,index)>=0) + 0.0001*(o(:,:,index)<0); % positive projection to avoid negative or zero intensity
% end
% op = 0;                                     % reconstruction of a complex-valued object
% [u0est,t,RMSE] = AL(o,u0,op,z1,dz,lambda,delta,delta,gamma,mufactor,[],ITER);
%
% Reference:            [1] A. Migukin, V. Katkovnik, and J. Astola, J. Opt. Soc. Am. 28, 993-1002 (2011). 
% Authors:              Artem Migukin
% Software version:     2.0, January 2013.
% Required functions:   TransferFunctionFDDT(), TransferFunctionASD()

% default parameters (if not specified)
dgamma = 10;         % coefficients (all the same) for the augmented Lagrangian penalties
dmufactor = 1e-2;    % regularization parameter for the transfer functions
da0 = 1;             % true constant object amplitude in case of PM
dITER = 100;         % number of iteration of the phase-retrieval algorithm
JTER = 40;           % number of iterations for the gradient descent AL-Ph   
beta = 0.1;          % updating step of the gradient descent AL-Ph
q = [1,1];           % the redundancy of the computantional size w.r.t. the y- and x-direction, qx>=1 and qy>=1 
                     % (by default ASD is used of the same size as the size of the object/sensor plane wave field dq=[1,1],
                     % for instance, for realization of ASD via double size zero-padding , please set dq = [2,2].)
                     % By default F-DDT uses the double size support.
dalpha = 1;          % updating step of the Lagrange multipliers (alpha is from [0,1]), if alpha = 0 no Largange multipliers are used:
                     % instead one replaces the calculated amplitude of the sensor plane wave field by the observed one
method = 0;          % method of the diffraction propagation and update of the sensor plane wave field
                     % if method = 0 - F-DDT with fitting of the sensor plane wave field to the observations                     
                     % if method = 1 - F-DDT with replacement of the sensor plane amplitude by the observation
                     % if method = 2 - ASD with fitting of the sensor plane wave field to the observations
                     % if method = 3 - ASD with replacement of the sensor plane amplitude by the observation
av = true;           % averaging for the F-DDT diffraction operators (true - with averaging, false - with no averaging)


% --------------- please do not change the code under this line -----------------------------------------
i = sqrt(-1); RMSE = []; tm=0;
if ~isreal(op) || ~isnumeric(op) || ceil(op)-floor(op)==1 
    disp('Error: Wrong setting of the operation, please choose op={0,1,2,3}.'); out = []; return 
end
[N,M,K]=size(o);                 % the size of the input set of the (noisy) intensity observations
if mod(N,2)~=0 || mod(M,2)~=0
   disp('Error: The size of the intensity observations must be even integers'); out = []; return 
end
% if K<2
%     disp('Error: It must be at least 2 intensity observations, otherwise phase retrieval doesnt make sense'); out = []; return 
% end
if (nargin<9 || isempty(gamma)),        gamma=dgamma*ones(K,1); end
if (nargin<14 || isempty(alpha)),       alpha=dalpha; end
if (nargin<11 || isempty(a0)), a0=da0; end
if any(~isreal(gamma)) || any(~isnumeric(gamma)) || any(gamma<=0)
    disp('Error: "gamma" must be a vector of K positive real values'); out = []; return 
end
if length(gamma)~=K
    disp(['Error: The length of "gamma" must be equal to K=' num2str(K)]);
    if length(gamma)==1
        gamma=gamma*ones(1,K); 
    end
end
if (nargin<10 || isempty(mufactor)),    mufactor=dmufactor; end
if ~isreal(mufactor) || ~isnumeric(mufactor) || mufactor<0
    disp('Error: "mufactor" must be a non-negative real value'); out = []; return 
end
if (nargin<12 || isempty(ITER)),  ITER=dITER; end
if q(1) < 1 || q(2) < 1, disp('Error: q must be not smaller than 1'); out = []; return ; end
Nya = ceil(q(1)*N/2)*2;  Nxa = ceil(q(2)*M/2)*2;     % the support of the computed Fourier images
cRMSE = false;                   % control parameter which control whether RMSE should be calculated
                                 % it the true u0 is given - RMSE is computed and the output variable 'out' is the calculated object 
                                 % otherwise - RMSE=[] is empty variable and 'out' is the calculated volumetric sensor plane wave field
if ~isempty(u0)                  % the input object is used to estimate the reconstruction accuracy in RMSE
    [Ny0,Nx0,tmp] = size(u0);  % the size of the input object
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
bourdery0 = Nya/2-Ny0/2+1:Nya/2+Ny0/2; bourderx0 = Nxa/2-Nx0/2+1:Nxa/2+Nx0/2; % center part of the object

% initialization (preallocation)
uzk_iter = zeros(Nya,Nxa,K);
% INITIAL GUESS
if mod(method,2)
    LAMBDA = 0; % the amplitude observation
    str = 'simple,';
else
    LAMBDA=zeros(Nya,Nxa,K);
    str = 'opt,';
end

% intial guess for the object: u0,0[k] = 1/2
u0est = zeros(Nya,Nxa);     a_=1/2;
u0est(bourdery0,bourderx0)=a_;        % phi = 0;

% figure(floor(op)+1),close(floor(op)+1),


%V=function_AL_DataSynthesis(u0,lambda,delta,z1,dz,K,[],1,1); 
[V,B,strs,tm] = TransformMatrix(method<2,z1,dz,lambda,deltaA,deltaB,Nya,Nxa,K,gamma,mufactor,av);

% 0. calculation of the wave fields at the sensor planes
for index = 1:K     
    uzk_iter(:,:,index)=ifft2((fft2(u0est)).*V(:,:,index));
    uzk_iter(coordy,coordx,index)=sqrt(o(:,:,index)).*exp(i*angle(uzk_iter(coordy,coordx,index)));
    
end

Izk_iter = uzk_iter;



% ----------------------------- AL ----------------------------------------
for jndex = 1:ITER
    tic;
    for index = 1:K
        if jndex>1
            uzk_iter(:,:,index)=ifft2((fft2(u0est)).*V(:,:,index));
            Izk_iter(:,:,index)=uzk_iter(:,:,index);
        end
        % minimization on ur
        if ~mod(method,2) % fitting to the observations / chosen with respect to the min value of the AL criterion            
            [cc]=CubicEquationAL(o(:,:,index),uzk_iter(coordy,coordx,index),LAMBDA(coordy,coordx,index),[gamma(index),mufactor],u0est(bourdery0,bourderx0));
            Izk_iter(coordy,coordx,index)=(uzk_iter(coordy,coordx,index)-LAMBDA(coordy,coordx,index))./cc;
        else % keep the amplitude, change the phase as in the conventional Gerchberg-Saxton-Fienup algorithm    
            Izk_iter(coordy,coordx,index) = uzk_iter(coordy,coordx,index)./abs(uzk_iter(coordy,coordx,index)).*sqrt(o(:,:,index));
        end
    end
    % minimization on u0
    switch floor(op)
        case 0  % complex-valued
            if jndex == 1, opp = ' complex'; end
            PP = SF(Izk_iter+LAMBDA,B,gamma);            u0est=ifft2(PP);
            out = u0est(bourdery0,bourderx0);
            
        case 1    % AM
            if jndex == 1, opp = ' AM'; end
            PP = SF(Izk_iter+LAMBDA,B,gamma);            u0est=ifft2(PP);
            out = real(u0est(bourdery0,bourderx0));            
            out = out.*(out>0)+(out<=0).*0.0001;% positive projection
            
        case 2     % PM with known amplitude a0, gradient descent AL
            if jndex == 1, opp = ' PM, known a0'; end
            ph = angle(u0est(bourdery0,bourderx0));              
            PP = SF(Izk_iter+LAMBDA,conj(V),gamma);
            
            for kndex=1:JTER
          
                XX = SF(a0*exp(i*ph),conj(V).*V,gamma);XX = ifft2(PP - XX);                
                YY = 2*imag(conj( conj(a0*exp(i*ph)) .*XX(coordy,coordx) ));                 
                ph = ph - beta*YY;                 
            end
            
            out = a0*exp(i*ph);
            
         case 3     % PM with unknown amplitude a_ , gradient descent AL
             if jndex == 1, opp = ' PM, unknown a0'; end
             ph = angle(u0est(bourdery0,bourderx0)); 
             PP = SF(Izk_iter+LAMBDA,conj(V),gamma); 
             for kndex=1:JTER
                 XX = SF(a_(end)*exp(i*ph),conj(V).*V,gamma);XX = ifft2(PP - XX);
                 YY = 2*imag(conj( conj(a_(end)*exp(i*ph)) .*XX(coordy,coordx) ));
                 ph = ph - beta*YY;
             end

             % recalculation of the constant amplitude
             S = 0; SS=0; ph = ph-mean(mean(ph)); u0est = exp(i*angle(u0est)); u0est(bourdery0,bourderx0) = exp(i*ph);
             PP = Izk_iter+LAMBDA;
             for index=1:K                 
                 %tmp_ = ifft2(RightHandSideTerm(u0est,1,conj(V(:,:,index))));
                 XX = ifft2(fft2(u0est).*V(:,:,index));                 
                 SS = SS + real( sum(sum( conj(XX(coordy,coordx)) .* PP(coordy,coordx,index) )) )./gamma(index);
                 S = S + sum(sum( conj(XX(coordy,coordx)) .* XX(coordy,coordx) ))./gamma(index);
             end
             a_=[a_ SS/(S+mufactor*N*M)];
             out = a_(end)*exp(i*ph);
        otherwise
            disp('Error: Wrong setting of the operation, please choose op={0,1,2,3}.'); out=[];t=0; return            
    end
     % zero padding (assume that there is nothing around the object in the object plane)
     u0est = zeros(Nya,Nxa); u0est(bourdery0,bourderx0)=out;
     tm = tm + toc;
     
     if ~mod(method,2)  % minimization on Lambda (we rewrite the values)           
         tic;
         for index=1:K
             LAMBDA(:,:,index)=LAMBDA(:,:,index)+alpha*(Izk_iter(:,:,index)-uzk_iter(:,:,index)); % center is different, bourders =0;
         end
         tm = tm + toc;
     end
% ------------------- visulaization of the result -------------------------           
     if cRMSE
         figure(4)
         v1 = abs(out) - abs(u0);  RMSE(1,jndex)=sqrt(mean(mean( ( v1 ).^2 )));
         x1 = angle(out) - angle(u0);   txmp = mean(x1(:)); x1 = x1-txmp;        RMSE(2,jndex)=sqrt(mean(mean( ( x1 ).^2 )));   
         subplot(2,3,1),    imshow(abs(out),[]), title('amplitude, AL'),xlabel(['RMSE=' num2str(RMSE(1,jndex))])
         subplot(2,3,2),    plot(RMSE(1,1:jndex)),  xlabel(['\gamma=' num2str(mean(gamma)) ',\alpha_r=' num2str(alpha) ',\mu=' num2str(mufactor)]) 
         if length(a_)>1
             subplot(2,3,2),hold on,plot(a_,'r--'),hold off
             legend('error','est. ampl.','Location','Best')
             axis([1,max(jndex,2),min(min(RMSE(1,1:jndex)),min(a_)),max(max(RMSE(1,1:jndex)),max(a_))]),title('convergence: amplitude'),ylabel('RMSE, a_0^{est}'),grid on
         else
             axis([1,max(jndex,2),min(min(RMSE(1,1:jndex)),0),max(max(RMSE(1,1:jndex)),eps)]),title('convergence: amplitude'),ylabel('RMSE'),grid on
         end
         subplot(2,3,4),imshow((angle(out)-txmp),[]), title('phase, AL'), xlabel(['RMSE=' num2str(RMSE(2,jndex))])
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


function [F,B,strs,t]=TransformMatrix(x,z1,dz,lambda,deltaA,deltaB,N,M,K,gamma,mufactor,av)
% calculation of the F-DDT or ASD transfer functions for the forward (F)
% and backward (B) diffraction propagation

F = zeros(N,M,K);  
if x,  QQ = zeros(N,M); strs = 'F-DDT '; 
else   QQ = 0;          strs = 'ASD ';
end
t = 0;    
for index = 1:K
    if x       [F(:,:,index),tmp] = TransferFunctionFDDT(z1+(index-1)*dz,lambda,deltaA,deltaB,N,M,av);         %function_TransferFunction_FDDT(z1+(index-1)*dz,lambda,deltaA,N,M,av);
    else       [F(:,:,index),tmp] = TransferFunctionASD(z1+(index-1)*dz,lambda,deltaA,deltaB,N,M); 
    end; t = t+tmp;
end
tic
for index = 1:K
    if x, QQ = QQ + abs(F(:,:,index)).^2 / gamma(index);
    else  QQ = QQ + 1/gamma(index);
    end
end
QQ = QQ + mufactor;
B = zeros(N,M,K);
for index = 1:K
    B(:,:,index) = conj(F(:,:,index)) ./ QQ;
end
t = t + toc;


function [Y] = SF(X,V,gamma)
[Nv,Mv,Kv] = size(V); [Nx,Mx,Kx] = size(X); 

if Nx == Nv && Mx == Mv && Kx == Kv
    Y = zeros(Nx,Mx);
    for index=1:Kx
        Y = Y + fft2(X(:,:,index)).*V(:,:,index)/gamma(index);
    end
end

if Nx < Nv && Mx < Mv && Kx == 1
    Y = zeros(Nv,Mv); XX = zeros(Nv,Mv);
    bourdery0 = Nv/2-Nx/2+1:Nv/2+Nx/2; bourderx0 = Mv/2-Mx/2+1:Mv/2+Mx/2;
    XX(bourdery0,bourderx0) = X;
    
    for index=1:Kv
        Y = Y + fft2(XX).*V(:,:,index)/gamma(index);
    end
end


function cc = CubicEquationAL(varargin)
% CubicEquation() serves to find the solution of the optimization problem min L(o_r,u_0,u_r,L_r) w.r.t. u_r,
% Here L = 1/2*abs(o_r-abs(u_r).^2).^2+1/gamma*abs(u_x-u_r).^2+2/gamma*real(conj(L_r).*(u_x-u_r)+mu*norm(u_0)^2,
% where     o_r         is the (noisy) intensity observation
%           u_0         is the object wave field
%           u_r         is the complex-valued wave field at a sensor plane
%           u_x         is the result of the diffraction propagation u_0 from the object to the sensor plane 
%                       (sensor plane estimate to be estimated) to be fitted to the given observation
%           L_r         is the Lagrange multiplier at the sensor plane
% cc = CubicEquation(o_r,u_r,L_r,[gamma,mufactor],u_0) results in cc = gamma*(abs(u_x_opt).^2-o_r)+1, 
% where u_x_opt is the optimazed sensor plane wave field to have the best fit to the observation
% u_x_opt is calculted using all these data: u_x, o_r, u_r, L_r and gamma


% -> input
Qr               = varargin{1};
Uzr              = varargin{2};
Lr               = varargin{3};
factor           = varargin{4};

u0est = varargin{5};  n_u0est = norm(abs(u0est)).^2;  clear u0est

% A.*X.^3 + B.*X.^2 +C.*X + D = 0 - general cubic equation, change the variables: X -> Y=X+B./3/A  (Tschirnhaus transformation), then
% P = (3*A.*C-B.B) ./ (3*A.A), Q = (2*B.^3 - 9*A*B*C + 27*A*A*D) ./ (27*A.^3) and we are forced to conclude that Y.^3 + P.*Y + Q = 0, 
% and solve it with the method due to Scipione del Ferro and Tartaglia, published by Gerolamo Cardano in 1545.

P = 1/factor(1)-Qr;  
Q = -1/factor(1).*abs(Uzr-Lr);  % # 1 cc>0

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


JJtmp = 10000*ones(NN,MM);JJtmp(Y1) = JR(Qr(Y1),Uzr(Y1),rX1(Y1),Lr(Y1),factor,n_u0est,sign(X1(Y1)));
JJ(:,:,1)=JJtmp;JJtmp = 10000*ones(NN,MM);JJtmp(Y2) = JR(Qr(Y2),Uzr(Y2),rX2(Y2),Lr(Y2),factor,n_u0est,sign(X2(Y2)));
JJ(:,:,2)=JJtmp;JJtmp = 10000*ones(NN,MM);JJtmp(Y3) = JR(Qr(Y3),Uzr(Y3),rX3(Y3),Lr(Y3),factor,n_u0est,sign(X3(Y3)));
JJ(:,:,3)=JJtmp;minJJ3 = min(JJ,[],3);rrY = double( minJJ3 == JJ(:,:,1)  ).*rX1+double(  minJJ3 == JJ(:,:,2)  ).*rX2+double(  minJJ3 == JJ(:,:,3)   ).*rX3;
cc = factor(1)*(abs(rrY).^2-Qr)+1;
 
   
function [out] = JR(Qr,Uzr,Vr,Lr,factor,n_u0est,sgn)
% this function just calculates the criterion to be minimized w.r.t. the current values of components
% -> Qr - observations, 
% -> Uzr - estimate of the sensor wave field denoted in the AL algorithm as ur, t+1/2 
% -> Vr - estimation (preliminary) of the sensor wave field denoted in the AL algorithm as ur, t+1 
% -> Lr - Lagrange multipliers
% -> n_u0est - a norm of the estimate of the object wave field
% -> factor = [gamma, mufactor] - only two parameters !!! why so see below
% -> sgn - a matrix of sign of the denuminator
% Vr must be real and positive
cc1 = factor(1)*sgn.*(abs(Vr).^2-Qr)+1;
Izk_iter=(Uzr-Lr)./cc1; 
out = 1/2*abs(Qr-abs(Izk_iter).^2).^2+1/factor(1)*abs(Izk_iter-Uzr).^2+2/factor(1)*real(conj(Lr).*(Izk_iter-Uzr))+factor(2)*n_u0est;

