function output = CG(o,u0,op,z1,dz,lambda,deltaA,deltaB,gamma,mu,a0,ITER,alpha)
dgamma = 0.1;
dmu = 5e-3;
da0 = 1;
dITER = 500;
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
if method < 2
    Nya = Nya+Ny0;Nxa = Nxa+Nx0;    
end
coordy = Nya/2-N/2+1:Nya/2+N/2; coordx = Nxa/2-M/2+1:Nxa/2+M/2; % center part of the sensor plane wave field
bourdery = Nya/2-Ny0/2+1:Nya/2+Ny0/2; bourderx = Nxa/2-Nx0/2+1:Nxa/2+Nx0/2; % center part of the object

gamma0=0.9;
gamma1=0.05;
gamma2=0.05;
gamma3=0.01;
gammar = gamma0/K*ones(K,1);
gamma_scaler = 2;

[V,B,strs,tm]=TransformMatrix(0,z1,dz,lambda,deltaA,deltaB,Nya,Nxa,K,gammar,mu,av);


u0_iter = zeros(Nya, Nxa); a_=1/2;
u0_iter=a_+u0_iter; 


for jndex = 1:ITER
    PP = zeros(Nya, Nxa);
    for index = 1:K
        PP1 = ifft2(fft2(u0_iter).*V(:,:,index));
        PP1(bourdery,bourderx) = PP1(bourdery,bourderx) - o(:,:,index);
        PP = PP + ifft2(fft2(PP1).*B(:,:,index));
    end
    PP1 = gamma0*PP + abs(gamma1*TV(TV(u0_iter,1),1)+gamma2*TV(TV(u0_iter,2),2));
    u0_iter = u0_iter-gamma3*PP1;
    out = u0_iter(bourdery,bourderx);
    if cRMSE
        figure(5)
         v3 = abs(out) - abs(u0);  RMSE(1,jndex)=sqrt(mean(mean( ( v3 ).^2 )));
         x1 = angle(out) - angle(u0);   txmp = mean(x1(:)); x1 = x1-txmp;        RMSE(2,jndex)=sqrt(mean(mean( ( x1 ).^2 )));   
         subplot(2,3,1),    imshow(abs(out),[]), title('amplitude, ADMM'),xlabel(['RMSE=' num2str(RMSE(1,jndex))])
         subplot(2,3,2),    plot(RMSE(1,1:jndex)),  xlabel(['\gamma=' num2str(mean(gamma)) ',\alpha_r=' num2str(alpha) ',\mu=' num2str(mu)]) 
         drawnow 
    end
    gamma1 = gamma1;
    gamma2 = gamma2;
    gamma3 = gamma3;
end
output=u0_iter;