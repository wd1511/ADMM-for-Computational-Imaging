function [out,t,RMSE] = MWPR_serial(o,u0,z1,lambda,deltaA,deltaB,ITER)

% default parameters (if not specified)
dITER = 100;         % number of iteration of the phase-retrieval algorithm
q = [1,1];           % the redundancy of the computantional size w.r.t. the y- and x-direction, qx>=1 and qy>=1 
av = true;           % (by default ASD is used of the same size as the size of the object/sensor plane wave field dq=[1,1],
                     % for instance, for realization of ASD via double size zero-padding , please set dq = [2,2].)

                     
% --------------- please do not change the code under this line -----------------------------------------
i = sqrt(-1); RMSE = []; t=0;
if (nargin<9 || isempty(ITER)),  ITER=dITER; end

[N,M,K]=size(o);                 % the size of the input set of the (noisy) intensity observations
if mod(N,2)~=0 || mod(M,2)~=0
   disp('Error: The size of the intensity observations must be even integers'); out = []; return 
end
cRMSE = false;                   % control parameter which control whether RMSE should be calculated
                                 % it the true u0 is given - RMSE is computed and the output variable 'out' is the calculated object 
                                 % otherwise - RMSE=[] is empty variable and 'out' is the calculated volumetric sensor plane wave field
if q(1) < 1 || q(2) < 1, disp('Error: q must be not smaller than 1'); out = []; return; end
Nya = ceil(q(1)*N/2)*2;  Nxa = ceil(q(2)*M/2)*2;     % the support of the computed Fourier images
% if K<2
%     disp('Error: It must be at least 2 intensity observations, otherwise phase retrieval doesnt make sense'); out = []; return 
% end
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

% Nya = Nya+Ny0;Nxa = Nxa+Nx0;

o = sqrt(o); % the amplitude observation
coordy = Nya/2-N/2+1:Nya/2+N/2; coordx = Nxa/2-M/2+1:Nxa/2+M/2; % center part of the sensor plane wave field
bourdery0 = Nya/2-Ny0/2+1:Nya/2+Ny0/2; bourderx0 = Nxa/2-Nx0/2+1:Nxa/2+Nx0/2; % center part of the object

% calculation of the transfer functions for the forward and backward ASD diffraction propagation to the distance dz
if nargout > 1,    tic;end

% Sf = TransferFunctionASD( dz,lambda,deltaA,deltaB,Nya,Nxa);
for index=1:K
    Sf(:,:,index) = TransferFunctionASD(z1,lambda(index),deltaA,deltaB,Nya,Nxa);
    Sb(:,:,index) = TransferFunctionASD(-(z1),lambda(index),deltaA,deltaB,Nya,Nxa);
    
%     Sf(:,:,index) = TransferFunctionFDDT(z1,lambda(index),deltaA,deltaB,Nya,Nxa,av);
%     Sb(:,:,index) = TransferFunctionFDDT(-z1,lambda(index),deltaA,deltaB,Nya,Nxa,av);
end

if nargout > 1,    t=t+toc; end

u0est=ones(Nya,Nxa);
for jndex = 1:ITER
     if nargout > 1,    tic; end
     for index=1:K
         Ar_u0(:,:)=ifft2((fft2(u0est)).*Sf(:,:,index));
         Ar_u0(coordy,coordx)=o(:,:,index).*exp(i*angle(Ar_u0(coordy,coordx)));
         u0est=ifft2((fft2(Ar_u0)).*Sb(:,:,index));
     end
%      u0est=mean(o_est,3);
    if nargout > 1,   t=t+toc; end
     if cRMSE
        if nargout > 1 && jndex==ITER,  tic;     end 
        out = u0est(bourdery0,bourderx0); % the reconstructed object to be estimated
        if nargout > 1 && jndex==ITER,  t=t+toc; end
        figure(1)
        v2=abs(u0) - abs(out); RMSE(1,jndex)=sqrt(mean(mean( v2.^2 )));
        x2=angle(u0) - angle(out);tmp=mean(x2(:));x2 = wrap(x2-tmp); RMSE(2,jndex)=sqrt(mean(mean( x2.^2 )));
        
        subplot(2,3,1),imshow(abs(out),[]), title('amplitude, serial'),xlabel(['RMSE=' num2str(RMSE(1,jndex))])
        subplot(2,3,2),plot(RMSE(1,1:jndex)), title('convergence: amplitude'),ylabel('RMSE'),xlabel('# of iterations');grid on 
        
        subplot(2,3,4),imshow(angle(out)+tmp,[]), title('phase, serial'), xlabel(['RMSE=' num2str(RMSE(2,jndex))])
        subplot(2,3,5),plot(RMSE(2,1:jndex)), title('convergence: phase'),ylabel('RMSE'), xlabel(['K=' num2str(K) ',iter=' num2str(jndex) '/' num2str(ITER)]),grid on
        
        subplot(2,3,3),    plot(1:Nx0,abs(out(Ny0/2+1,:)),'b',1:Nx0,abs(u0(Ny0/2+1,:)),'r'),title(['amplitude: cross-section along the ' num2str(Ny0/2+1) '-th line'])
        axis([1,Nx0,0,max( max(abs(out(Ny0/2+1,:))), max(abs(u0(Ny0/2+1,:))) )])
        subplot(2,3,6),    plot(1:Nx0,angle(out(Ny0/2+1,:))+tmp,'b',1:Nx0,angle(u0(Ny0/2+1,:)),'r'),title(['phase: cross-section along the ' num2str(Ny0/2+1) '-th line'])
        axis([1,Nx0, min( min(angle(out(Ny0/2+1,:))+tmp), min(angle(u0(Ny0/2+1,:))) )    ,max( max(angle(out(Ny0/2+1,:))+tmp), max(angle(u0(Ny0/2+1,:))) )])
        
        drawnow
    end
end
   
fprintf('elapsed time = %8.4f seconds\n',t);

function y=wrap(x)
% wrapping of the phase
y = atan2(sin(x),cos(x));    y(y==pi)=-pi; %  alternative > y = mod(x+pi,2*pi)-pi;

% 2010 - 2013 (c) TUT v. 2.0