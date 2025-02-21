% This demo is made for the object wave field reconstruction and a qualitative comparison of reconstruction
% imaging of the TUT logo test-image (amplitude-only object,  [1, Fig. 2]). The object estimates can be
% obtained by the successive SBMIR or parallel AL phase-retrieval algorithms 
% [1] A. Migukin, V. Katkovnik, and J. Astola,  J. Opt. Soc. Am. 28, 993-1002 (2011). 
% ------------------- default setup parameters ----------------------------
% test-image (NxM):                           TUT logo (256x256), binary
% object type:                                amplitude-only (AM)
% wavelength :                                lambda = 532 nm
% pixel size (delta):                         6.7 mkm (100% fill factor)
% number of iterations:                       ITER = 100
% number of measurements:                     K = 5 
% noise level at all senor planes:            sigma = 0
% distance to the first sensor plane:         z1 = 64.8 mm
% distance between sensor planes):            dz = 2 mm
% penalty coefficient for all sensor planes:  gamma = 10
% regularization parameter:                   mufactor = 0.005 (recommended value)

clc; clear; close all

u0 = im2double(imread('monarch.bmp'));     % read test-image
if size(u0,3)>1,    u0=rgb2gray(u0);end
u0 = u0+0.1;    u0 = u0/max(u0(:));         % make all pixels positive in [0.(09),1]
% u0=padarray(u0,[50,50]);
lambda = 532e-9;                            % wavelength [m]
delta = 3.1e-6;                             % square pixel size with 100% fill factor [m]
ITER = 100;                                 % the number of iterations
K = 2;                                      % number of measurements
sigma = 0.1*ones(K,1);                        % noise level of the intnesity observations

[N,M]=size(u0);                             % object size
Nz = N; Mz = M;                             % sesnor size (be default, the same as the object size)

z1 = 2* delta*delta*min([Nz Mz])/lambda;    % distance to the first sensor plane
% z1=2e-3;
dz = 1e-3;                                  % distance between measurement planes [m]

gamma = 10*ones(K,1);                       % penalty coefficient of the AL base algorithms

% -------------------- synthesis of the inensity observations ------------- 
uz = zeros(Nz,Mz,K); o = zeros(Nz,Mz,K);noise = randn(Nz,Mz,K);
for index = 1:K
    if 1 % calculation by F-DDT
        uz(:,:,index) = FDDT(u0,0,z1+(index-1)*dz,lambda,delta,delta,[Nz/N,Mz/M]);
    else % calculation by ASD
        xx = zeros(Nz,Mz); coordy = Nz/2-N/2+1:Nz/2+N/2; coordx = Mz/2-M/2+1:Mz/2+M/2; xx(coordy,coordx)=u0;
        S = TransferFunctionASD(z1+(index-1)*dz,lambda,delta,delta,Nz,Mz);
        uz(:,:,index) = ifft2(fft2(xx).*S);
    end
    o(:,:,index)= abs(uz(:,:,index)).^2 + sigma(index)*noise(:,:,index);
    o(:,:,index)=o(:,:,index).*(o(:,:,index)>=0) + 0.0001*(o(:,:,index)<0); % positive projection to avoid negative or zero intensity
end


% opSBMIR = 2;                                % forward-forward 
% [u0SBMIR,tSBMIR,eSBMIR] = SBMIR(o,u0,opSBMIR,z1,dz,lambda,delta,delta,ITER);
% 
%  
% [u0APR,tAPR,eAPR] = APR3(o,u0,z1,dz,lambda,delta,delta,ITER);
% 
% 
% a=0.7;b=0.5;
% [u0APRDF,tAPRDF,eAPRDF]=APRDF3(o,u0,z1,dz,lambda,delta,delta,ITER,a,b);


% mufactor = 0.02;                           % regularization parameter
% opAL = 0;                                   % if op = 0 - complex-valued object, set opAL = 1 to use a priori information that the object is amplitude-only
% [u0AL,tAL,eAL] = AL(o,u0,opAL,z1,dz,lambda,delta,delta,gamma,mufactor,[],ITER);


opADMM = 0;                                   % if op = 0 - complex-valued object, set opAL = 1 to use a priori information that the object is amplitude-only
% num=50;
% o=padarray(o,[num,num]);
% u0=padarray(u0,[num,num]);
[u0ADMM,tADMM,eADMM] = ADMM(o,u0,opADMM,z1,dz,lambda,delta,delta,[]);