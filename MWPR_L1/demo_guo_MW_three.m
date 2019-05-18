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

addpath('./function');

u0 = im2double(imread('foreman.bmp'));     % read test-image
if size(u0,3)>1,    u0=rgb2gray(u0);end
u0 = u0+0.1;    u0 = u0/max(u0(:));         % make all pixels positive in [0.(09),1]
lambda = [488e-9 532e-9 632.8e-9];                            % wavelength [m]
delta = 1.34e-6;                             % square pixel size with 100% fill factor [m]
ITER = 100;                                 % the number of iterations
K = 3;                                      % number of measurements
sigma = 0.05*ones(K,1);                        % noise level of the intnesity observations

[N,M]=size(u0);                             % object size
Nz = N; Mz = M;                             % sesnor size (be default, the same as the object size)

z1 = 2* delta*delta*min([Nz Mz])/lambda(3);    % distance to the first sensor plane


% -------------------- synthesis of the inensity observations ------------- 
uz = zeros(Nz,Mz,K); o = zeros(Nz,Mz,K);noise = randn(Nz,Mz,K);
for index = 1:K
    if 0 % calculation by F-DDT
        uz(:,:,index) = FDDT(u0,0,z1,lambda(index),delta,delta,[Nz/N,Mz/M]);
        
%         xx = zeros(Nz,Mz); coordy = Nz/2-N/2+1:Nz/2+N/2; coordx = Mz/2-M/2+1:Mz/2+M/2; xx(coordy,coordx)=u0;
%         S = TransferFunctionFDDT(z1,lambda(index),delta,delta,Nz,Mz,true);
%         uz(:,:,index) = ifft2(fft2(xx).*S);
        
    else % calculation by ASD
%         xx = zeros(Nz,Mz); coordy = Nz/2-N/2+1:Nz/2+N/2; coordx = Mz/2-M/2+1:Mz/2+M/2; xx(coordy,coordx)=u0;
        S = TransferFunctionASD(z1,lambda(index),delta,delta,Nz,Mz);
        uz(:,:,index) = ifft2(fft2(u0).*S);
    end
    o(:,:,index)= abs(uz(:,:,index)).^2 + sigma(index)*noise(:,:,index);
    o(:,:,index)=o(:,:,index).*(o(:,:,index)>=0) + 0.0001*(o(:,:,index)<0); % positive projection to avoid negative or zero intensity
end

% o=padarray(o,[N/2,M/2]);
% u0=padarray(u0,[N/2,M/2]);

[u0_average,tAPR_average,eAPR_average] = MWPR_average(o,u0,z1,lambda,delta,delta,ITER);

[u0APR_serial,tAPR_serial,eAPR_serial] = MWPR_serial(o,u0,z1,lambda,delta,delta,ITER);

%[u0APR_pc,tAPR_pc,eAPR_pc] = MWPR_phasecontrast(o,u0,z1,lambda,delta,delta,ITER);

opADMM = 0;
gamma = 0.1;
mu=5e-3;
a0=1;
alpha=1;
[u0ADMM,tADMM,eADMM] = ADMM_MWPR(o,u0,opADMM,z1,lambda,delta,delta,gamma,mu,a0,ITER,alpha);

[u0ADMM,tADMM,eADMM] = ADMMTV_MWPR(o,u0,opADMM,z1,lambda,delta,delta,gamma,mu,a0,ITER,alpha);