close all;
clear all;clc;
addpath('./Functions');
addpath('./FISTA-master');
addpath('./FISTA-master');
addpath('./FISTA-master/utils');
addpath('./FISTA-master/proj');
addpath('./BM');

str1='cell';
str2='.jpg';
str=strcat(str1,str2);
o = imread(str);
o = rgb2gray(o);
o = imresize(o,[500,500]);
o = 1-im2double(o);
o_compare = 1-o;
%o = imnoise(o,'gaussian',0,0.001);

% inout signal f 
f=o;
%figure;imshow(1-abs(f),[],'border','tight');


%% Paprameters (1)
nx=size(f,1);  % data size
ny=size(f,2);
nz=1;

%lambda=0.532;  % wavelength (um)
lambda=0.469;
%lambda=0.5;
k = 2*pi/lambda;
detector_size=4;  % pixel pitch (um)
sensor_size=nx*detector_size;  % detector size (um)
z=12000;  % distance from detector to first reconstructed plane (um)
deltaX=detector_size;
deltaY=detector_size;
Nx=nx;
Ny=ny*nz*2;
Nz=1;

%% Propagation kernel (2)
Phase=MyMakingPhase(nx,ny,z,lambda,deltaX,deltaY);
%figure;imagesc(plotdatacube(angle(Phase)));title('Phase of kernel');axis image;drawnow;
%axis off; colormap(hot); colorbar;
E0=ones(nx,ny);  % illumination light
E=MyFieldsPropagation(E0,nx,ny,nz,Phase);  % propagation of illumination light

%% Field measurement and backpropagation (3)
cEs=zeros(nx,ny,nz);
Es=f.*E;
for i=1:nz
    cEs(:,:,i)=fftshift(fft2(Es(:,:,i)));
end
cEsp=sum(cEs.*Phase,3);
S=(ifft2(ifftshift(cEsp)));


f1 = ones(nx,ny);
Es1=f1.*E;
for i=1:nz
    cEs1(:,:,i)=fftshift(fft2(Es1(:,:,i)));
end
cEsp1=sum(cEs1.*Phase,3);
S1=(ifft2(ifftshift(cEsp1)));

% squared field
s=(S+1).*conj(S+1);
s1=(S1+1).*conj(S1+1);
%  diffracted field
g = s./s1; % normalized 
g = im2double(g);
%g=imnoise(g,'salt & pepper',0.001);

%figure;imshow(abs(g),[]);title('Diffracted field')

g=MyC2V(g(:));
transf=MyAdjointOperatorPropagation(g,E,nx,ny,nz,Phase);


%transf=reshape(MyV2C(transf),nx,ny,nz);
%figure;imshow(plotdatacube(abs(transf)),[],'border','tight')

%% Propagation operator (4)
A = @(f_twist) MyForwardOperatorPropagation(f_twist,E,nx,ny,nz,Phase);  % forward propagation operator
AT = @(g) MyAdjointOperatorPropagation(g,E,nx,ny,nz,Phase);  % backward propagation operator

tau = 0.005; 
piter = 4;
tolA = 1e-6;
iterations = 200;

Phi = @(f) MyTVphi(f,Nx,Ny,Nz);
Psi = @(f,th) MyTVpsi(f,th,0.05,piter,Nx,Ny,Nz);

opts.pos = true;
opts.lambda = lambda;
opts.check_grad = 0;
opts.verbose = true;
f_reconstruct = fista_row_sparsity1(g, A, AT,Phi,Psi,transf, opts,Nx,Ny,Nz);

f_reconstruct=reshape(MyV2C(f_reconstruct),nx,ny,nz);
re=(abs(f_reconstruct));
re = im2double(re);
%figure,imshow(re,[])
re1 = 1-im2double(re);
%re1 = wiener2(re1,[5,5]);
transf1=reshape(MyV2C(transf),nx,ny,nz);
transf1=1-im2double(abs(transf1));
Ssim=ssim(re1,o_compare);
[peaksnr,snr] =psnr(re1,o_compare);
fprintf('%d\n',Ssim);
fprintf('%d\n',peaksnr);

%Ssim=ssim(transf1,o_compare);
%[peaksnr,snr] =psnr(transf1,o_compare);
%fprintf('%d\n',Ssim);
%fprintf('%d\n',peaksnr);
%figure;imshow(transf1,[])
figure;imshow(re1,[])
str3=strcat(strcat(str1,'3'),str2);
imwrite(re1,strcat('./result/',str3));

