close all;
clear all;clc;
addpath('./Functions');
addpath('./FISTA-master');
addpath('./FISTA-master');
addpath('./FISTA-master/utils');
addpath('./FISTA-master/proj');
addpath('./BM3D');

str1='lena256';
str2='.bmp';
str=strcat(str1,str2);
o = imread(str);
%o = rgb2gray(o);
o = imresize(o,[500,500]);
o = 1-im2double(o);
o_compare = 1-o;
%o = imnoise(o,'gaussian',0,0.001);

% inout signal f 
f=o;
%figure;imshow(1-abs(f),[],'border','tight');

nx=size(f,1);  % data size
ny=size(f,2);
nz=1;
sigma1 = 20; 
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

%加噪声，生成孪生像
[g,E,Phase]=noise(f,nx,ny,nz,lambda,z,deltaX,deltaY,sigma1);
str3=strcat(strcat(str1,'1'),str2);
imwrite(1-abs(g),strcat('./result/',str3));
g1=MyC2V(g(:));
tf1=MyAdjointOperatorPropagation(g1,E,nx,ny,nz,Phase);
tf1 = abs(reshape(MyV2C(tf1),nx,ny,nz));


%去除孪生像和噪声

g = im2uint8(g);
[ssss,g]=BM3D(1,g,sigma1,'np',1);
g = im2double(mat2gray(g));
g=MyC2V(g(:));

transf=MyAdjointOperatorPropagation(g,E,nx,ny,nz,Phase);

%figure;imshow(tf1,[]);title('观测到的图像')
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
re = 1-im2double(re);
re = im2uint8(re);
[ssss,re]=BM3D(1,re,sigma1,'np',1);
re = im2double(mat2gray(re));

tf1 = 1-tf1;
Ssim=ssim(tf1,o_compare);
[peaksnr,snr] =psnr(tf1,o_compare);
fprintf('恢复前的Ssim：%d\n',Ssim);
fprintf('恢复前的psnr：%d\n',peaksnr);

Ssim=ssim(re,o_compare);
[peaksnr,snr] =psnr(re,o_compare);
fprintf('恢复后的Ssim：%d\n',Ssim);
fprintf('恢复后的psnr：%d\n',peaksnr);

figure;imshow(tf1,[]);title('观测到的图像')
figure;imshow(re,[]);title('恢复后的图像')

str3=strcat(strcat(str1,'2'),str2);
imwrite(tf1,strcat('./result/',str3));

str3=strcat(strcat(str1,'3'),str2);
imwrite(re,strcat('./result/',str3));