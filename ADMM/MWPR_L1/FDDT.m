function [out,t,Vo]=FDDT(in,op,z,lambda,deltaA,deltaB,q,av,reg,V)
% FDDT() realizes the pixel-to-pixel discrete diffraction propagation in the frequency domain described in [1]. 
% This function calculates the sensor plane wave field using the input (complex-valued) object or
% the object reconstruction from the given/computed complex-valued observation by the Tikhonov regularized inverse
% or via the recursive regularized inverse (iterative procedure computed the boundary problem). The (rectangular) 
% pixels at the object and senor planes are assumed to be of the same size deltaA x deltaB. 
% By default FDDT uses double size zero-padded object in toder to satisfy the support of the discrete convolution.
%
% Notation:
% [out,t,Vo] = FDDT(in,op,z,lambda,deltaA,deltaB,q,av,reg,V)
% ------------------------------- input -----------------------------------
% in                    the input (object or sensor plane) wave field distribution depending on
%                       the requied operation op of the forward/backward diffraction propagation
%                       The size of the "in" must be even integers!                       
% op                    this parameter determines the perfomed operation:
%                       If op = 0, FDDT() realizes the forward wave field propagation and results in the sensor 
%                       plane wave field. In this case the input parameter must be the object wave field. 
%                       If op = 1, FDDT() realizes the Tikhonov regularize inverse and results in the estimate
%                       of the object wave field. Here the input parameter must be the sensor plane wave field. 
%                       If op = 2, FDDT() realizes the recursive regularized inverse typically by 10 iterations. 
%                       The input parameter must be the sensor plane wave field.
% z                     propagation distance between the object and sensor planes [m]
% lambda                wavelength [m]
% deltaA                pitch size w.r.t. rows of the discrete object/sensor plane wave field [m]
% deltaB                pitch size w.r.t. columns of the discrete object/sensor plane wave field [m] 
% q = [qy,qx]           redundancy of the sensor (the ratio of the sesor plane wave field distribution 
%                       w.r.t. the object one) in the y- and x-direction, respectively (by default for the same 
%                       size of the object and sensor plane wave fields qy = 1, qx = 1, must be qy >= 1, qx >= 1).
% av                    this parameters determines which transfer function to use in the F-DDT model (optional): 
%                       If av = 1 (by default) - the transfer function with averaging.
%                       If av = 0 - the transfer function with no averaging.
% reg                   If op = 1 or op = 2 and the regularized inverse is about to be performed, then
%                       reg must be the regularization parameter (positive scalar, by default reg=...)                       
% V                     the F-DDT transfer function used for accelaration of calculation (optional)
% ------------------------------- output ----------------------------------
% out                   the computed wave field: if op = 0 - the sensor plane wave field 
%                                                if op = 1 or op = 2 - the object estimate 
% t                     the computational time in seconds (optional)
% Vo                    the F-DDT transfer function for the chosen parameters
%
% Example #1:
% Calculation of the sensor plane wave field uz of the size M1 x M2 using the 
% amplitude object u0 of the size N1 x N2. Here we use rectangular pixels at the
% object and sensor planes of the same size deltaA x deltaB.
%
% u0 = im2double(imread('Lena256.png'));                 % read the test-image
% [N1,N2,N3] = size(u0);                                 % image size
% lambda = 6328e-10;                                     % wavelength (red laser)
% deltaA = 6e-6; deltaB = 8e-6;                          % the size of the rectangular pixels at the object 
% M1 =  N1; M2 = N2;                                     % and sensor planes are assumed to be even
% z = deltaA*deltaA*M1/lambda;                           % propagation distance
% uz = FDDT(u0,0,z,lambda,deltaA,deltaB);
%
% Example #2:
% Simulation of the F-DDT pixel-to-pixel diffraction propagation of rectangular pixel and images of the same sizes at the object and sensor planes.  
% The discrete amplitude-only object u00 is of the size 512x256 with the pixels of the size 5x8 (microns).
% The sesnor plane discrete wave field uz (the result of the forward propagation) is calculated of the size 616x334.
% The reconstruction of the object (backward diffraction propagation) is performed by the (recursive) regularized inverse.
% u0 = im2double(imread('Baboon512.png'));
% u00 = u0(:,129:384);                                   % rectangular amplitude-only object of the size 512x256
% [N1,N2,N3] = size(u00);
% lambda = 532e-9;                                       % wavelength (green laser)
% deltaB = 8e-6; deltaA = 5e-6;                          % rectangular pixel size for the object
% q = [1.2,1.3];                                         % redundancy of the senor to the size 616x334
% M1 =  q(1)*N1; M2 = q(2)*N2;                           % the sensor plane wave field is taken of the same size as the object size
% z = deltaA*deltaA*M1/lambda;                           % propagation distance
% uz = FDDT(u00,0,z,lambda,deltaA,deltaB,q); % the estimate of the sensor plane wave field
% u0F = FDDT(uz,1,z,lambda,deltaA,deltaB,q,[],1e-3);     % object reconstruction by the reqularized inverse 
%                                                        % with the regularization parameter equal to 1^(-3)
% u0FRRI = FDDT(uz,2,z,lambda,deltaA,deltaB,q,[],16e-3); % object reconstruction by the recursive (10 iterations) regularized inverse 
%                                                        % with the regularization parameter equal to 16^(-3)
%
% Reference:            [1] V. Katkovnik, J. Astola, and K. Egiazarian, Appl. Opt. 47, 3481-3493 (2008).  
% Authors:              Artem Migukin, V. Katkovnik
% Software version:     2.0, January 2013.
% Required functions:   TransferFunctionFDDT()

ITER = 10;                         % number of iteration for the recursive regularized inverse of F-DDT

%% --------------- please do not change the code under this line -----------------------------------------
if (nargin<7 || isempty(q)),    q=[1,1]; end
if length(q)~=2 || ~isreal(q(1)) || ~isnumeric(q(1)) 
    disp('Error: Wrong setting of the parameter q=[qy,qx], Use the default q=[1,1]!'); q = [1,1]; 
end
if q(1) < 1 || q(2) < 1, disp('Error: qy or qx must be not smaller than 1'); q = [1,1]; end
if (nargin<8 || isempty(av)),  av=true; end          % by default the averaging for transform matrix is switched on 

if op == 0                                                                  % forward propagation
    [Ny0_,Nx0_,tmp] = size(in);                                             % in is the object
    if mod(Ny0_,2)~=0 || mod(Nx0_,2)~=0
        disp('Error: the size of the input signal must be an even integer!');  out=[]; t=0; Vo=[]; return
    end
    Ny0 = ceil(Ny0_/2)*2;Nx0 = ceil(Nx0_/2)*2;
    tmp = zeros(Ny0,Nx0); tmp(Ny0/2-Ny0_/2+1:Ny0/2+Ny0_/2,Nx0/2-Nx0_/2+1:Nx0/2+Nx0_/2)=in; in = tmp;
    Nyz = ceil(Ny0*q(1)/2)*2; Nxz = ceil(Nx0*q(2)/2)*2;                     % make the size of the sensor plane wave field even
    Nxa = Nx0 + Nxz;  Nya = Ny0 + Nyz;                                      % the support of the discrete convolution +1 (to use even values)
else
    [Nyz_,Nxz_,tmp]=size(in);                                               % in is the senser plane wave field
    Nyz = ceil(Nyz_/2)*2;Nxz = ceil(Nxz_/2)*2;               
    tmp = zeros(Nyz,Nxz); tmp(Nyz/2-Nyz_/2+1:Nyz/2+Nyz_/2,Nxz/2-Nxz_/2+1:Nxz/2+Nxz_/2)=in; in = tmp;
    Ny0 = floor(Nyz/q(1)/2)*2; Nx0 = floor(Nxz/q(2)/2)*2; % find even size of the object for the given q=[qy,qx]
    if (nargin<9 || isempty(reg)),  reg= 1e-7; end
    if length(reg)~=1 || reg(1)<=0 || ~isreal(reg(1)) || ~isnumeric(reg(1))
        disp('Error: Wrong setting of the regularization parameter. reg should be a positive real-valued scalar'); 
        out=[]; t=0; Vo=[]; return
    end
    
    if ~(nargin<10 || isempty(V))
        [Nya,Nxa,tmp]=size(V); 
    else
        Nxa = Nx0 + Nxz;  Nya = Ny0 + Nyz;                                  % the support of the discrete convolution +1 (to use even values)
    end
end


bourdery0 = Nya/2-Ny0/2+1:Nya/2+Ny0/2; bourderx0 = Nxa/2-Nx0/2+1:Nxa/2+Nx0/2; % center part of the object
bourderyz = Nya/2-Nyz/2+1:Nya/2+Nyz/2; bourderxz = Nxa/2-Nxz/2+1:Nxa/2+Nxz/2; % center part of the sensor plane wave field

if nargin<10 || isempty(V)
    [Vo,t] = TransferFunctionFDDT(z,lambda,deltaA,deltaB,Nya,Nxa,av);
else  tic;  Vo = V; t=toc; clear V
end

if op == 0
    if nargout > 1,    tic;end                                              % start a stopwatch timer 
    u0 = zeros(Nya,Nxa);  u0(bourdery0,bourderx0) = in;                     % 'in' here is the object plane wave field
    uz = ifft2(fft2(u0).*Vo); out = uz(bourderyz,bourderxz);
    if nargout > 1,     t=t+toc; % stop the stopwatch timer - obtain the time of calculation
    end
elseif op == 1
    if nargout > 1,    tic;end                                              % start a stopwatch timer 
    rVo = conj(Vo)./(abs(Vo).^2 + reg);                                     % inverse FDDT transfer function 
    uz = zeros(Nya,Nxa);  uz(bourderyz,bourderxz) = in;                     % 'in' is the sensor plane wave field 
    u0est = ifft2(rVo.*fft2(uz)); out = u0est(bourdery0,bourderx0);
    if nargout > 1,     t=t+toc; % stop the stopwatch timer - obtain the time of calculation
    end
elseif op == 2
    fprintf('Please wait...');
    if nargout > 1,    tic;end                                              % start a stopwatch timer 
    rVo = conj(Vo)./(abs(Vo).^2 + reg);                                     % inverse FDDT transfer function 
    uz = zeros(Nya,Nxa);  uz(bourderyz,bourderxz) = in;                     % 'in' is the sensor plane wave field 
    if ITER<1
        disp('It must be at least one iteration of the reconstruction procedure!'); 
    end
    
    tmp = ifft2(rVo.*fft2(uz)); out = tmp(bourdery0,bourderx0); t = t+toc;
    if ITER>1        
        for index=1:ITER-1
            tic
            tmp = zeros(Nya,Nxa); tmp(bourdery0,bourderx0) = out;
            uz = ifft2(fft2(tmp).*Vo); uz(bourderyz,bourderxz) = in; 
            tmp = ifft2(rVo.*fft2(uz)); out = tmp(bourdery0,bourderx0);            
            t=t+toc; fprintf(' . %d ',index);
        end
        fprintf('\n');
    end 
else 
    disp('Error: Wrong setting of the operation, please choose op={0,1,2}.'); out=[]; return
end
% 2008 - 2013 (c) TUT v. 2.0