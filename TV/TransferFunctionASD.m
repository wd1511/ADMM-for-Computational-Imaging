function [S,t] = TransferFunctionASD(z,lambda,deltaA,deltaB,M1,M2,c,NM1,NM2,pad)
% TransferFunctionASD() calculates the conventional transfer function of the angular spectrum deconvolution (ASD) 
% S[fy,fx]=exp((i*2*pi*z/lambda)*sqrt(1-lambda^2*((fx/deltaA/M1)^2+(fy/deltaB/M2)^2))).
% 
% Notation:
% [S,t] = TransferFunctionASD(z,lambda,deltaA,deltaB,M1,M2,c,NM1,NM2,pad)
% ------------------------------- input -----------------------------------
% z                     propagation distance between the object and sensor planes [m]
% lambda                wavelength [m]
% deltaA                pitch size w.r.t. rows of the discrete object/sensor plane wave field [m]
% deltaB                pitch size w.r.t. columns of the discrete object/sensor plane wave field [m] 
% M1                    the highest spatial frequency in the y-direction of the computing transfer function
% M2                    the highest spatial frequency in the x-direction of the computing transfer function
% c                     switching on/off of the centering of the coordinate grid (optional):
%                       if c = 1 (default value) - coordinate grid in the interval -M1/2...M1/2-1,  -M2/2...M2/2-1
%                       if c = 0 - coordinate grid in the interval 0...M1-1,  0...M2-1
% NM1 x NM2             the size of the output transfer function (optional, 
%                       must be not smaller than M1 x M2, by default it is equal to M1 x M2)
% pad                   padding of the calculated tranfer function (optional) if NM1 x NM2 is larger than M1 x M2:
%                       if pad = 0 (default value) - by zero-padding, 
%                       if pad = 1 - calculating ASD for the extended grid to the size NM1 x NM2.
% Note! M1, M2, NM1 and NM2 must be positive even integers!
% ------------------------------- output ----------------------------------
% S                     the output ASD transfer function of the size N1 x N2 
% t                     the computational time in seconds (optional)  
%
% Example #1:
% Calculation of the sensor plane wave field uz using the amplitude object u0 of the size N1 x N2. 
% Here we use rectangular pixels at the object and sensor planes of the same size deltaA x deltaB.
%
% u0 = im2double(imread('Lena256.png'));                 % read the test-image
% [N1,N2,N3] = size(u0);                                 % image size
% lambda = 6328e-10;                                     % wavelength (red laser)
% deltaA = 6e-6; deltaB = 8e-6;                          % rectangular pixel size of the object !!!  
% z = deltaA*deltaA*N1/lambda;                           % propagation distance
% S = TransferFunctionASD(z,lambda,deltaA,deltaB,N1,N2);
% uz = ifft2(fft2(u0).*S);
%
% Example #2:
% Simulation of the diffraction propagation using the 4f optical system with an optical filter encoded onto 
% the spatial light modulator placed at the Fourier plane between two lenses with the focal distance f. 
% Here the parameters for the Fourier plane are calculated via the settings of the optical system 
% 
% u0 = im2double(imread('Baboon512.png'));               % read the test-image
% [N1,N2,N3] = size(u0);                                 % image size
% lambda = 532e-9;                                       % wavelength (green laser)
% delta1 = 3.45e-6;                                      % the size of the square pixels at the object and sensor planes  
% delta2 = 8e-6;                                         % the size of the square pixels at the Fourier planes
% f = 58e-3;                                             % focal distance of the lenses
% NM = ceil(lambda*f/delta1/delta2/2)*2;                 % our computational size of the wave field distribution which 
%                                                        % obeys the sampling condition for a perfect reconstruction (even number)
%                                                        % Note that in this case we calculate the fourier transform for the sparial 
%                                                        % frequencies {fx/lambda/f, fy/lambda/f} with the pitch size delta2 
%                                                        % because 1/delta1/NM = delta2/lambda/f.  
% u0_ = zeros(NM); u0_(NM/2-N1/2+1:NM/2+N1/2,NM/2-N2/2+1:NM/2+N2/2) = u0;
% z = 0.1;                                               % arbitrary propagation distance (set e.g. z = 0 to test the 
%                                                        % 4f propagation model with no optical mask)
% S = TransferFunctionASD(z,lambda,delta1,delta1,NM,NM);% ASD with the object plane pixel size as the input parameters
% uz_ = -exp(i*8*pi/lambda*f)*fft2(fft2(u0_).*S)/NM/NM;
% uz = uz_(NM/2-N1/2+1:NM/2+N1/2,NM/2-N2/2+1:NM/2+N2/2); % the output estimate is the resulting central part
%
% Reference:            [1] J. W. Goodman, Introduction to Fourier Optics, 2nd. ed., (New York: McGraw-Hill, 1996). 
% Authors:              Artem Migukin
% Software version:     2.0, January 2013.

% --------------- please do not change the code under this line -----------------------------------------
if nargin < 7 || isempty(c),       c = 1;end                                % centering of the coordinate grid
if (M1<1 || M2<1) || (ceil(M1)-floor(M1)==1 || ceil(M2)-floor(M2)==1)
   disp('Error: the higher normalized frequency of the transfer function must be a positive integer!');   S = [];t=0; return
end
if mod(M1,2)~=0 || mod(M2,2)~=0
    disp('Error: the higher normalized frequency of the transfer function must be an even integer!');  S = [];t=0; return
end
if nargin < 8 || isempty(NM1),     NM1 = M1;end
if nargin < 9 || isempty(NM2),     NM2 = M2;end

if (NM1<1 || NM2<1) || (ceil(NM1)-floor(NM1)==1 || ceil(NM2)-floor(NM2)==1)
   disp('Error: the size of the transfer function must be a positive integer!');   S = [];t=0; return
end
if mod(NM1,2)~=0 || mod(NM2,2)~=0
    disp('Error: the size of the transfer function must be an even integer!');  S = [];t=0; return
end
if NM1<M1, disp('Error: output size NM1 must be not smaller than M1. Use the size M1'); NM1 = M1;end
if NM2<M2, disp('Error: output size NM2 must be not smaller than M2. Use the size M2'); NM2 = M2;end
if nargin < 10 || isempty(pad), pad = 0; end                                % default zero-padding if NM1 > M1 and/or NM2 > M2 
if z ~=0 && abs(z)<lambda
   disp('Error: too small propagation distance, "z" must be larger "lambda"!') 
   S = []; return
elseif z==0,   S = ones(NM1,NM2); return 
end

i = sqrt(-1);
if nargout > 1,    tic;end                                                  % start a stopwatch timer 
if pad==0                                                                   % zero-padding
    S = zeros(NM1,NM2);
    % sampling grid
    if c==1,  [x,y]=meshgrid(-M2/2:M2/2-1,-M1/2:M1/2-1);
    else      [x,y]=meshgrid(0:M2-1,0:M1-1);    end
    SS = exp((i*2*pi*z/lambda)*sqrt(1-(lambda/deltaA*y/M1).^2-(lambda/deltaB*x/M2).^2));
    S(NM1/2-M1/2+1:NM1/2+M1/2,NM2/2-M2/2+1:NM2/2+M2/2) = SS;
else                                                                        % extended grid
    % sampling grid
    if c==1,  [x,y]=meshgrid(-NM2/2:NM2/2-1,-NM1/2:NM1/2-1);
    else      [x,y]=meshgrid(0:NM2-1,0:NM1-1);    end
    S = exp((i*2*pi*z/lambda)*sqrt(1-(lambda/deltaA*y/M1).^2-(lambda/deltaB*x/M2).^2));    
end

if c==1                                                                     % the transfer function is shifted to obtain the 0th frequency
                                                                            % in the upper left corner (for the correct FFT)
    S = circshift(S,-round([NM1,NM2]/2));
end
if nargout > 1
    t=toc; % stop the stopwatch timer - obtain the time of calculation
end
% 2010-2013 (c) TUT v. 2.0