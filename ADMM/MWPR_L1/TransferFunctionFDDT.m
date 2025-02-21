function [V,t,c] = TransferFunctionFDDT(z,lambda,deltaA,deltaB,NM1,NM2,av)
% TransferFunctionFDDT() calculates the averaged or non-averaged F-DDT [1] 
% transfer function used for the precise forward wave field propagation from 
% the object to the parallel transverse sensor plane at the distance z performed 
% in the frequency domain. The Fresnel approximation of the Rayleigh-Sommerfeld
% integral is used. The (rectangular) pixels at the object and senor planes are
% assumed to be of the same size deltaA x deltaB. 
%
% Notation:
% [V,t,c] = TransferFunctionFDDT(z,lambda,deltaA,deltaB,NM1,NM2,av)
% ------------------------------- input -----------------------------------
% z                     propagation distance between the object and sensor planes [m]
% lambda                wavelength [m]
% deltaA                pitch size w.r.t. rows of the discrete object/sensor plane wave field [m]
% deltaB                pitch size w.r.t. columns of the discrete object/sensor plane wave field [m] 
% NM1                   number of rows of the computing transfer function
% NM2                   number of columns of the computing transfer function 
% av                    switch on/off (av = 1/0) the averaging for the F-DDT 
%                       transfer function (optional, by default av = 1). 
%                       If av = 0, the non-averaged transfer function result 
%                       in the conventional Fresnel discrete diffraction transfrom.
% NM1 and NM2 must be even positive integers!
% ------------------------------- output ----------------------------------
% V                     the averaged or non-averaged F-DDT transfer function of the size NM1 x NM2 
% t                     the computational time in seconds (optional)  
% c                     the conditioning number of the transfer function (optional)
%
% Example:
% Calculation of the sensor plane wave field uz of the size M1 x M2 using the 
% amplitude object u0 of the size N1 x N2. Here we use rectangular pixels at the
% object and sensor planes of the same size deltaA x deltaB.
%
% u0 = im2double(imread('Lena256.png'));                 % read the test-image
% [N1,N2,N3] = size(u0);                                 % image size
% lambda = 6328e-10;                                     % wavelength (red laser)
% deltaA = 6e-6; deltaB = 8e-6;                          % the size of the rectangular pixels at the object and sensor planes
% M1 =  N1; M2 = N2; NM1 = N1 + M1; NM2 = N2 + M2;       % NM1, NM2 are assumed to be even
% u0_ = zeros(NM1,NM2); u0_(NM1/2-N1/2+1:NM1/2+N1/2,NM2/2-N2/2+1:NM2/2+N2/2) = u0;
% z = deltaA*deltaA*M1/lambda;                           % propagation distance
% V = TransferFunctionFDDT(z,lambda,deltaA,deltaB,NM1,NM2);
% uz_ = ifft2(fft2(u0_).*V);
% uz = uz_(NM1/2-M1/2+1:NM1/2+M1/2,NM2/2-M2/2+1:NM2/2+M2/2); % the output estimate is the resulting central part
%
% Reference:            [1] V. Katkovnik, J. Astola, and K. Egiazarian, Appl. Opt. 47, 3481-3493 (2008).
% Authors:              Artem Migukin, Vladimir Katkovnik
% Software version:     2.0, January 2013.

% --------------- please do not change the code under this line -----------------------------------------
if (nargin<7 || isempty(av)),    av=true;      end                          % by default the averaging for
                                                                            % transform matrix is switched on 
if (NM1<1 || NM2<1) || (ceil(NM1)-floor(NM1)==1 || ceil(NM2)-floor(NM2)==1)
    disp('Error: the size of the transfer function must be a positive integer!');  V = [];t=0;c=[]; return
end
if mod(NM1,2)~=0 || mod(NM2,2)~=0
    disp('Error: the size of the transfer function must be an even integer!');  V = [];t=0;c=[]; return
end

if z ~=0 && abs(z)<lambda
   disp('Error: too small propagation distance, "z" must be larger "lambda"!');  V = [];t=0;c=[]; return
elseif z==0
   V = ones(NM1,NM2);t=0;c=1; return 
end

i = sqrt(-1);
if nargout > 1,    tic;end                                                  % start a stopwatch timer 

t=i*pi/z/lambda;                                                            % common exponent
kf = -i*exp(i*2*pi*z/lambda)/(z*lambda)*deltaA*deltaB;                      
xint=-NM2/2:NM2/2-1;yint=-NM1/2:NM1/2-1;                                    % the calculation is performed for a single row
[X,Y]=meshgrid(xint,yint);                                                  % sampling grid
a_ks=zeros(NM1,NM2);
fx = zeros(1,length(xint));fy = zeros(1,length(yint));
for q=1:length(xint)                                                        % we calculate a row of all kernel values
    if av
        fx(q)=quad(@(x)kernel(x,t,xint(q),deltaB),-1,1);                     % averaged values via a numerical integration
    else
        fx(q)=exp(t*deltaB^2*xint(q)^2);                                     % non averaged values
    end
end

if (NM1 == NM2) && (deltaA == deltaB)                                                                     % acceleration is the matrix is square
    smoothed=fx(X+NM2/2+1).*fx(Y+NM1/2+1);    
else
    for q=1:length(yint)
        if av
            fy(q)=quad(@(x)kernel(x,t,yint(q),deltaA),-1,1);
        else
            fy(q)=exp(t*deltaA^2*yint(q)^2);
        end
    end
    smoothed=fx(X+NM2/2+1).*fy(Y+NM1/2+1); 
end   

y1=kf.*smoothed;                                                            % combine all 3 components of the Fresnel approx. 
% the convolution support of two signals of the length N and M is N+M-1
% we construct the transfer function of the size N+M and set the last
% component equalt to zero
a_ks(2:end,2:end)=y1(2:end,2:end);

V = fft2(circshift(a_ks,-1*round([NM1,NM2]/2)));

if nargout > 1
    t=toc; % stop the stopwatch timer - obtain the time of calculation
end
if nargout > 2
    c = max(max(abs(V)))/min(min(abs(V)));                                  % the conditioning for the transfer function
end

% instrumental inline functions for numerical integration 
function y = kernel(x,t,d,delta)
    y = (1-abs(x)).*exp(t*delta^2*(d+x).^2);
% 2008 - 2013 (c) TUT v. 2.0