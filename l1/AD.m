function output=AD(o,u0,K,z1,dz,lambda1,lambda2,lambda3,lambda4,lambda5,lambda6,alpha0,alpha1,beta,maxiter)

[m,n]=size(u0);
u  = zeros(m,n);
du = zeros(m,n,2);
p  = zeros(m,n,2);
x  = zeros(m,n,2);
y  = zeros(m,n,4);
v1 = zeros(m,n,2);
v2 = zeros(m,n,4);

u=o(:,:,1);
du(:,:,1) = Dx(u);
du(:,:,2) = Dy(u);

method = 0;
av = true;mu = 20e-3;
deltaA = 3.1e-6;deltaB = 3.1e-6;
lambda = 532e-9;
gammar = 1.0/K *ones(K,1);

[A,B,strs,tm]=TransformMatrix(method<1,z1,dz,lambda,deltaA,deltaB,m,n,K,gammar,mu,av);

a=alpha1*lambda2;
b=alpha0*lambda4;

D1 = abs(psf2otf([-1,1],[m, n])).^2;
D2 = abs(psf2otf([-1;1],[m, n])).^2;
d1 = a.*(D1+D2)+beta*sum(A.*B,3);
d2 = a+b*(D1+.5*D2);
d3 = a+b*(D2+.5*D1);
d4 = psf2otf([1, -1],[m, n]);
d5 = psf2otf([1; -1],[m, n]);
d6 = d5(:,1)*d5(:,1)';

d4 = -a*d4;
d5 = -a*d5;
d6 = .5*d6;

d4t = conj(d4);
d5t = conj(d5);
d6t = conj(d6);

denom = d1.*d2.*d3+d4t.*d6t.*d5+d4.*d5t.*d6...
    -d2.*d5.*d5t-d1.*d6t.*d6-d3.*d4.*d4t;

for i = 1:maxiter
    %solve x,y
    x = lambda1*shrink2(du-p+v1,lambda1/lambda2);
    Epp = Ep(p);
    y = lambda1*shrink2(Epp+v2,lambda3/lambda4);
    
    %solve u,p
    Ar_u0=B.*fft2(o);
    FB1 = a*(fft2(Dxt(x(:,:,1)-v1(:,:,1))) + fft2(Dyt(x(:,:,2)-v1(:,:,2))))...
        +beta*sum(Ar_u0,3);
    FB2 = a*fft2(v1(:,:,1)-x(:,:,1))+b*(fft2(Dxt(y(:,:,1)-v2(:,:,1)))...
        +fft2(Dyt(y(:,:,3)-v2(:,:,3))));
    FB3 = a*fft2(v1(:,:,2)-x(:,:,2))+b*(fft2(Dyt(y(:,:,2)-v2(:,:,2)))...
        +fft2(Dxt(y(:,:,3)-v2(:,:,3))));
    
    RHS1 = (d2.*d3-d6.*d6t).*FB1-(d3.*d4t-d6.*d5t).*FB2...
        +(d4t.*d6t-d2.*d5t).*FB3;
    u = ifft2(RHS1./denom);
    u = abs(u);
    du(:,:,1) = Dx(u);
    du(:,:,2) = Dy(u);
    
    RHS2 = (d5.*d6t-d3.*d4).*FB1+(d1.*d3-d5.*d5t).*FB2...
        +(d4.*d5t-d1.*d6t).*FB3;
    p(:,:,1) = ifft2(RHS2./denom);
    RHS3 = (d4.*d6-d2.*d5).*FB1+(d4t.*d5-d1.*d6).*FB2...
        +(d1.*d2-d4.*d4t).*FB3;
    p(:,:,2) = ifft2(RHS3./denom);
    
    %solve v1,v2
    v1 = v1 + lambda5*(du-p-x);
    v2 = v2 + lambda6*(Ep(p)-y);
    fprintf('%8.4f\n',ssim(u,u0));
    a=a/1.01;b=b/1.01;
end
output = u;