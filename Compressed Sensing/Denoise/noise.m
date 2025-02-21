function [g,E,Phase] = noise(f,nx,ny,nz,lambda,z,deltaX,deltaY,sigma1)
    f = im2uint8(f);
    f = double(f)+sigma1 * randn(size(f));
    f = mat2gray(f);
    Phase=MyMakingPhase(nx,ny,z,lambda,deltaX,deltaY);
    %figure;imagesc(plotdatacube(angle(Phase)));title('Phase of kernel');axis image;drawnow;
    %axis off; colormap(hot); colorbar;
    E0=ones(nx,ny);  % illumination light
    E=MyFieldsPropagation(E0,nx,ny,nz,Phase); 

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
    %g = im2double(g);
    % º”‘Î…˘
    
    g = im2uint8(g);
    g = double(g)+sigma1 * randn(size(g));
    g = mat2gray(g);
    %[ssss,g]=BM3D(o,g,sigma1);
    %g = BM3D_Gray(g, 0, sigma1, 1); 
    %g = im2double(g);