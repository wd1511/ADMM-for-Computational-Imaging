function z = Ep(p)
[m,n,q]=size(p);
z = zeros(m,n,4);
z(:,:,1) = Dx(p(:,:,1));
z(:,:,2) = Dy(p(:,:,2));
z(:,:,3) = (Dy(p(:,:,1))+Dx(p(:,:,2)))./2;
z(:,:,4) = z(:,:,3);

