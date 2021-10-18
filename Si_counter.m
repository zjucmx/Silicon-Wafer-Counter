function [si,img_out] = Si_counter(img_in)
%SI_COUNTER �Թ�Ƭ���м���
img=img_in;
img=double(img(:,:,1))/255;
[m,n]=size(img);
% �������˲�
h = fspecial('gaussian',[floor(m/100),floor(n/100)],1);
img = imfilter(img,h);

% ̬ͬ�˲������Ⱦ��⻯
img1=log(0.01+img);
imgfft=fftshift(fft2(img1));
h1=zeros(m,n);
for i=1:m
    for j=1:n
        h1(i,j)=1-exp(-((i-m/2)^2+(j-n/2)^2)/100^2);
    end
end
h1=0.2+2*h1;
imgfft=imgfft.*h1;
img2=exp(real(ifft2(ifftshift(imgfft))))-0.01;
img2(img2>1)=1;
img2(img2<0)=0;

% Sobel���ӣ���ȡ��ֱ��Ե
img2=imfilter(img2,fspecial('gaussian',5,0.5));
img3=imfilter(img2,[1,0,-1;2,0,-2;1,0,-1]);
% OSTU����ֵ������->�㽻�棩
img3=img3>graythresh(img3);

% Hough�任ȷ��ֱ�߷���
BW=img3(floor(m*3/8):ceil(m*5/8),floor(n*3/8):ceil(n*5/8));
[H,T,R]=hough(BW,'RhoResolution',1,'ThetaResolution',0.2);
P=houghpeaks(H,10,'threshold',ceil(0.3*max(H(:))));
X=T(P(:,2));
Y=R(P(:,1));
lines = houghlines(BW,T,R,P,'FillGap',5,'MinLength',m/8);

% ��תͼ��
T_av=mean(X);
img=imrotate(img,T_av,'crop');

% ̬ͬ�˲������Ⱦ��⻯�����⣺���������������
img1=log(0.01+img);
imgfft=fftshift(fft2(img1));
h1=zeros(m,n);
for i=1:m
    for j=1:n
        h1(i,j)=1-exp(-((i-m/2)^2+(j-n/2)^2)/100^2);
    end
end
h1=0.2+2*h1;
imgfft=imgfft.*h1;
img2=abs(exp(real(ifft2(ifftshift(imgfft))))-0.01);
img2(img2>1)=1;
img2(img2<0)=0;

% Ƶ���˲������������ߣ����⣺����α��Ե��
fft=fftshift(fft2(img2));
h2=zeros(m,n);
for i=1:m
    h2(i,:)=exp(-(i-m/2)^2);
end
fft=fft.*h2;
img2=abs(real(ifft2(ifftshift(fft))));

% Sobel���ӣ���ȡ��ֱ��Ե
% img2=imfilter(img2,fspecial('gaussian',5,0.5));
img3=imfilter(img2,[1,0,-1;2,0,-2;1,0,-1]);
% ����
img3=abs(img3);
img3=imclose(img3,ones(5,5));
img3=img3.*(1-img2);
img3=img3/max(max(img3));
img3(img3<0.1)=0;

% ��Ϣѹ��
img4=zeros(m,n);
for i=1:n
    img4(m/4:m*3/4,i)=mean(img3(:,i));
end
% img4(img4<0.1)=0;
img5 = img4;
img5=img5/max(max(img5));
img5=imclose(img5>0.1,ones(5,5));

% ͨ���������½����ж���Ƭλ��
dif_img = diff(mean(img5));
[~,p_loc] = findpeaks(dif_img); % ΢�ַ�ֵΪ������
[~,v_loc] = findpeaks(-dif_img);% ΢�ֹ�ֵΪ�½���
% ����½��غ�������֮���ƽ���Ҷ��ж��Ƿ�Ϊ��Ƭ
si.ind = [];
si.start = [];
si.end = [];
for i=1:length(v_loc)
    start_pix = v_loc(i);
    end_pix = min(p_loc(p_loc>start_pix));
    mean_gray = mean(mean(img(:,start_pix:end_pix)));
    if mean_gray>mean(mean(img))*0.5
        si.ind = [si.ind, length(si.ind)+1];
        si.start = [si.start, start_pix];
        si.end = [si.end, end_pix];
    end
end
img_out = img;
end

