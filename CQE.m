%Panetta 2012
function out=CQE(Image)
%Image is color image
c1=0.4358;
c2=0.1722;
c3=0.3920;
out=c1*colorfulness1(Image)+c2*sharpness1(Image)+c3*contrast1(Image);
out=-1+out*2;
if out>1
    out=1;
elseif out<-1
    out=-1;
end
end
function out=contrast1(Image)
Y=rgb2ycbcr(Image);
G=(Y(:,:,1));
out=AME(im2double(G),3);
end
function AME_out = AME(image,wsz)
if ~isa(image,'double')
    image = double(image);
end

[r,c]=size(image);

if mod(r,wsz)==0&&mod(c,wsz)==0
    k1=r/wsz;
    k2=c/wsz;    
    dimg=image;
else
    k1=ceil(r/wsz);
    k2=ceil(c/wsz);
    rp=wsz*k1-r;
    cp=wsz*k2-c;
    dimg=padarray(image,[rp cp]);
end
k=k1*k2;

AME_tmp=0;

for i=1:k1
   for j=1:k2
      block=dimg((i-1)*wsz+1:i*wsz,(j-1)*wsz+1:j*wsz);
      %A=nonzeros(block);
      Imax=max(max(block));
      Imin=min(min(block));
      if Imin == Imax || Imin == -Imax || isnan(Imax) || isnan(Imin)||Imax==0||Imin==0
        k=k-1;
      else
          tmp=((Imax+Imin)/(Imax-Imin));
          AME_tmp_plus=(log10(tmp))^(-0.5);
          AME_tmp=AME_tmp+AME_tmp_plus;
      end
   end
end
AME_out=AME_tmp/k;
end
function out=colorfulness1(Image)
if ~isa(Image,'double')
   Image= double(Image);
end
R=Image(:,:,1);
G=Image(:,:,2);
B=Image(:,:,3);

alpha=R-G;
beta=0.5*(R+G)-B;

Mean_alpha=mean2(alpha);
Mean_beta=mean2(beta);

Var_alpha=mean2((alpha.^2-Mean_alpha^2));
Var_beta=mean2((beta.^2-Mean_beta^2));

out=0.02*log(Var_alpha/(abs(Mean_alpha))^0.2)*log(Var_beta/(abs(Mean_beta))^0.2);
end
function out=sharpness1(Image)
if ~isa(Image,'double')
   Image= double(Image);
end
R=Image(:,:,1);
G=Image(:,:,2);
B=Image(:,:,3);
RS=SW1(R);
GS=SW1(G);
BS=SW1(B);
out=0.299*EMEN(RS.*R,3)+0.587*EMEN(GS.*G,3)+0.114*EMEN(BS.*B,3);
end
function EME_out=EMEN(image,ws)
if ~isa(image,'double')
    image = double(image);
end
image = double(image) ; 
[m,n] = size(image);
map = zeros(m,n);

Imin = ordfilt2(image,1,ones(ws,ws),'symmetric');
Imax = ordfilt2(image,ws^2,ones(ws,ws),'symmetric');
map=log(Imax./Imin);

index = (Imax == 0) | (Imin == 0)| Imin==Imax;
map(index) = 0 ; 
%EME_tmp_plus=log(tmp);
EME_out=2*mean2(map);
end

function map=SW(x)
    x = double(x) ; 
    [m,n] = size(x);
    map = zeros(m,n);

    gv = imfilter(x,(fspecial('sobel'))','replicate');
    gh = imfilter(x,(fspecial('sobel')),'replicate');
    mv = imfilter(x,abs(fspecial('sobel'))','replicate');
    mh = imfilter(x,abs(fspecial('sobel')),'replicate');
    
    % zero background luminance
    index = (mv == 0) & (mh == 0) ;
    map(index) = 0 ; 
    % zero vertical background luminance
    index = (mv == 0) & (mh ~= 0); 
    map(index) = abs(gh(index)./mh(index))/2;
    % zero horizontal background luminance
    index = (mv ~= 0) & (mh == 0); 
    map(index) = abs(gv(index)./mv(index))/2;
    % else
    index = (mv ~= 0) & (mh ~= 0); 
    map(index) = (abs(gv(index)./mv(index))+abs(gh(index)./mh(index)))/2;
end

% function EME_out=EME(image,wsz)
% if ~isa(image,'double')
%     image = double(image);
% end
% 
% [r,c]=size(image);
% padding image with zeros if the row or column is not multiple of window
% size
% if mod(r,wsz)==0&&mod(r,wsz)==0
%     k1=r/wsz;
%     k2=c/wsz;    
%     dimg=image;
% else
%     k1=ceil(r/wsz);
%     k2=ceil(c/wsz);
%     rp=wsz-rem(r,wsz);
%     cp=wsz-rem(c,wsz);
%     dimg=padarray(image,[rp cp]);
% end
% k=k1*k2;
% 
% EME_tmp=0;
% 
% for i=1:k1
%    for j=1:k2
%       block=dimg((i-1)*wsz+1:i*wsz,(j-1)*wsz+1:j*wsz);
%       A=nonzeros(block);
%       Imax=max(A);
%       Imin=min(A);
%       if isempty(A) || isnan(Imax) || isnan(Imin)
%           k=k-1;
%       else
%           tmp=(Imax/Imin);
%           EME_tmp_plus=log(tmp);
%           EME_tmp=EME_tmp+EME_tmp_plus;
%       end
%    end
% end
% 
% EME_out=2*EME_tmp/k;
% end
% function map=SW1(x)
%     x = double(x) ; 
%     [m,n] = size(x);
%     map = zeros(m,n);
% 
%     gv = imfilter(x,(fspecial('sobel'))'/8,'replicate');
%     gh = imfilter(x,(fspecial('sobel'))/8,'replicate');
%     mv = imfilter(x,abs(fspecial('sobel'))'/8,'replicate');
%     mh = imfilter(x,abs(fspecial('sobel'))/8,'replicate');
%     map = (abs(gv)+abs(gh))./(mv+mh);
%     % zero background luminance
%     index = (mv == 0) & (mh == 0) ;
%     map(index) = 0 ; 
% 
% end
function out=SW1(x);
    gh = imfilter(x,[1 2 1 ; 0 0 0; -1 -2 -1],'replicate');
    gv = imfilter(x,[1 0 -1;2 0 -2; 1 0 -1],'replicate');
    %mv = imfilter(x,abs(fspecial('sobel'))'/8,'replicate');
    %mh = imfilter(x,abs(fspecial('sobel'))/8,'replicate');
    S=max(abs(gh),abs(gv));
    out = S./max(S(:));
end
%Panetta 2012
% function out=CQE(Image)
%Image is color image
% % c1=0.4358;
% % c2=0.1722;
% % c3=0.3920;
% %Gaussian Blur
% c1=0.5002;
% c2=0.2448;
% c3=0.2549;
% out=c1*colorfulness1(Image)+(c2+c3)*contrast1(Image);%+c2*sharpness1(Image);
% out=-1+out*2;
% if out>1
%     out=1;
% elseif out<-1
%     out=-1;
% end
% end
% function out=contrast1(Image)
% Y=rgb2ycbcr(Image);
% G=(Y(:,:,1));
% out=AME(im2double(G),3);
% end
% function AME_out = AME(image,wsz)
% if ~isa(image,'double')
%     image = double(image);
% end
% 
% [r,c]=size(image);
% 
% if mod(r,wsz)==0&&mod(c,wsz)==0
%     k1=r/wsz;
%     k2=c/wsz;    
%     dimg=image;
% else
%     k1=ceil(r/wsz);
%     k2=ceil(c/wsz);
%     rp=wsz*k1-r;
%     cp=wsz*k2-c;
%     dimg=padarray(image,[rp cp]);
% end
% k=k1*k2;
% 
% AME_tmp=0;
% 
% for i=1:k1
%    for j=1:k2
%       block=dimg((i-1)*wsz+1:i*wsz,(j-1)*wsz+1:j*wsz);
%       %A=nonzeros(block);
%       Imax=max(max(block));
%       Imin=min(min(block));
%       if Imin == Imax || Imin == -Imax || isnan(Imax) || isnan(Imin)||Imax==0||Imin==0
%         k=k-1;
%       else
%           tmp=((Imax+Imin)/(Imax-Imin));
%           AME_tmp_plus=(log10(tmp))^(-0.5);
%           AME_tmp=AME_tmp+AME_tmp_plus;
%       end
%    end
% end
% AME_out=AME_tmp/k;
% end
% function out=colorfulness1(Image)
% if ~isa(Image,'double')
%    Image= double(Image);
% end
% R=Image(:,:,1);
% G=Image(:,:,2);
% B=Image(:,:,3);
% 
% alpha=R-G;
% beta=0.5*(R+G)-B;
% 
% Mean_alpha=mean2(alpha);
% Mean_beta=mean2(beta);
% 
% Var_alpha=mean2((alpha.^2-Mean_alpha^2));
% Var_beta=mean2((beta.^2-Mean_beta^2));
% 
% out=0.02*log(Var_alpha/(abs(Mean_alpha))^0.2)*log(Var_beta/(abs(Mean_beta))^0.2);
% end