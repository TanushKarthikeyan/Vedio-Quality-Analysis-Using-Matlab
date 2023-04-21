clear all;
close all;
clc;
Vptr = VideoReader('G:\Image_Quality_SIH\SS584-Video-Quality\Ref-videos\Video001.mp4');
Vptr1 = VideoReader('G:\Image_Quality_SIH\SS584-Video-Quality\M1-Videos\Video001M1001.ts');

Nrows = Vptr.height;
Ncols = Vptr.width;

f1 = im2double(read(Vptr,25*46));
f2 = read(Vptr1,[25*4 25*5+1200]);
min_err=50000;
co=0;
for i = 1:1200
mae1=sum(sum(sum(abs(f1-im2double(f2(:,:,:,i)))*255)))/(Nrows*Ncols*3);
[mssim, ssim_map] = ssim(f1,im2double(f2(:,:,:,i)));
if mae1<min_err
    min_err=mae1;
%     figure,imshow(f2(:,:,:,i));
    fx=i;
end
MAE(i,1)=mae1;
SSIM(i,1)=mssim;
end
t = table(MAE,SSIM);
writetable(t,'M.csv');
