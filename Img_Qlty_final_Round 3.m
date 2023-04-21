clear all;
close all;
clc;
Vptr = VideoReader('E:\SS584-Round-1\Video010M7001.ts');        
% f1 = read(Vptr,[25*46 25*46+2000]);
totalFrames = Vptr.NumFrames
% f1 = read(Vptr,[25*46 25*46+2000]);
f1 = read(Vptr,[1 totalFrames]);
co=0;
for i=1:750:totalFrames     
        co=co+1;
        VQM(co,1)=CQE((f1(:,:,:,i)));
end
t = table(VQM);
writetable(t,'ss584_25045_runnumber_VideoQC.csv');
tic
pause(1)
fileID = fopen('ss584_25045_runnumber_time','w');
fprintf(fileID,'%4f',toc);
% clear all;
% close all;
% clc;
% Vptr = VideoReader('E:\SS584_Video_Quality\SS584-Video-Quality\M1-Videos\Video001M1001.ts');        
% % f1 = read(Vptr,[25*46 25*46+2000]);
% totalFrames = Vptr.NumFrames;
% % d = duration(H,MI,S);
% % f1 = read(Vptr,[25*46 25*46+2000]);
% f1 = read(Vptr,[25 totalFrames]);
% co=0;
% %VQM(co,1)=0;
% for i=1:750:totalFrames
%         co=co+1;
%         VQM(co,1)=CQE((f1(:,:,:,i)));
% end
% t = table(VQM);
% writetable(t,'ss584_25045_runnumber_VideoQC.csv');
% tic
% pause(1)
% fileID = fopen('ss584_25045_runnumber_time','w');
% fprintf(fileID,'%4f',toc);
clear all;
close all;
clc;
Vptr = VideoReader('E:\SS584_Video_Quality\SS584-Video-Quality\Ref-videos\Video004.mp4');        
% f1 = read(Vptr,[25*46 25*46+2000]);
totalFrames = Vptr.NumFrames;
% f1 = read(Vptr,[25*46 25*46+2000]);
f1 = read(Vptr,[25 totalFrames]);
co=0;
% while hasFrame(Vptr)
%     frame = readFrame(Vptr);
%     pause(1/Vptr.totalFrame);
% end

for i=1:750:totalFrames     
        co=co+1;
        VQM(co,1)=CQE((f1(:,:,:,i)));
end
t = table(VQM);
writetable(t,'ss584_25045_runnumber_VideoQC.csv');
tic
pause(1)
fileID = fopen('ss584_25045_runnumber_time','w');
fprintf(fileID,'%4f',toc);