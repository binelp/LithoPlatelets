clear
close all

% RA is a imref2d object which is used to store the position at which the
% image is located. One has to run inspectStereoImage with the camA image,
% put a breakpoint where newImg_W and newImg_H are caculated, then run the
% followig snippet to save RA
% % RA = imref2d(size(newImg),[(2448-newImg_W)/2, (2448+newImg_W)/2],-[-(2048-newImg_H)/2, -(2048+newImg_H)/2])
% % save('RA.mat','RA')
% Then repeat the procedure with the camB image and save RB this time
% RA and RB together allow creating a combo image which will perfectly
% superimpose with both contours of camera A and B

load('RB');
load('RA');
a = imread('camA.png');
b = imread('camB.png');
overimg = imshowpair(a,RA,b,RB,'method','diff');
% newImg = rgb2gray(overimg.CData);
newImg = overimg.CData;
newImg = uint8(255) - newImg;
low = 0.3;
high = 0.7;
newImg = 0.5.*imadjust(newImg,[0.5 0.98]); % I is double

close all
pStruct = inspectStereoImage('/imageData/binelp/20210903_200x100x30_TH33/03092021-160516/data_camB/B-03092021-160516-60.csv',...
'imgdir','/imageData/binelp/20210903_200x100x30/03092021-160516/data_camB',...
'pid',1,...
'newImg',newImg);