%% To test the toolbox on the test image "test_data.tif":
%% download the folder and cd in it, then run this script
%% we recommend reading the README file.
clear
clc
close all

current_fold=pwd;

% TO PERFORM ESTIMATION OF OPTIMAL K uncomment lines 14-15-16
% [120 370 120 370 30 100] is a representative crop for
% test_data.tif
Kopt=6;
% mkdir([pwd filesep 'estK' filesep]);
% Kopt = senpai_estimateK([pwd filesep],'test_data.tif',[pwd filesep 'estK' filesep],[120 370 120 370 30 100]);
% cd(current_fold)

% segmentation step
senpai_seg_core_v4([current_fold filesep],'test_data.tif',[current_fold filesep 'res_test_full' filesep],0,[256 256 149],[],[],Kopt);
load('senpai_final.mat', 'senpai_final')
load('senpai_final.mat', 'cIM')
cd ..
% load inputs for the parcellation step: these include markers of somas 
% and additional markers manually set with senpai_prune 
load('somas.mat', 'somas')
load('markers.mat', 'markers')
senpai_separator(senpai_final,cIM,somas | markers);
load('senpai_separator.mat', 'parcel_final');
% visualization of a selection of close-by neurons
sel=[73 30 43 12 13 11 16 10 8 6];
figure;
for ss=sel
[p,v]=isosurface(parcel_final==ss,0.2);
hold on;patch('Faces',p,'Vertices',v,'FaceColor',rand(1,3),'EdgeColor','none','FaceLighting','gouraud');
end
axis equal;box off;axis off;material dull;camlight headlight;camlight headlight;
title('selection of close-by neurons')
view(60,-20)

% eventually, check neurons one by one and perform manual corrections on the parcellation 
% manual correction of the parcellation
%load('senpai_separator.mat', 'parcel_final')
%senpai_prune(parcel_final,10,somas | markers)