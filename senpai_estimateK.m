function Kopt = senpai_estimateK(path_in, im_in, path_out, crop_idxs)
% senpai_estimateK:
%     produces an estimate for the optimal K parameter for the
%     senpai_seg_core.m function. The estimation is performed by running
%     the segmentation multiple times with different K,
%     until one K satisfies the criteria in the manuscript.
%     When no K in [2 10] satisfies the criteria, 
%     Kopt=10 is provided and a warning is displayed.
%
%     Execute the function in the command window:
%     Syntax:
%
%       senpai_estimateK(path_in, im_in, path_out):
%
%       inputs
%
%       path_in: path for the input file
%
%       im_in: filename for the input file
%
%       path_out: output path for the segmentation tests
%
%       crop_idxs: start and end indexes for a representative crop on which
%           to run the estimation, e.g., [10 200 30 220 1 50] to run the
%           estimation on a crop starting at x=10, y=30 and z=1 and ending 
%           with x=200, y=220 and z=50.

startK = 2;  % minimum K
stopK = 10;  % maximum K

info1 = imfinfo([path_in im_in]);
Nz=length(info1)-3;
Ny=info1(1).Width;
Nx=info1(1).Height;
bitl=info1(1).BitDepth;
if bitl==16
    tp='single';
elseif bitl==8
    tp='uint8';
else
    error('input image must be 8 or 16 bit')
end
cIMc=zeros(crop_idxs(2)-crop_idxs(1)+1,crop_idxs(4)-crop_idxs(3)+1,crop_idxs(6)-crop_idxs(5)+1,tp);
for zz=1:size(cIMc,3)
    cIMc(:,:,zz) = imread([path_in im_in],crop_idxs(5)+zz-1,...
        'PixelRegion',{[crop_idxs(1) crop_idxs(2)],[crop_idxs(3) crop_idxs(4)]});
end
imwrite(cIMc(:,:,1),[path_out filesep 'tmp_crp.tif'],'tif');
for sl=2:size(cIMc,3)
    imwrite(cIMc(:,:,sl),[path_out filesep 'tmp_crp.tif'] ,'WriteMode','append');
end

%% SENPAI LOOP
for ff = startK:stopK
    path_out_cl = [path_out 'k_for_' num2str(ff) '_classes'];
    
    senpai_seg_core_v4(path_out, 'tmp_crp.tif', path_out_cl, 0, [crop_idxs(2)-crop_idxs(1)+1 crop_idxs(4)-crop_idxs(3)+1 crop_idxs(6)-crop_idxs(5)+1], 1, 0, ff);
    load([path_out_cl filesep 'sl1_1_' mat2str(crop_idxs(2)-crop_idxs(1)+1) '_1_' mat2str(crop_idxs(2)-crop_idxs(1)+1) '.mat'], 'TOT_KM1', 'GxxKt', 'GyyKt', 'GzzKt')

    % Analisi distribuzioni
    senpai_KM_lv1 = TOT_KM1;  
    for kk = 2:ff    % Creo matrice contenete sulle righe le derivate e sulle colonne la percentuale di distribuzione positiva
        dist_xx = zeros(1,kk);
        dist_yy = zeros(1,kk);
        dist_zz = zeros(1,kk);
        for tt = 1:kk
            dist_xx(1,tt) = sum(sum(sum(GxxKt(senpai_KM_lv1(:)==tt)>0)))/sum(sum(sum(senpai_KM_lv1(:)==tt,2)));  % Numero delle derivate (per ogni classe) a valore positivo / numero delle derivate della classe 
            dist_yy(1,tt) = sum(sum(sum(GyyKt(senpai_KM_lv1(:)==tt)>0)))/sum(sum(sum(senpai_KM_lv1(:)==tt,3)));
            dist_zz(1,tt) = sum(sum(sum(GzzKt(senpai_KM_lv1(:)==tt)>0)))/sum(sum(sum(senpai_KM_lv1(:)==tt,4)));
        end

        th = 0.75;  % Soglia per capire qunado ci va bene
        dist_tot = [dist_xx; dist_yy; dist_zz];
        save dist.mat dist_tot
        dist_seg = dist_tot > th;
    end

    val_class = sum(sum(dist_seg,1)>0);  % Somma sulle righe (>0 perchè ci va bene che siano anche più di una)
    val_der = sum(sum(dist_seg,2)>0);  % Somma sulle colonne
    if val_class >= 3 && val_der == 3  % Se ci viene 3 in tutte e due vuol dire che abbiamo almeno una distribuzione positiva per ciascuna classe
       Kopt = ff;
       break
    elseif ff == stopK
           Kopt = stopK;
           warning('no K in the range [2 10] satisfied the criteria. Kopt was set to 10.')
    end
end