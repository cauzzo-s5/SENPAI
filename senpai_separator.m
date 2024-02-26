function parcel_final=senpai_separator(senpai_seg,cIM,somas)
    % senpai_separator:
    %   takes the result of a k-means clustering performed by
    %   senpai_seg_core and produces separated segmentations of single 
    %   neurons. The program is intended to be used when the segmentation 
    %   produced on low-magnification images fails to correctly separate 
    %   neurons.
    %
    %   Execute the function in the command window:
    %   Syntax:
    %       parcel_final = senpai_neurogrow;
    %       parcel_final = senpai_neurogrow(senpai_KM_lv1,cIM,somas)
    %
    %   Inputs:
    %       senpai_seg: logical matrix with the final segmentation produced by senpai_seg_core.m.
    %
    %       cIM:   numeric matrix. Image that generated the senpai_KM_lv1 segmentation
    %              (e.g., senpai_KM_lv1 = senpai_seg_core(path_in,im_in,varargin)
    %
    %       somas:  logical matrix encoding a gross binary segmentation of the somas in the image
    %
    %   Output:
    %       parcel_final: numeric matrix with the parcelation of the final
    %                     Kmeans segmentation. Every value of the matrix correspond to a
    %                     neuron.
    
    %take the negative of the image filtered with a median filter
    disp('Separating neurons...')
    senpai_seg=logical(senpai_seg);
    db=max(cIM(:));
    cIM_inv=db-medfilt3(cIM);
    cIM_inv(~senpai_seg)=db;
    clear cIM
    %impose minima in the mask of somas
    cIM_inv=imimposemin(cIM_inv,somas);
    %watershed transform
    ww=uint16(watershed(cIM_inv)); %must be uint16
    %provide final parcellation
    parcel_final=ww.*uint16(senpai_seg);
    neuLst=1:max(ww(:));
    save senpai_separator.mat ww
    clear ww senpai_seg
    %remove non-connected pieces
    disp('Pruning non-connected branches...')
    for vv=neuLst
        bb=bwconncomp(parcel_final==vv,6);
        [~, ii]=max(cellfun(@length, bb.PixelIdxList));
        clusLst=1:bb.NumObjects;
        parcel_final(cell2mat(bb.PixelIdxList(clusLst(clusLst~=ii))'))=0;
    end
    save senpai_separator.mat parcel_final -append
    disp('Done!')
end