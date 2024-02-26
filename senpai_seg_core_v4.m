function senpai_seg_core_v4(path_in,im_in,varargin)

    % senpai_seg_core:
    %     segments a 3D .tiff image by using a topology-informed
    %     k-means based clustering approach
    %     
    %     Execute the function in the command window:
    %     Syntax:
    %     
    %       senpai_seg_core(path_in,im_in): if you don't choose the
    %                                       path_out, a window will be opened where you can select the
    %                                       final folder for the outputs.
    %      
    %       senpai_seg_core(path_in,im_in,path_out)
    %     
    %       senpai_seg_core(path_in,im_in,path_out,sig_G)
    %     
    %       senpai_seg_core(path_in,im_in,path_out,sig_G,size_lim)
    %     
    %       senpai_seg_core(path_in,im_in,path_out,sig_G,size_lim,verbmem)
    %
    %       senpai_seg_core(path_in,im_in,path_out,sig_G,size_lim,verbmem,paralpool)
    %
    %       senpai_seg_core(path_in,im_in,path_out,sig_G,size_lim,verbmem,paralpool,clusters)
    %     
    %     
    %     Mandatory inputs:
    %       path_in:  the path at which the file im_in is found.
    %     
    %       im_in:    is a string specifying the file name of the input
    %                 image (.tiff file type, either 8 or 16 bits)
    %     
    %     Elective inputs:
    %       path out (char):    path in which partial and final results are saved.
    %                           Default is ./senpaiseg/
    %     
    %       sig_G (numerical):  scalar or 1x2 numerical array.
    %                           The segmentation is repeated length(sig_G) times.
    %                           Each time, derivatives are computed on im_in smoothed
    %                           with a 3d Gaussian kernel having sigma = sig_G(p).
    %                           Default is sig_G=[0], advised for 40x images;
    %                           sig_G=[0 3] advised for 93x images.
    %                           If present, sig_G(2) must be > 0
    %                           sig_G(1)=-1 applies a 3x3x3 median filter
    %     
    %       size_lim (numerical 1x3 array):   maximum crop size on which the
    %                                         algorithm is allowed to work.
    %                                         The image will be divided in crops to save memory.
    %                                         size_lim should be set according to the specifics of your machine.
    %                                         Default is size_lim=[1024 1024 10];
    %     
    %        verbmem (boolean):    if 1, keep partial results (slabs-specific), else,
    %                              delete them. Default is 0.
    %
    %        paralpool (boolean):  if 1, starts a parallel pool while
    %                              doing the kmeans. Default is 1.
    %
    %        clusters (numerical): number of classes for kmeans clustering.
    %                               Default is 6
    %     
    %     Outputs:
    %        senpai_final:   file .mat including:
    %                     
    %                             senpai_final: logical matrix with the final segmentation.
    %             
    %                             senpai_KM: numeric matrix with integer values in range (1:K). It is the
    %                                       product of a clustering performed with k-means by the
    %                                       senpai_ algorithm.
    %             
    %                              cIM:   numeric matrix. Image that generated the senpai_KM segmentation
    %                                    (e.g., senpai_KM = senpai_seg_core(path_in,im_in,varargin)
    %             
    %                              Gxx, Gyy, Gzz:   second order derivatives computed on cIM along
    %                                               the three main axes
    %                                           


    % check input arguments
        if nargin <1
            try
                [im_in,path_in]=uigetfile;
                path_out = uigetdir;             
                catch
                warning('no dataset selected')
                return
            end
        end
    
    if nargin<0
        error('not enough input arguments!')
    else
        disp(['input file: ' path_in im_in])
        if ~isfile([path_in im_in])
            error('input file not found')
        end
    
        % collect information from image header
        info1 = imfinfo([path_in im_in]);
        Nz=length(info1);
        Ny=info1(1).Width;
        Nx=info1(1).Height;
        
        % collect voxel size info 
        if ~isempty(info1(1).XResolution)
            xres=1/info1(1).XResolution;
            yres=1/info1(1).YResolution;
            if isfield(info1(1),'ImageDescription')==1
                strinfo=info1(1).ImageDescription;
                um=strfind(strinfo,'PhysicalSizeZ=');
                cut_s=15;
                if isempty(um)
                    um=strfind(strinfo,'spacing=');
                    cut_s=8;
                end
                zres=str2double(regexp(strinfo(um+cut_s:end),'\w.\w+','match','once'));
            else
                zres=1;
            end
            xres=1;yres=1;
        end
        
        
        bitl=info1(1).BitDepth;
        if bitl==16
            tp='single';
        elseif bitl==8
            tp='uint8';
        else
            error('input image must be 8 or 16 bit')
        end
        
    
        %defaults
        if nargin==2
            try
                path_out = uigetdir;
            catch
                path_out=[pwd filesep 'senpaiseg' filesep];
            end
        end

        sig_G=0;
        if Nx*Ny*Nz<20*10^6
            size_lim=[Nx Ny Nz];
        else
            size_lim=[min(1024,Nx) min(1024,Ny) max(10,floor(10^7/(min(1024,Nx)*min(1024,Ny))))];
        end
        verbmem=0;
        paralpool=1;
        clusters=6;
    end
    
    % parse optional inputs
    if nargin>2
        path_out=varargin{1};
    end
    
    if nargin>3
        sig_G=varargin{2};
        if ~isnumeric(sig_G)
            error('non numeric sigma value for variable sig_G')
        elseif length(sig_G)>2
            error('variable sig_G must be a 1 or 2 elements array')
        elseif length(sig_G)==2 && sig_G(2)<=0
            error('second element in sig_G should be greater than 0')
        end
    end
    
    if nargin>4
        size_lim=varargin{3};
        if ~isnumeric(size_lim)
            throw error
        elseif sum(size_lim<0)>0
            error('no negative values are allowed for the maximum slab size')
        end
        %if limit exceeds size or limit is close to size, limit=size.
        size_lim(1)=min(Nx,size_lim(1));if abs(size_lim(1)-Nx)<10;size_lim(1)=Nx;end
        size_lim(2)=min(Ny,size_lim(2));if abs(size_lim(2)-Ny)<10;size_lim(2)=Ny;end
        size_lim(3)=min(Nz,size_lim(3));if abs(size_lim(3)-Nz)<10;size_lim(3)=Nz;end
    end
    
    if nargin>5
        verbmem=varargin{4};
        if ~islogical(verbmem)
            verbmem=verbmem==1;
        end
    end

    if nargin>6
        paralpool=varargin{5};
        if ~islogical(paralpool)
            paralpool=paralpool==1;
        end
    end

    if nargin>7
        clusters=varargin{6};
        if ~isnumeric(clusters)
            clusters=clusters==2;
        end
    end
    
    [dirsucc,dirmess]=mkdir(path_out);
    if dirsucc==0
        error(dirmess)
    else
        cd(path_out)
    end
    

    %% Core of segmentation algorithm
    % define crop size (heuristic)
    nxiL = 1:size_lim(1):Nx;
    nxeL = min(Nx,size_lim(1):size_lim(1):Nx+size_lim(1)-1);
    nyiL = 1:size_lim(2):Ny;
    nyeL = min(Ny,size_lim(2):size_lim(2):Ny+size_lim(2)-1);
    th_back=0.02*(2^bitl);
    % other crop parameters
    win     = size_lim(3);  % crop slice
    safe    = 3;            % z margin
    safe_xy = 3;
    k_seq   = 1:win:Nz;     % Nz number of slice of the initial image
    
    % start segmentation
    % x axis loop
    for nxin=1:length(nxiL)
        % define indexes for the crop
        nxi  = nxiL(nxin);
        nxe  = nxeL(nxin);
        nxiS = max(1,nxi-safe_xy);
        nxeS = min(Nx,nxe+safe_xy);
        NxC  = nxe-nxi+1;
    
        % y axis loop
        for nyin=1:length(nyiL)
            % define indexes for the crop
            nyi  = nyiL(nyin);
            nye  = nyeL(nyin);
            nyiS = max(1,nyi-safe_xy);
            nyeS = min(Ny,nye+safe_xy);
            NyC  =nye-nyi+1;
    
            % z axis loop
            for k=k_seq
                % define indexes for the crop
                xinit = max(1,k-safe);
                xend  = min(Nz,k+win+safe);
                NzC =min(Nz,k+win-1)-k+1;


                %check if whether crop has already been processed, in that
                %case, continue...
                if isfile([path_out filesep 'sl' mat2str(k) '_' mat2str(nxi) '_' ...
                    mat2str(nxe) '_' mat2str(nyi) '_' mat2str(nye) '.mat'])
                    disp('...already done!')
                    continue
                end
    
                % display status message
                disp(['number of clusters is:' mat2str(clusters)]);
                disp(['crop is x(' mat2str(nxi) ':' mat2str(nxe) '), y(' mat2str(nyi) ':' mat2str(nye) ')']);
                disp(['starting with slice ' mat2str(k) ' over ' mat2str(Nz) ', window set to ' mat2str(win)]);
    
                % read crop of image (i) allows handlig superbig images; (ii)
                % at first define actual size of crop, then (iii) iterates on slices
                nx   = length(nxiS:nxeS);
                ny   = length(nyiS:nyeS);
                nz   = xend-xinit+1;
                cIMc = zeros(nx,ny,nz,tp);
    
                % main loop
                for zz=1:size(cIMc,3)
                    cIMc(:,:,zz) = imread([path_in im_in],xinit+zz-1,...
                        'PixelRegion',{[nxiS nxeS],[nyiS nyeS]});
                end
                
                % median filter if sig_G=-1
                if sig_G(1)==-1
                    disp('applying 3x3x3 median filter...')
                    cIMc=medfilt3(cIMc);
                end
    
                % give some feedback to user
                disp(['1st of ' mat2str(length(sig_G)) ' passes...'])
    
                % 1st-order derivatives
                if sig_G(1)>0
                    [Gx2, Gy2, Gz2] = imgradientxyz(imgaussfilt3(single(cIMc),[sig_G(1) sig_G(1)*xres/yres sig_G(1)*xres/zres]));
                else
                    [Gx2, Gy2, Gz2] = imgradientxyz(single(cIMc));
                end
    
                % 2nd-order derivatives
                [Gxx2, ~, ~] = imgradientxyz(Gx2);
                [~, Gyy2, ~] = imgradientxyz(Gy2);
                [~, ~, Gzz2] = imgradientxyz(Gz2);
    
                % crop out the borders (that were included to avoid edge
                % effects on derivatives computation)
    
                % on cIMc
                cIMc  = cIMc(min(safe_xy+1,nxi):min(NxC+min(safe_xy+1,nxi)-1,size(cIMc,1)),...
                    min(safe_xy+1,nyi):min(NyC+min(safe_xy+1,nyi)-1,size(cIMc,2)),...
                    min(safe+1,k):min(win+min(safe+1,k)-1,size(cIMc,3)));
    
                % on GxxKt
                GxxKt = Gxx2(min(safe_xy+1,nxi):min(NxC+min(safe_xy+1,nxi)-1,...
                    size(Gxx2,1)),min(safe_xy+1,nyi):min(NyC+min(safe_xy+1,nyi)-1,...
                    size(Gxx2,2)),min(safe+1,k):min(win+min(safe+1,k)-1,size(Gxx2,3)));
    
                % on GyyKt
                GyyKt = Gyy2(min(safe_xy+1,nxi):min(NxC+min(safe_xy+1,nxi)-1,size(Gyy2,1)),...
                    min(safe_xy+1,nyi):min(NyC+min(safe_xy+1,nyi)-1,...
                    size(Gyy2,2)),min(safe+1,k):min(win+min(safe+1,k)-1,size(Gyy2,3)));
    
                % and on GzzKt
                GzzKt = Gzz2(min(safe_xy+1,nxi):min(NxC+min(safe_xy+1,nxi)-1,...
                    size(Gzz2,1)),min(safe_xy+1,nyi):min(NyC+min(safe_xy+1,nyi)-1,...
                    size(Gzz2,2)),min(safe+1,k):min(win+min(safe+1,k)-1,size(Gzz2,3)));
    
                % save current image crop
                save([path_out filesep 'sl' mat2str(k) '_' mat2str(nxi) '_' ...
                    mat2str(nxe) '_' mat2str(nyi) '_' mat2str(nye) '.mat'],...
                    'cIMc','GxxKt','GyyKt','GzzKt')
                
                % mask close-to-zero-intensity background to reduce
                % memory usage
                mask_back=find(cIMc(:)>=th_back);
                resiz_km=ones(length(cIMc(:)),1,'uint8');
                % define features space for k-means
                km_in1 = [single(cIMc(mask_back)) GxxKt(mask_back) GyyKt(mask_back) GzzKt(mask_back)];
    
                % clean some memory
                clear G* cIMc
    
                % kmeans clustering options
                options = statset('UseParallel',paralpool);
    
                % kmeans clustering
                if size(km_in1,1)<10
                    TOT_KM1=ones(NxC,NyC,NzC,tp);
                else
                    TOT_KM1=kmeans(single(km_in1),clusters,'Replicates',10,'Options',...
                    options,'MaxIter',1000);
                    resiz_km(mask_back)=uint8(reorderKM(TOT_KM1,km_in1(:,1)));
                    TOT_KM1=reshape(resiz_km,...
                    [nxe-nxi+1 nye-nyi+1 length(resiz_km)/((nxe-nxi+1)*(nye-nyi+1))]);
                end
    
                % save results of k-means for the current crop
                save([path_out filesep 'sl' mat2str(k) '_' mat2str(nxi) '_' ...
                    mat2str(nxe) '_' mat2str(nyi) '_' mat2str(nye) '.mat'],...
                    'TOT_KM1','-append')
    
                % clean some memory
                clear km_in1 TOT_KM1 resiz_km
    
                % eventually, go for a second round
                if length(sig_G)>1
                    disp(['2nd of ' mat2str(length(sig_G)) ' rounds...'])
    
                    % initialize output
                    cIMc=zeros(nx,ny,nz,tp);
    
                    % iterate
                    for zz=1:size(cIMc,3)
                        cIMc(:,:,zz)=imread([path_in im_in],xinit+zz-1,...
                            'PixelRegion',{[nxiS nxeS],[nyiS nyeS]});
                    end
    
                    % 1st-order derivative of smoothed image
                    [Gx4, Gy4, Gz4] = imgradientxyz(imgaussfilt3(single(cIMc),[sig_G(2) sig_G(2)*xres/yres sig_G(2)*xres/zres]));
    
                    % 2nd-order derivative
                    [Gxx4, ~, ~]=imgradientxyz(Gx4);
                    [~, Gyy4, ~]=imgradientxyz(Gy4);
                    [~, ~, Gzz4]=imgradientxyz(Gz4);
    
                    % remove borders
                    cIMc   = cIMc(min(safe_xy+1,nxi):min(NxC+min(safe_xy+1,nxi)-1,...
                        size(cIMc,1)),min(safe_xy+1,nyi):min(NyC+min(safe_xy+1,nyi)-1,...
                        size(cIMc,2)),min(safe+1,k):min(win+min(safe+1,k)-1,size(cIMc,3)));
    
                    GxxKt2 = Gxx4(min(safe_xy+1,nxi):min(NxC+min(safe_xy+1,nxi)-1,size(Gxx4,1)),...
                        min(safe_xy+1,nyi):min(NyC+min(safe_xy+1,nyi)-1,size(Gxx4,2)),...
                        min(safe+1,k):min(win+min(safe+1,k)-1,size(Gxx4,3)));
    
                    GyyKt2 = Gyy4(min(safe_xy+1,nxi):min(NxC+min(safe_xy+1,nxi)-1,size(Gyy4,1)),...
                        min(safe_xy+1,nyi):min(NyC+min(safe_xy+1,nyi)-1,size(Gyy4,2)),...
                        min(safe+1,k):min(win+min(safe+1,k)-1,size(Gyy4,3)));
    
                    GzzKt2 = Gzz4(min(safe_xy+1,nxi):min(NxC+min(safe_xy+1,nxi)-1,...
                        size(Gzz4,1)),min(safe_xy+1,nyi):min(NyC+min(safe_xy+1,nyi)-1,...
                        size(Gzz4,2)),min(safe+1,k):min(win+min(safe+1,k)-1,size(Gzz4,3)));
    
                    % append to saved file
                    save([path_out filesep 'sl' mat2str(k) '_' mat2str(nxi) '_' ...
                        mat2str(nxe) '_' mat2str(nyi) '_' mat2str(nye) '.mat'],...
                        'GxxKt2','GyyKt2','GzzKt2','-append')
                    
                    
                    resiz_km=ones(length(cIMc(:)),1,'uint8');
    
                    % define features space for k-means
                    km_in2=[single(cIMc(mask_back)) GxxKt2(mask_back) GyyKt2(mask_back) GzzKt2(mask_back)];
    
                    % clean memory
                    clear G* cIMc
    
                    % kmeans options
                    options = statset('UseParallel',1);
    
                    % run kmeans
                    if size(km_in2,1)<10
                        TOT_KM2=ones(NxC,NyC,NzC,tp);
                    else
                        TOT_KM2=kmeans(single(km_in2),clusters,'Replicates',10,'Options',options,'MaxIter',1000);
                        resiz_km(mask_back)=uint8(reorderKM(TOT_KM2,km_in2(:,1)));
                        TOT_KM2=reshape(resiz_km,[nxe-nxi+1 nye-nyi+1 length(resiz_km)/((nxe-nxi+1)*(nye-nyi+1))]);
                    end
    
                    % save
                    save([path_out filesep 'sl' mat2str(k) '_' mat2str(nxi) '_' ...
                        mat2str(nxe) '_' mat2str(nyi) '_' mat2str(nye) '.mat'],...
                        'TOT_KM2','-append')
                end %end of if length(sig_G)>1
            end %end of for k=k_seq
        end %end of for nyin=1:length(nyiL)
    end %end of for nxin=1:length(nxiL)
    
    % recompose crops
    senpai_recompose(sig_G,Nx,Ny,Nz,size_lim,path_out,verbmem,tp);

    
    
    function senpai_recompose(sig_G,Nx,Ny,Nz,size_lim,path_out,verbmem,tp)
    
    % this function recomposes processed slabs into a single matrix of
    % k-means clusters having the same size as the original image

    disp('composing segmentations...') 
    % step 1: recompose crops into single images
    senpai_KM_lv1 = zeros(Nx,Ny,Nz,'uint8');
    cIM           = zeros(Nx,Ny,Nz,tp);
    
    % for multiple runs
    if length(sig_G)>1
        senpai_KM_lv2=zeros(Nx,Ny,Nz,'uint8');
    end
    
    for k=1:win:Nz
        for x=1:size_lim(1):Nx
            for y=1:size_lim(2):Ny
                tmpcrp=load(['sl' mat2str(k) '_' mat2str(x) '_' mat2str(min(x+size_lim(1)-1,Nx)) '_' mat2str(y) '_' mat2str(min(y+size_lim(2)-1,Ny)) '.mat']);
                senpai_KM_lv1(x:min(x+size_lim(1)-1,Nx),y:min(y+size_lim(2)-1,Ny),k:min(k+win-1,Nz))=tmpcrp.TOT_KM1;
                cIM(x:min(x+size_lim(1)-1,Nx),y:min(y+size_lim(2)-1,Ny),k:min(k+win-1,Nz))=tmpcrp.cIMc;
                if length(sig_G)>1
                    senpai_KM_lv2(x:min(x+size_lim(1)-1,Nx),y:min(y+size_lim(2)-1,Ny),k:min(k+win-1,Nz))=tmpcrp.TOT_KM2;
                end
            end
        end
    end
    
    bitlr=ceil(log2(double(max(cIM(:)))));
    save([path_out filesep 'senpai_final.mat'],'senpai_KM_lv1','cIM')
    if length(sig_G)>1
        save([path_out filesep 'senpai_final.mat'],'senpai_KM_lv2','-append')
    end
    
    
    % step 2: define classes for neural structures: we save some vectors of
    % indexes (senpai_KM_lv1_msk,senpai_KM_lv2_msk) specifying which classes
    % are to be considered of interest
    
    senpai_KM_lv1_msk=zeros(Nx,Ny,Nz,'uint8');
    
    % here we compute the mean value assumed by second order derivatives
    % within each k-means cluster (LEVEL 1)
    for k=1:win:Nz
        for x=1:size_lim(1):Nx
            for y=1:size_lim(2):Ny
                tmpcrp=load(['sl' mat2str(k) '_' mat2str(x) '_' mat2str(min(x+size_lim(1)-1,Nx)) '_' mat2str(y) '_' mat2str(min(y+size_lim(2)-1,Ny)) '.mat'],'GxxKt','GyyKt','GzzKt');
                Gxx_cl=zeros(1,clusters);
                Gyy_cl=zeros(1,clusters);
                Gzz_cl=zeros(1,clusters);
                for clu = 1:clusters
                    Gxx_cl(clu) = mean(tmpcrp.GxxKt(senpai_KM_lv1(x:min(x+size_lim(1)-1,Nx),y:min(y+size_lim(2)-1,Ny),k:min(k+win-1,Nz))==clu));
                    Gyy_cl(clu) = mean(tmpcrp.GyyKt(senpai_KM_lv1(x:min(x+size_lim(1)-1,Nx),y:min(y+size_lim(2)-1,Ny),k:min(k+win-1,Nz))==clu));
                    Gzz_cl(clu) = mean(tmpcrp.GzzKt(senpai_KM_lv1(x:min(x+size_lim(1)-1,Nx),y:min(y+size_lim(2)-1,Ny),k:min(k+win-1,Nz))==clu));
                end
                senpai_KM_lv1_msk(x:min(x+size_lim(1)-1,Nx),y:min(y+size_lim(2)-1,Ny),k:min(k+win-1,Nz))=(find((Gxx_cl<0 & Gyy_cl<0 & Gzz_cl<0)==0,1,'last')+1).*ones(length(x:min(x+size_lim(1)-1,Nx)),length(y:min(y+size_lim(2)-1,Ny)),length(k:min(k+win-1,Nz)));
            end
        end
    end
    
    
    % update saved file
    save([path_out filesep 'senpai_final.mat'],'senpai_KM_lv1_msk',...
        'Gxx_cl','Gyy_cl','Gzz_cl','-append')
    
    % here we compute the mean value assumed by second order derivatives
    % within each k-means cluster (LEVEL 2)
    if length(sig_G)>1
        senpai_KM_lv2_msk=zeros(Nx,Ny,Nz,'uint8');
        for k=1:win:Nz
            for x=1:size_lim(1):Nx
                for y=1:size_lim(2):Ny
                    tmpcrp=load(['sl' mat2str(k) '_' mat2str(x) '_' mat2str(min(x+size_lim(1)-1,Nx)) '_' mat2str(y) '_' mat2str(min(y+size_lim(2)-1,Ny)) '.mat'],'GxxKt2','GyyKt2','GzzKt2');
                    Gxx_cl=zeros(1,clusters);
                    Gyy_cl=zeros(1,clusters);
                    Gzz_cl=zeros(1,clusters);
                    for clu = 1:clusters
                        Gxx_cl(1,clu) = mean(tmpcrp.GxxKt2(senpai_KM_lv2(x:min(x+size_lim(1)-1,Nx),y:min(y+size_lim(2)-1,Ny),k:min(k+win-1,Nz))==clu));
                        Gyy_cl(1,clu) = mean(tmpcrp.GyyKt2(senpai_KM_lv2(x:min(x+size_lim(1)-1,Nx),y:min(y+size_lim(2)-1,Ny),k:min(k+win-1,Nz))==clu));
                        Gzz_cl(1,clu) = mean(tmpcrp.GzzKt2(senpai_KM_lv2(x:min(x+size_lim(1)-1,Nx),y:min(y+size_lim(2)-1,Ny),k:min(k+win-1,Nz))==clu));
                    end
                    senpai_KM_lv2_msk(x:min(x+size_lim(1)-1,Nx),y:min(y+size_lim(2)-1,Ny),k:min(k+win-1,Nz))=(find((Gxx_cl<0 & Gyy_cl<0 & Gzz_cl<0)==0,1,'last')+1).*ones(length(x:min(x+size_lim(1)-1,Nx)),length(y:min(y+size_lim(2)-1,Ny)),length(k:min(k+win-1,Nz)));
                end
            end
        end
        save([path_out '/senpai_final.mat'],'senpai_KM_lv2_msk','Gxx_cl','Gyy_cl','Gzz_cl','-append')
    end
    
    % delete temporary files in case the users chose so..
    if ~verbmem
        delete sl*
    end
    
    %step 3: in case length(sig_G)>1, produce the final "OR" segmentation
    if length(sig_G)>1
        %create logical matrix with final segmentation
        seg = zeros(size(senpai_KM_lv1),'logical');
    
        % assign values
        seg(senpai_KM_lv1>=senpai_KM_lv1_msk)=1;
        seg(senpai_KM_lv2>=senpai_KM_lv2_msk)=1;
    
        %fill holes in the segmentation
        segf=zeros(size(seg),'logical');for zz=1:size(seg,3);segf(:,:,zz)=imfill(seg(:,:,zz),'holes');end
    
        %clean the segmentation from too small clusters
        BW2=zeros(size(cIM),'logical');for zz=1:size(cIM,3);BW2(:,:,zz) = bwpropfilt((segf(:,:,zz)==1 & cIM(:,:,zz)>0.125*(2^(bitlr))) | (segf(:,:,zz)==0 & cIM(:,:,zz)>0.55*(2^(bitlr))),'Area',[7 inf]);end
    
        senpai_final=zeros(size(BW2),'logical');
        for zz=1:size(cIM,3)
            senpai_final(:,:,zz)=imfill(BW2(:,:,zz),'holes');
        end
        for zz=1:size(cIM,2)
            senpai_final(:,zz,:)=imfill(senpai_final(:,zz,:),'holes');
        end
        for zz=1:size(cIM,1)
            senpai_final(zz,:,:)=imfill(senpai_final(zz,:,:),'holes');
        end
    
        % update saved file
        save([path_out filesep 'senpai_final.mat'],'senpai_final','-append')
    else
        senpai_final = zeros(size(senpai_KM_lv1),'logical');
        
        % assign values
        senpai_final(senpai_KM_lv1>=senpai_KM_lv1_msk)=1;

        %new part to take somas
        for zz=1:size(senpai_final,3);senpai_final(:,:,zz)=imfill(senpai_final(:,:,zz) | (senpai_final(:,:,zz)==0 & cIM(:,:,zz)>0.55*(2^(bitlr))),'holes');end
        for zz=1:size(cIM,3)
            senpai_final(:,:,zz)=imfill(senpai_final(:,:,zz),'holes');
        end
        for zz=1:size(cIM,2)
            senpai_final(:,zz,:)=imfill(senpai_final(:,zz,:),'holes');
        end
        for zz=1:size(cIM,1)
            senpai_final(zz,:,:)=imfill(senpai_final(zz,:,:),'holes');
        end
        %end new part to take somas
       
    
        % update saved file
        save([path_out filesep 'senpai_final.mat'],'senpai_final','clusters','-append')
    end %end di if length(sig_G>1)
    disp([im_in ' : DONE.']) 
    end %end di function senpai_recompose
    end %end of function senpai_seg_core
    
    
    function orderedKM=reorderKM(KMclass,weight)
    % function for sorting kmean classes based on image average intensity
    % within each class
        KM_N=length(unique(KMclass));
        KMc_mean=zeros(KM_N,1);
        for kmm=1:KM_N
            KMc_mean(kmm)=mean(weight(KMclass==kmm));
        end
        [~, ikm]=sort(KMc_mean,'ascend');
        orderedKM=zeros(size(KMclass));
        for kmm=1:KM_N
            orderedKM(KMclass==kmm)=KM_N+find(ikm==kmm);
        end
        orderedKM=orderedKM-KM_N;
    end