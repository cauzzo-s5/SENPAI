
function senpai_prune(parcellation,varargin)

% senpai_prune 
%       This function runs a GUI that displays, one by one, all the neurons
%       in the input parcellation. The GUI is intended to allow the user to 
%       assess eventual misassignments. It is possible to mark unrecognized
%       neural cores that can be used in a new run of senpai_separator to 
%       produce a corrected parcellation
%
    %   Execute the function in the command window:
    %   Syntax:
    %       senpai_prune(parcellation);
    %       senpai_prune(parcellation,nn);
    %       senpai_prune(parcellation,nn,markers);
    %
    %   Inputs:
    %       parcellation: uint8/uint16 3D matrix containing a
    %           parcellation as obtained with senpai_separator.m
    %       nn: label index of the neuron that is displayed first (default:
    %           nn=1)
    %       markers: logical 3D matrix (same size as parcellation)
    %           containing the masks of neural cores
    %
    %   Usage of the GUI:
    %   the user may rotate the 3D reconstruction of the neuron and zoom on 
    %   details of interest. To place a marker, place first the data tip,
    %   then click on the "Mark this branch!" button. Navigate through the
    %   neurons of the parcellation with the arrows below the 3D view. Save
    %   before exiting.

N=max(parcellation(:));
tmp=parula(101);
if nargin>1
    nn=varargin{1};
else
    nn=1;
end
neuron=parcellation==nn;
if nargin>2
    markers=varargin{2};
    mbw=bwconncomp(markers,6);
else
    markers=zeros(size(parcellation),'logical');
    mbw=cell(0);
end

[mgx, mgy, mgz]=meshgrid(1:size(markers,1),1:size(markers,2),1:size(markers,3));

nf=figure('WindowState','maximize');
annotation('textbox','String',['SENPAI Prune' newline 'Mark branches to be pruned.'],...
    'Position',[0.75 0.85 0.2 0.1],'FontWeight','bold','FontSize',15,...
    'HorizontalAlignment','right','EdgeColor','none');
%current index
nnvis=uicontrol('Parent',nf,'Style','text',...
    'String',['Currently visualizing' newline 'neuron ' mat2str(nn) ' out of ' mat2str(N)],...
    'FontWeight','bold','FontSize',12,'Units','normalized',...
    'Position',[0.75 0.4 0.2 0.1],'Visible','on');
%Add pushbutton to put marker
ButtonH=uicontrol('Parent',nf,'Style','pushbutton',...
    'String','Mark branch for pruning','FontWeight','bold','Units','normalized',...
    'Position',[0.75 0.3 0.2 0.1],'Visible','on',...
    'Callback',@senpai_prune_press);
%Add pushbutton to go to next neuron
ButtonNxt=uicontrol('Parent',nf,'Style','pushbutton',...
    'String','>','FontWeight','bold','Units','normalized',...
    'Position',[0.5 0.05 0.05 0.05],'Visible','on',...
    'Callback',@senpai_next);
%Add pushbutton to go to previous neuron
ButtonPrv=uicontrol('Parent',nf,'Style','pushbutton',...
    'String','<','FontWeight','bold','Units','normalized',...
    'Position',[0.1 0.05 0.05 0.05],'Visible','on',...
    'Callback',@senpai_prev);
%neuron index
neuIdx=uicontrol('Parent',nf,'Style','edit',...
    'String',mat2str(nn),'Units','normalized',...
    'Position',[0.3 0.05 0.05 0.05]);
%go to neuron index
ButtonGo=uicontrol('Parent',nf,'Style','pushbutton',...
    'String','Go!','FontWeight','bold','Units','normalized',...
    'Position',[0.35 0.05 0.05 0.05],'Visible','on',...
    'Callback',@senpai_go);
%add pushbutton to save and exit
ButtonEx=uicontrol('Parent',nf,'Style','pushbutton',...
    'String','Save','FontWeight','bold','Units','normalized',...
    'Position',[0.75 0.2 0.2 0.1],'Visible','on',...
    'Callback',@senpai_prune_done);
ax1=subplot(1,4,1:3);

update_fig
axis equal;
box off
xticks([]);yticks([]);zticks([]);
set(gca,'Position',[0 0.1 0.8 0.9])
set(gca,'Visible','off')



    function senpai_prune_press(src,event)
        d = datacursormode(nf);
        vals = getCursorInfo(d);
        [xn, yn, zn]=ind2sub(size(neuron),find(neuron==1));
        coords=vals.Position;
        [~,minc]=min(sum(([yn xn zn]-coords).^2,2));
        %newsph=zeros(size(segmentation),'logical');
        [~, newsph]=regionGrowing(100.*neuron, [xn(minc) yn(minc) zn(minc)],1,5);
        %newsph((mgx-coords(1)).^2+(mgy-coords(2)).^2+(mgz-coords(3)).^2<8)=1;
        markers(newsph)=1;
        [p,v] = isosurface(newsph,0.1);
        subplot(ax1);hold on;patch('Faces',p, 'Vertices',v,'FaceColor','red',...
            'EdgeColor','red','FaceAlpha',0.8);
        disp(['added marker at coordinates ' mat2str(coords(1)) ', ' mat2str(coords(2)) ', ' mat2str(coords(3))])
    end
    function senpai_prune_done(src,event)
        save markers.mat markers
        disp('new markers saved!')
    end
    function senpai_next(src,event)
        nn=min(nn+1,N);
        nnvis.String=['Currently visualizing' newline 'neuron ' mat2str(nn) ' out of ' mat2str(N)];
        neuron=parcellation==nn;
         update_fig
    end
    function senpai_prev(src,event)
        nn=max(nn-1,1);
        nnvis.String=['Currently visualizing' newline 'neuron ' mat2str(nn) ' out of ' mat2str(N)];
        neuron=parcellation==nn;
        update_fig
    end
    function senpai_go(src,event)
        newnn=min(max(str2num(neuIdx.String),1),N);
        neuIdx.String=mat2str(newnn);
        nn=newnn;
        nnvis.String=['Currently visualizing' newline 'neuron ' mat2str(nn) ' out of ' mat2str(N)];
        neuron=parcellation==nn;
        update_fig
    end
    function update_fig
        ff=find(neuron(:)==1);
        if ~isempty(mbw)
            mcheck=find(cellfun(@(x) sum(ismember(ff,x))>0, mbw.PixelIdxList));
        else
            mcheck=[];
        end
        [p,v] = isosurface(neuron,0.5);
        comps=pca(v);
        projs=dot(repmat(comps(:,3)',size(v,1),1),v,2);
        projs=projs-min(projs(:));
        projs=1+projs.*100./max(projs(:));
        clr=tmp(round(projs),:);
        cla(ax1)
        patch(ax1,'Faces',p, 'Vertices',v,'FaceVertexCData',clr,'FaceColor','flat','EdgeColor','none')
        shading interp
        if ~isempty(mcheck)
        mold=zeros(size(neuron),'logical');
        mold(cell2mat(mbw.PixelIdxList(mcheck)'))=1;
        [p1,v1] = isosurface(mold,0.1);
        subplot(ax1);hold on;patch('Faces',p1, 'Vertices',v1,'FaceColor','green',...
            'EdgeColor','green','FaceAlpha',0.8);
        end
        camlight;
        camlight(0, -100);
        lighting phong;
    end
end