function [t,swc]=senpai_skeletonize(cIM,neuron,somas)
    % senpai_skeletonize:
    %     produces a matlab tree structure and an swc-formatted matrix for
    %     the skeleton of a binary segmentation
    %
    %     Execute the function in the command window:
    %     Syntax:
    %
    %       [t,pred,swc,xn,yn,zn]=senpai_skeletonize(cIM,neuron,somas):
    %
    %       inputs
    %
    %       cIM = original image stack on which the segmetation has been
    %       produced
    %       neuron = binary segmentation of a single neuron
    %       somas = binary mask of one (or multiple) somas
    %
    %       outputs
    %       
    %       t = minimum spanning tree (matlab structure) of the neuron
    %       swc = swc-like matrix encoding the skeleton
    %

somas=somas>0;

%include soma mask in the segmentation
neuron=neuron>0 | somas;
%keep only the correct neuron+soma
neusel=bwconncomp(neuron,26);
[~,neuselN]=max(cellfun(@length,neusel.PixelIdxList));
neuron=neuron.*0;
neuron(neusel.PixelIdxList{neuselN})=1;
%fill holes
neuron=imfill(neuron==1,'holes');
somas=somas.*neuron;
%distance transform
bwd=bwdist(~neuron);
%start from simple matlab skeletonization
sk=bwskel(neuron,'MinBranchLength', 3);
%clean skeleton in soma
%sk(somas==1)=0;
% points = indeces of voxels in the skeleton
points=find(sk);
% xn yn zn: coordinates of the voxels in the skeleton
[xn, yn, zn]=ind2sub(size(neuron),points);

%this loop converts the matlab skeleton to a graph structure
done=zeros(1,length(points));
ni=[]; %edge origin
ne=[]; %edge end
s=[]; %weight
count=1;
while sum(done)<length(points) % fin che non li ho visti tutti
    fprintf('building skeleton: %g nodes left...\n',length(points)-sum(done));
    % define seed
    seedpoint=find(~done,1);
    seed=[xn(seedpoint) yn(seedpoint) zn(seedpoint)];
    % find neighbor
    nextmatch=find(sum(([xn yn zn]-seed).^2,2)<=3);
    nextmatch(nextmatch==seedpoint)=[];
    nextmatch(ismember(nextmatch,find(done)))=[];
    % define connection
    for vv=1:length(nextmatch)
        coupl=sort([seedpoint nextmatch(vv)]);
        ni(count)=coupl(1);
        ne(count)=coupl(2);
        s(count)=cIM(points(nextmatch(vv)));
        count=count+1;
    end
    done(seedpoint)=1;
end
%build graph from vectors ni, ne, s
graphtree=graph(sparse(ni,ne,s,length(points),length(points),length(ni)),'upper');

%convert graph to tree: check for cycles and remove them
G_ac=graphtree;
[~,edgecycles] = allcycles(G_ac,'MaxNumCycles',1);
while ~isempty(edgecycles)
    %length(edgecycles)
    cc=1;
    edgetmp=edgecycles{cc};
    % cost function for cycle cut: width+intensity
    intdiff=min(bwd(points(G_ac.Edges.EndNodes(:,2))),bwd(points(G_ac.Edges.EndNodes(:,1))))+min(cIM(points(G_ac.Edges.EndNodes(:,2))),cIM(points(G_ac.Edges.EndNodes(:,1))))./255;
    [~, mididx]=min(intdiff(edgetmp));
    G_ac=rmedge(G_ac,edgetmp(mididx));
    [~,edgecycles] = allcycles(G_ac,'MaxNumCycles',1);
end

%check if graph is connected (single apical dendrite) or disconnected
%(multiple dendrites departing from the soma). In the second case, connect
%all the dendrites to the soma
if length(conncomp(G_ac))>1
    %in this case I have multiple disconnected graphs and I have to
    %produce a single one.
    %find leaves that are close to soma.
    xyzn=sub2ind(size(cIM),xn,yn,zn);
    roots_lists_1=G_ac.Edges.EndNodes((somas(xyzn(G_ac.Edges.EndNodes(:,1)))-somas(xyzn(G_ac.Edges.EndNodes(:,2))))~=0,:);
    mask_roots=somas(xyzn(G_ac.Edges.EndNodes((somas(xyzn(G_ac.Edges.EndNodes(:,1)))-somas(xyzn(G_ac.Edges.EndNodes(:,2))))~=0,:)));
    roots_lists=xyzn(roots_lists_1(~mask_roots(:)));
    roots_lists_opp=xyzn(roots_lists_1(mask_roots(:)==1));
    somacenter=regionprops3(somas==1,'centroid');
    xn=[xn;round(somacenter.Centroid(2))];yn=[yn;round(somacenter.Centroid(1))];zn=[zn;round(somacenter.Centroid(3))];
    points=[points; sub2ind(size(neuron),round(somacenter.Centroid(1,1)),round(somacenter.Centroid(1,2)),round(somacenter.Centroid(1,3)))];
    for rr=find(ismember(sub2ind(size(neuron),xn,yn,zn),roots_lists))'
        G_ac=addedge(G_ac,length(xn),rr,cIM(sub2ind(size(cIM),xn(end),yn(end),zn(end))));
    end
end

%define soma node
xyzn=sub2ind(size(cIM),xn,yn,zn);
G_ac_new=rmedge(G_ac,find(somas(xyzn(G_ac.Edges.EndNodes(:,1)))==1 & somas(xyzn(G_ac.Edges.EndNodes(:,2)))==1));
root_cand=length(xn);
[t, pred]=minspantree(G_ac_new,'Type','tree','Root',root_cand);
tmpr=zeros(size(pred));tmpr(root_cand)=1;
t=rmnode(t,find(isnan(pred)));
tmpr(isnan(pred))=[];
root_cand=find(tmpr);
xn(isnan(pred))=[];
yn(isnan(pred))=[];
zn(isnan(pred))=[];
points(isnan(pred))=[];
[t, pred]=minspantree(t,'Type','tree','Root',root_cand);
%define leaf and bifurcation nodes
% nodi che appaiono 1 volta nella lista degli edge sono leaf, 3 volte sono
% biforcazioni. 2 volte è normale
[gc,gr]=groupcounts(t.Edges.EndNodes(:));
bifurc=gr(gc>2);
leaves=gr(gc==1);
leaves(leaves==root_cand)=[];

%build swc-like matrix
pred(pred==0)=-1;
swc=zeros(length(pred),7);
swc(:,1)=1:length(pred);
swc(:,2)=0;
swc(pred==0,2)=-1;
swc(leaves,2)=6;swc(bifurc,2)=5;
swc(:,3)=yn;
swc(:,4)=xn;
swc(:,5)=zn;
swc(:,6)=bwd(points);
swc(:,7)=pred;

% to produce a .swc file use:

% writematrix(swc,'neuron_skel.txt','Delimiter',' ');
% movefile neuron_skel.txt neuron_skel.swc

end