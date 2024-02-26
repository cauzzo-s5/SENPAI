function senpai_strahlerord(swc_mat, neuron, path_out, filename)
% senpai_strahlerord:
%     produces a strahler ordering and a series of related statistics.
%     the output is stored in an output .mat file
%
%     Execute the function in the command window:
%     Syntax:
%
%       senpai_strahlerord(swc_mat, path_out, filename):
%
%       inputs
%
%       swc = swc-like matrix encoding the skeleton, as produced
%       by senpai_skeletonize.m
%       path_out = path where to store output .m file
%       filename = name of output .m file
%
%       output statistics stored in the output .mat file:
%       numSegSO (number of segments per SO)
%       numSegSOnorm (number of segments per SO, normalized by the total number of segments)
%       segLAve (average segment length per SO)
%       PSnum (coefficients of the linear fit for the normalized number of segments per SO)
%       segDAve (average segment diameter per SO) 
%       TOTL (total length of segments)
%       TopoSubLAve (topological subtree size per SO) 
%       numBrSO (number of branches per SO)
%       numBrSOnorm (number of branches per SO, normalized by the total number of branches)
%       brDAve (average branch diameter per SO)
%       brLAve (average branch length per SO)
%       PBnum (coefficients of the linear fit for the normalized number of branches per SO) 
%       normTotL (total dendritic length per SO)


%read from swc
xn=round(swc_mat(:,3));yn=round(swc_mat(:,4));zn=round(swc_mat(:,5));zn(zn<1)=1;xn(xn<1)=1;
diam=swc_mat(:,6);root_cand=find(swc_mat(:,7)==-1);
bifurc=find(swc_mat(:,2)==5);
leaves=find(swc_mat(:,2)==6);
pred=swc_mat(:,7);
ni=pred;
ne=swc_mat(:,1);
t=graph(sparse([ni(ni>0);ne(ni>0)],[ne(ni>0);ni(ni>0)],[diam(ni>0);diam(ni>0)],length(ni),length(ni),2*length(ni)));
t=minspantree(t,'Type','tree','Root',root_cand);

%produce strahler ordering
% copy tree. assign new weights:
% 0.5 to all edges in contact with a bifurcation
% 0 to all the others
% compute minimum path length between soma and all the other nodes in
% the tree
% assign to each node the weigth of the minimum path length, rounded to
% the floor
t_n=t;
t_n.Edges.Weight(:)=0;
t_n.Edges.Weight(union(find(ismember(t_n.Edges.EndNodes(:,1),bifurc)),find(ismember(t_n.Edges.EndNodes(:,2),bifurc))))=0.5;
[~, D]= shortestpathtree(t_n,find(pred==-1),'Method','positive');
D=floor(D);
D(D==inf)=max(D(isfinite(D)));
n1col=zeros(size(neuron));n1col(sub2ind(size(neuron),xn,yn,zn))=D;

% find branches length (schematize tree first)
%each edge connected to a root, leaf or bifurcation has weight 0.5, each
%other edge has weight zero. Minpathlength gives the number of segments
%that have been crossed.
t_d=t_n;
t_d.Edges.Weight(union(find(ismember(t_n.Edges.EndNodes(:,1),leaves)),find(ismember(t_n.Edges.EndNodes(:,2),leaves))))=0.5;
t_d.Edges.Weight(union(find(ismember(t_n.Edges.EndNodes(:,1),root_cand)),find(ismember(t_n.Edges.EndNodes(:,2),root_cand))))=0.5;
branches=cell(0);
towr=1;
for bn=[root_cand union(leaves',bifurc')]
    targets=unique([root_cand union(leaves',bifurc')]);
    targets(targets==bn)=[];
    [TRbif, Dbif]= shortestpathtree(t_d,bn,targets,'OutputForm','cell');
    % trova i Dbif==1
    for rr=find(Dbif<=1)
        %list of nodes
        branches{towr,1}=TRbif{rr};
        %order
        branches{towr,2}=mode(n1col(sub2ind(size(neuron),xn(TRbif{rr}),yn(TRbif{rr}),zn(TRbif{rr}))));
        %length
        tmpd=[xn(TRbif{rr}),yn(TRbif{rr}),zn(TRbif{rr})];
        branches{towr,3}=sum(sqrt(sum(diff(tmpd).^2)));
        %diameter
        branches{towr,5}=mean(diam(TRbif{rr}));
        towr=towr+1;
    end
    %for the node under consideration, set a high weight for all adiacent edges
    t_d.Edges.Weight(union(find(ismember(t_d.Edges.EndNodes(:,1),bn)),find(ismember(t_d.Edges.EndNodes(:,2),bn))))=10;
end
% "branches" columns: 1) list of nodes; 2) strahler order; 3) length; 
% 4) []; 5) diameter; 6) []; 7) [] 
%% strahler
branches(:,4)={0};
branches(:,6)={-1};
branches(cellfun(@(x) sum(ismember(x,leaves))>0,branches(:,1)),[4 7])={1};
tosee=find(cell2mat(branches(:,4))==0);
actual=2;
prev_lts=length(tosee);
fprintf('%g segments left...\n',length(tosee))
while ~isempty(tosee)
    for nn=tosee'
        if length(tosee)~=prev_lts
            fprintf('%g segments left...\n',length(tosee))
            prev_lts=length(tosee);
        end
        %find unseen connected segments
        apex=branches{nn,1};
        newbranches1=find(cellfun(@(x) sum(ismember(x,apex(1)))>0,branches(:,1)));
        newbranches_val1=cell2mat(branches(newbranches1,4));
        newbranches1=newbranches1(newbranches_val1>0);
        newbranches_val1=newbranches_val1(newbranches_val1>0);
        newbranches2=find(cellfun(@(x) sum(ismember(x,apex(end)))>0,branches(:,1)));
        newbranches_val2=cell2mat(branches(newbranches2,4));
        newbranches2=newbranches2(newbranches_val2>0);
        newbranches_val2=newbranches_val2(newbranches_val2>0);
        if length(newbranches_val1)>1 %if on one side there are two segments already seen
            if newbranches_val1(1)==newbranches_val1(2)
                branches{nn,4}=newbranches_val1(1)+1;
                branches{nn,6}=-1;
            else
                [branches{nn,4}, hh]=max(newbranches_val1);
                branches{nn,6}=newbranches1(hh);
            end
            %add subtree_length
            branches{nn,7}=branches{newbranches1(1),7}+branches{newbranches1(2),7}+1;
        elseif length(newbranches_val2)>1 %if on the other side there are two segments already seen
            if newbranches_val2(1)==newbranches_val2(2)
                branches{nn,4}=newbranches_val2(1)+1;
                branches{nn,6}=-1;
            else
                [branches{nn,4}, hh]=max(newbranches_val2);
                branches{nn,6}=newbranches2(hh);
            end
            branches{nn,7}=branches{newbranches2(1),7}+branches{newbranches2(2),7}+1;
        end
    end
    tosee=find(cell2mat(branches(:,4))==0);
    actual=actual+1;
end

%% statistics
ords=cell2mat(branches(:,4));
ccl=parula(max(ords));
figure;
a=subplot(2,2,1);
[p,v]=isosurface(neuron,0.7);
patch('Faces',p,'Vertices',v,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none','FaceAlpha',0.5);
for ccv=1:length(branches)
    hold on
    ccb=branches{ccv,1};
    plot3(xn(ccb),yn(ccb),zn(ccb),'*','MarkerEdgeColor',ccl(branches{ccv,4},:))
end
axis equal;axis off;set(a,'Position',[0 0.5 0.5 0.5])
hold on;plot3(xn(root_cand),yn(root_cand),zn(root_cand),'*','MarkerEdgeColor','blue','MarkerSize',15,'Linewidth',10)

% summary stats to save
TOTL=sum(cell2mat(branches(:,3)));
numSegSO=histcounts(ords,max(ords));
numSegSOnorm=numSegSO./sum(numSegSO);
numNod=cellfun(@length,branches(:,1));
numNodSO=arrayfun(@(x) sum(numNod(ords==x)),1:length(numSegSO));
segLAve=arrayfun(@(x) sum(cell2mat(branches(ords==x,3))),1:length(numSegSO))./numSegSO;
segDtmp=cell2mat(branches(:,5)).*numNod;
segDAve=arrayfun(@(x) sum(segDtmp(ords==x)),1:length(numSegSO))./numNodSO;
TopoSubLAve=arrayfun(@(x) mean(cell2mat(branches(ords==x,7))./(length(branches))),2:length(numSegSO));
% populations to plot
subplot(2,2,2);boxplot(cell2mat(branches(:,3)),ords,'whisker',0);title('segment length');xlabel('SO')
subplot(2,2,3);boxplot(cell2mat(branches(:,5)),ords,'whisker',0);title('segment diameter');xlabel('SO')
subplot(2,2,4);plot(1:length(numSegSO),numSegSOnorm);title('#segment (norm)');xlabel('SO');
%fit
PSnum=polyfit(log2(1:length(numSegSO)),log2(numSegSO./sum(numSegSO)),1); %log2 or normalization does not change slope, only intercept does.
hold on;plot(1:length(numSegSO),2.^(-(1:length(numSegSO))));
set(gca, 'YScale', 'log')
title(['normalized segment count (k=' mat2str(-PSnum(1),2) ,')'])
save([path_out filename '.mat'],'numSegSO','numSegSOnorm','segLAve','PSnum','segDAve','TOTL','TopoSubLAve');
%% branches (no segments)
toseevec=zeros(size(branches,1),1);
deps=cell2mat(branches(:,6));
count=1;
for bb=max(ords):-1:1
    tmpSO=find(ords==bb & toseevec==0);
    for bbb=tmpSO'
        if ismember(bbb,deps);continue;end
        tmpNodes=[];
        tmpseg=bbb;
        while tmpseg>0
            tmpNodes=[tmpNodes branches{tmpseg,1}];
            toseevec(tmpseg)=1;
            tmpseg=branches{tmpseg,6};
        end
        BB{count,1}=tmpNodes;
        BB{count,2}=bb;

        %length
        tmpd=[xn(tmpNodes),yn(tmpNodes),zn(tmpNodes)];
        BB{count,3}=sum(sqrt(sum(diff(tmpd).^2)));
        %diameter
        BB{count,4}=mean(diam(tmpNodes));
        count=count+1;
    end
end
ordsBB=cell2mat(BB(:,2));
ccl=hot(max(ordsBB));
figure
numBrSO=histcounts(ordsBB,max(ordsBB));
numBrSOnorm=numBrSO./sum(numBrSO);
numNodBB=cellfun(@length,BB(:,1));
brDtmp=cell2mat(BB(:,4)).*numNodBB;
brDAve=arrayfun(@(x) sum(brDtmp(ordsBB==x)),1:length(numBrSO))./numNodSO;
brDAve=rescale(brDAve);
brLAve=arrayfun(@(x) sum(cell2mat(BB(ordsBB==x,3))),1:length(numBrSO))./numBrSO;
subplot(2,2,1);boxplot(cell2mat(BB(:,3)),ordsBB,'whisker',0);title('branches length');xlabel('SO')
subplot(2,2,2);boxplot(cell2mat(BB(:,4)),ordsBB,'whisker',0);title('branches norm diameter');xlabel('SO')
subplot(2,2,3);plot(1:length(numBrSO),numBrSOnorm);title('#branch (norm)');xlabel('SO');
PBnum=polyfit(log2(1:length(numBrSO)),log2(numBrSO./sum(numBrSO)),1); %log2 or normalization does not change slope, only intercept does.
hold on;plot(1:length(numBrSO),4.^(1-(1:length(numBrSO))));
set(gca, 'YScale', 'log')
title(['normalized segment count (k=' mat2str(-PBnum(1),2) ,')'])
%total dendritic length
normTotL=arrayfun(@(x) sum(cell2mat(BB(ordsBB==x,3))),1:length(numBrSO))./(sum(cell2mat(BB(:,3))));
subplot(2,2,4);plot(1:length(numBrSO),normTotL);title('total dendritic length');xlabel('SO');
set(gca, 'YScale', 'log')
save([path_out filename '.mat'],'numBrSO','numBrSOnorm','brDAve','brLAve','PBnum','normTotL','-append');

end