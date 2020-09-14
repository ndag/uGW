addpath(genpath('uGW-Project'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% uGW phylo tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('uSLBInfluenzaTrees.mat')
L=length(uSLB)

min1=0;
max1=ceil(max(uSLB(:)));
m= linkage(squareform(uSLB(1:(L/2),1:(L/2))));
m2= linkage(squareform(uSLB((L/2+1):L,(L/2+1):L)));

test = m(160:200);

[H,T,outperm] =dendrogram(m,(L/2))
[H,T,outperm2] =dendrogram(m2,(L/2))
outperm2 = outperm2+200
uSLB(1:(L/2),:)=uSLB(outperm,:)
uSLB(:,1:(L/2))=uSLB(:,outperm)

uSLB((L/2+1):L,:)=uSLB(outperm2,:)
uSLB(:,(L/2+1):L)=uSLB(:,outperm2)

fig=figure; 
caxis manual
caxis([min1,max1]);
imagesc(uSLB,[min1,max1])
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colorbar;
axis square
saveas(fig,'uLSBphylo.pdf')
!pdfcrop uLSBphylo.pdf uLSBphylo.pdf

% Single Linkage dendrogram
cols = ones(400,1);
cols(201:end) = 2 ; 

fig=figure; 
msl= linkage(squareform(uSLB),'single');
h=dendrogram(msl,0,'Labels',num2str(cols))
set(h,'LineWidth',0.5)
set(h, 'Color', 'black');
ax = gca; % get the axes handle
lab = ax.XAxis.TickLabels; % get all the labels
loc = ax.XAxis.TickValues; % get all labels location
[ulab,~,lab_ind] = unique(lab); % find all unique labels
clr = lines(numel(ulab)); % make a color map
X = get(ax.Children,'XData'); % Get x values of all lines
Y = get(ax.Children,'YData'); % Get y values of all lines

for k = 1:numel(X)
    if X{k}(1)==fix(X{k}(1))
        line(ax,X{k}(1:2),Y{k}(1:2),'Color',clr(lab_ind(X{k}(1)),:),...
            'LineWidth',0.5);
    end
    if X{k}(2)==fix(X{k}(2))
       line(ax,X{k}(2:3),Y{k}(2:3),'Color',clr(lab_ind(X{k}(2)),:),...
       'LineWidth',0.5);
    end
    if X{k}(3)==fix(X{k}(3))
        line(ax,X{k}(3:4),Y{k}(3:4),'Color',clr(lab_ind(X{k}(3)),:),...
            'LineWidth',0.5);
    end
end
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gca,'XTick',[])
saveas(fig,'uLSBsingle.pdf')
!pdfcrop uLSBsingle.pdf uLSBsingle.pdf
% Complete Linkage dendrogram

mcl= linkage(squareform(uSLB),'complete');
fig=figure;
h=dendrogram(mcl,0,'Labels',num2str(cols))
set(h,'LineWidth',0.5)
set(h, 'Color', 'black');
ax = gca; % get the axes handle
lab = ax.XAxis.TickLabels; % get all the labels
loc = ax.XAxis.TickValues; % get all labels location
[ulab,~,lab_ind] = unique(lab); % find all unique labels
clr = lines(numel(ulab)); % make a color map
X = get(ax.Children,'XData'); % Get x values of all lines
Y = get(ax.Children,'YData'); % Get y values of all lines

for k = 1:numel(X)
    if X{k}(1)==fix(X{k}(1))
        line(ax,X{k}(1:2),Y{k}(1:2),'Color',clr(lab_ind(X{k}(1)),:),...
            'LineWidth',0.5);
    end
    if X{k}(2)==fix(X{k}(2))
       line(ax,X{k}(2:3),Y{k}(2:3),'Color',clr(lab_ind(X{k}(2)),:),...
       'LineWidth',0.5);
    end
    if X{k}(3)==fix(X{k}(3))
        line(ax,X{k}(3:4),Y{k}(3:4),'Color',clr(lab_ind(X{k}(3)),:),...
            'LineWidth',0.5);
    end
end
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gca,'XTick',[])
saveas(fig,'uLSBcomplete.pdf')
!pdfcrop uLSBcomplete.pdf uLSBcomplete.pdf
%mds
Y = cmdscale(uSLB,2)
fig = figure;
plot(Y(:,1),Y(:,2),'ro')
hold on
plot(Y(1:(L/2),1),Y(1:(L/2),2), 'blueo')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
saveas(fig,'uLSBphylomds.pdf')
!pdfcrop uLSBphylomds.pdf uLSBphylomds.pdf


fig=figure;
Y = cmdscale(uSLB,3)
fig = figure;
plot3(Y(:,1),Y(:,2),Y(:,3),'ro')
hold on
plot3(Y(1:(L/2),1),Y(1:(L/2),2),Y(1:(L/2),3), 'blueo')
direction = [0 1 0];
set(findall(gcf,'-property','FontSize'),'FontSize',12)
saveas(fig,'uLSBmds3d.pdf')
!pdfcrop uLSBmds3d.pdf uLSBmds3d.pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dGW phylo tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('SLBInfluenzaTrees.mat')
L=length(SLB)
cols = ones(400,1);
cols(201:end) = 2 ; 
%min1=0;
%max1=max(SLB(:));
m= linkage(squareform(SLB(1:(L/2),1:(L/2))));
m2= linkage(squareform(SLB((L/2+1):L,(L/2+1):L)));

test2 = m(160:200);


[H,T,outperm] =dendrogram(m,(L/2))
[H,T,outperm2] =dendrogram(m2,(L/2))
outperm2 = outperm2+200
SLB(1:(L/2),:)=SLB(outperm,:)
SLB(:,1:(L/2))=SLB(:,outperm)

SLB((L/2+1):L,:)=SLB(outperm2,:)
SLB(:,(L/2+1):L)=SLB(:,outperm2)

fig=figure; 
caxis manual
caxis([min1,max1]);
imagesc(SLB,[min1,max1])
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colorbar;
axis square
saveas(fig,'SLBphylo.pdf')
!pdfcrop SLBphylo.pdf SLBphylo.pdf



% Single Linkage dendrogram
fig=figure;
msl= linkage(squareform(SLB),'single');
h=dendrogram(msl,0,'Labels',num2str(cols))
set(h,'LineWidth',0.5)
set(h, 'Color', 'black');
ax = gca; % get the axes handle
lab = ax.XAxis.TickLabels; % get all the labels
loc = ax.XAxis.TickValues; % get all labels location
[ulab,~,lab_ind] = unique(lab); % find all unique labels
clr = lines(numel(ulab)); % make a color map
X = get(ax.Children,'XData'); % Get x values of all lines
Y = get(ax.Children,'YData'); % Get y values of all lines

for k = 1:numel(X)
    if X{k}(1)==fix(X{k}(1))
        line(ax,X{k}(1:2),Y{k}(1:2),'Color',clr(lab_ind(X{k}(1)),:),...
            'LineWidth',0.5);
    end
    if X{k}(2)==fix(X{k}(2))
       line(ax,X{k}(2:3),Y{k}(2:3),'Color',clr(lab_ind(X{k}(2)),:),...
       'LineWidth',0.5);
    end
    if X{k}(3)==fix(X{k}(3))
        line(ax,X{k}(3:4),Y{k}(3:4),'Color',clr(lab_ind(X{k}(3)),:),...
            'LineWidth',0.5);
    end
end
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gca,'XTick',[])
saveas(fig,'SLBsingle.pdf')
!pdfcrop SLBsingle.pdf SLBsingle.pdf

% Complete Linkage dendrogram
fig=figure;
mcl= linkage(squareform(SLB),'complete');
h=dendrogram(mcl,0,'Labels',num2str(cols))
set(h,'LineWidth',0.5)
set(h, 'Color', 'black');
ax = gca; % get the axes handle
lab = ax.XAxis.TickLabels; % get all the labels
loc = ax.XAxis.TickValues; % get all labels location
[ulab,~,lab_ind] = unique(lab); % find all unique labels
clr = lines(numel(ulab)); % make a color map
X = get(ax.Children,'XData'); % Get x values of all lines
Y = get(ax.Children,'YData'); % Get y values of all lines

for k = 1:numel(X)
    if X{k}(1)==fix(X{k}(1))
        line(ax,X{k}(1:2),Y{k}(1:2),'Color',clr(lab_ind(X{k}(1)),:),...
            'LineWidth',0.5);
    end
    if X{k}(2)==fix(X{k}(2))
       line(ax,X{k}(2:3),Y{k}(2:3),'Color',clr(lab_ind(X{k}(2)),:),...
       'LineWidth',0.5);
    end
    if X{k}(3)==fix(X{k}(3))
        line(ax,X{k}(3:4),Y{k}(3:4),'Color',clr(lab_ind(X{k}(3)),:),...
            'LineWidth',0.5);
    end
end
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gca,'XTick',[])
saveas(fig,'SLBcomplete.pdf')
!pdfcrop SLBcomplete.pdf SLBcomplete.pdf

%mds
Y = cmdscale(SLB,2)
fig = figure;
plot(Y(:,1),Y(:,2),'ro')
hold on
plot(Y(1:(L/2),1),Y(1:(L/2),2), 'blueo')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
saveas(fig,'SLBphylomds.pdf')
!pdfcrop SLBphylomds.pdf SLBphylomds.pdf

Y = cmdscale(SLB,3)
fig = figure;
plot3(Y(:,1),Y(:,2),Y(:,3),'ro')
hold on
plot3(Y(1:(L/2),1),Y(1:(L/2),2),Y(1:(L/2),3), 'blueo')
set(findall(gcf,'-property','FontSize'),'FontSize',12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Colijn Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% distmat plot
load('colijn_distmat.mat')
cols = ones(400,1);
cols(201:end) = 2 ; 
%max1=max(SLB(:));
m= linkage(squareform(colijn_distmat(1:(L/2),1:(L/2))));
m2= linkage(squareform(colijn_distmat((L/2+1):L,(L/2+1):L)));



[H,T,outperm] =dendrogram(m,(L/2));
[H,T,outperm2] =dendrogram(m2,(L/2));
outperm2 = outperm2+200
colijn_distmat(1:(L/2),:)=colijn_distmat(outperm,:)
colijn_distmat(:,1:(L/2))=colijn_distmat(:,outperm)

colijn_distmat((L/2+1):L,:)=colijn_distmat(outperm2,:)
colijn_distmat(:,(L/2+1):L)=colijn_distmat(:,outperm2)


fig=figure; 
caxis manual
caxis([min1,max1]);
imagesc(colijn_distmat)
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colorbar;
axis square
saveas(fig,'colijn.pdf')
!pdfcrop colijn.pdf colijn.pdf


fig=figure; 
caxis manual
caxis([min1,max1]);
imagesc(colijn_distmat)
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colorbar;
axis square
% Single Linkage dendrogram
    
msl= linkage(squareform(colijn_distmat),'single');
fig=figure; 
h=dendrogram(msl,0,'Labels',num2str(cols))
set(h,'LineWidth',0.5)
set(h, 'Color', 'black');
ax = gca; % get the axes handle
lab = ax.XAxis.TickLabels; % get all the labels
loc = ax.XAxis.TickValues; % get all labels location
[ulab,~,lab_ind] = unique(lab); % find all unique labels
clr = lines(numel(ulab)); % make a color map
X = get(ax.Children,'XData'); % Get x values of all lines
Y = get(ax.Children,'YData'); % Get y values of all lines

for k = 1:numel(X)
    if X{k}(1)==fix(X{k}(1))
        line(ax,X{k}(1:2),Y{k}(1:2),'Color',clr(lab_ind(X{k}(1)),:),...
            'LineWidth',0.5);
    end
    if X{k}(2)==fix(X{k}(2))
       line(ax,X{k}(2:3),Y{k}(2:3),'Color',clr(lab_ind(X{k}(2)),:),...
       'LineWidth',0.5);
    end
    if X{k}(3)==fix(X{k}(3))
        line(ax,X{k}(3:4),Y{k}(3:4),'Color',clr(lab_ind(X{k}(3)),:),...
            'LineWidth',0.5);
    end
end
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gca,'XTick',[])
saveas(fig,'colijnsingle.pdf')
!pdfcrop colijnsingle.pdf colijnsingle.pdf
% Complete Linkage dendrogram

mcl= linkage(squareform(colijn_distmat),'complete');
fig=figure; 
h=dendrogram(mcl,0,'Labels',num2str(cols))
set(h,'LineWidth',0.5)
set(h, 'Color', 'black');
ax = gca; % get the axes handle
lab = ax.XAxis.TickLabels; % get all the labels
loc = ax.XAxis.TickValues; % get all labels location
[ulab,~,lab_ind] = unique(lab); % find all unique labels
clr = lines(numel(ulab)); % make a color map
X = get(ax.Children,'XData'); % Get x values of all lines
Y = get(ax.Children,'YData'); % Get y values of all lines

for k = 1:numel(X)
    if X{k}(1)==fix(X{k}(1))
        line(ax,X{k}(1:2),Y{k}(1:2),'Color',clr(lab_ind(X{k}(1)),:),...
            'LineWidth',0.5);
    end
    if X{k}(2)==fix(X{k}(2))
       line(ax,X{k}(2:3),Y{k}(2:3),'Color',clr(lab_ind(X{k}(2)),:),...
       'LineWidth',0.5);
    end
    if X{k}(3)==fix(X{k}(3))
        line(ax,X{k}(3:4),Y{k}(3:4),'Color',clr(lab_ind(X{k}(3)),:),...
            'LineWidth',0.5);
    end
end
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gca,'XTick',[])
saveas(fig,'colijncomplete.pdf')
!pdfcrop colijncomplete.pdf colijncomplete.pdf
% mds

Y = cmdscale(colijn_distmat,2)
fig = figure;
plot(Y(:,1),Y(:,2),'ro')
hold on
plot(Y(1:(L/2),1),Y(1:(L/2),2), 'blueo')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
saveas(fig,'colijnmds.pdf')
!pdfcrop colijnmds.pdf colijnmds.pdf

Y = cmdscale(colijn_distmat,3)
fig = figure;
plot3(Y(:,1),Y(:,2),Y(:,3),'ro')
hold on
plot3(Y(1:(L/2),1),Y(1:(L/2),2),Y(1:(L/2),3), 'blueo')
set(findall(gcf,'-property','FontSize'),'FontSize',12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% uGW weighted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('uSLBInfluenzaTreesweighted.mat')
L=length(uSLB)

min1=0;
max1=(max(uSLB(:)));
m= linkage(squareform(uSLB(1:(L/2),1:(L/2))));
m2= linkage(squareform(uSLB((L/2+1):L,(L/2+1):L)));

test3 = m(160:200);


[H,T,outperm] =dendrogram(m,(L/2))
[H,T,outperm2] =dendrogram(m2,(L/2))
outperm2 = outperm2+200
uSLB(1:(L/2),:)=uSLB(outperm,:)
uSLB(:,1:(L/2))=uSLB(:,outperm)

uSLB((L/2+1):L,:)=uSLB(outperm2,:)
uSLB(:,(L/2+1):L)=uSLB(:,outperm2)

fig=figure; 
caxis manual
caxis([min1,max1]);
imagesc(uSLB,[min1,max1])
set(findall(gcf,'-property','FontSize'),'FontSize',12)
colorbar;
axis square
saveas(fig,'uLSBphyloWeighted.pdf')
!pdfcrop uLSBphyloWeighted.pdf uLSBphyloWeighted.pdf

% Single Linkage dendrogram
cols = ones(400,1);
cols(201:end) = 2 ; 

fig=figure; 
msl= linkage(squareform(uSLB),'single');
h=dendrogram(msl,0,'Labels',num2str(cols))
set(h,'LineWidth',0.5)
set(h, 'Color', 'black');
ax = gca; % get the axes handle
lab = ax.XAxis.TickLabels; % get all the labels
loc = ax.XAxis.TickValues; % get all labels location
[ulab,~,lab_ind] = unique(lab); % find all unique labels
clr = lines(numel(ulab)); % make a color map
X = get(ax.Children,'XData'); % Get x values of all lines
Y = get(ax.Children,'YData'); % Get y values of all lines

for k = 1:numel(X)
    if X{k}(1)==fix(X{k}(1))
        line(ax,X{k}(1:2),Y{k}(1:2),'Color',clr(lab_ind(X{k}(1)),:),...
            'LineWidth',0.5);
    end
    if X{k}(2)==fix(X{k}(2))
       line(ax,X{k}(2:3),Y{k}(2:3),'Color',clr(lab_ind(X{k}(2)),:),...
       'LineWidth',0.5);
    end
    if X{k}(3)==fix(X{k}(3))
        line(ax,X{k}(3:4),Y{k}(3:4),'Color',clr(lab_ind(X{k}(3)),:),...
            'LineWidth',0.5);
    end
end
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gca,'XTick',[])
saveas(fig,'uLSBsingleWeighted.pdf')
!pdfcrop uLSBsingleWeighted.pdf uLSBsingleWeighted.pdf
% Complete Linkage dendrogram

mcl= linkage(squareform(uSLB),'complete');
fig=figure;
h=dendrogram(mcl,0,'Labels',num2str(cols))
set(h,'LineWidth',0.5)
set(h, 'Color', 'black');
ax = gca; % get the axes handle
lab = ax.XAxis.TickLabels; % get all the labels
loc = ax.XAxis.TickValues; % get all labels location
[ulab,~,lab_ind] = unique(lab); % find all unique labels
clr = lines(numel(ulab)); % make a color map
X = get(ax.Children,'XData'); % Get x values of all lines
Y = get(ax.Children,'YData'); % Get y values of all lines

for k = 1:numel(X)
    if X{k}(1)==fix(X{k}(1))
        line(ax,X{k}(1:2),Y{k}(1:2),'Color',clr(lab_ind(X{k}(1)),:),...
            'LineWidth',0.5);
    end
    if X{k}(2)==fix(X{k}(2))
       line(ax,X{k}(2:3),Y{k}(2:3),'Color',clr(lab_ind(X{k}(2)),:),...
       'LineWidth',0.5);
    end
    if X{k}(3)==fix(X{k}(3))
        line(ax,X{k}(3:4),Y{k}(3:4),'Color',clr(lab_ind(X{k}(3)),:),...
            'LineWidth',0.5);
    end
end
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gca,'XTick',[])
saveas(fig,'uLSBcompleteWeighted.pdf')
!pdfcrop uLSBcompleteWeighted.pdf uLSBcompleteWeighted.pdf
%mds
Y = cmdscale(uSLB,2)
fig = figure;
plot(Y(:,1),Y(:,2),'ro')
hold on
plot(Y(1:(L/2),1),Y(1:(L/2),2), 'blueo')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
saveas(fig,'uLSBphylomdsWeighted.pdf')
!pdfcrop uLSBphylomdsWeighted.pdf uLSBphylomdsWeighted.pdf


fig=figure;
Y = cmdscale(uSLB,3)
fig = figure;
plot3(Y(:,1),Y(:,2),Y(:,3),'ro')
hold on
plot3(Y(1:(L/2),1),Y(1:(L/2),2),Y(1:(L/2),3), 'blueo')
direction = [0 1 0];
set(findall(gcf,'-property','FontSize'),'FontSize',12)
saveas(fig,'uLSBmds3dWeighted.pdf')
!pdfcrop uLSBmds3dWeighted.pdf uLSBmds3dWeighted.pdf

