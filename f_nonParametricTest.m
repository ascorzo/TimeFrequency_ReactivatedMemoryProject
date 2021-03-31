function stats = f_nonParametricTest(ROI,xMeassure,yMeassure,condition1,condition2,data1,data2,Xt,color,start_point,end_sample)

c1 = data1; % Vector de Condicion 1
c2 = data2; %Vector de Condicion 2

data = {data1'; data2'};
paired = 'on';

method = 'bootstrap';
perms = 10000;

[stats, ~, pvals, ~] = statcond( data,'paired', paired,...
    'method',method,'naccu',perms);

[~, pmask] = fdr(pvals(:),0.05,'nonParametric');
%pmask = pvals;
cluster = 'on';

grupo = ROI;
yaxis = [nanmean(nanmean([c1;c2]))-7*nanstd(nanstd([c1;c2])), nanmean(nanmean([c1;c2]))+7*nanstd(nanstd([c1;c2]))];        % define yaxis limits

%% Plot configuration
%EjeX; % Time o EjeX

figure1 = grupo;
%figure%('visible',display_fig,'position', [0, 0, 1000, 500]); hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1 = shadedErrorBar(Xt,mean(data1,1),std(data1,[],1)/sqrt(size(data1,1)),'lineprops', color);hold on;
p2 = shadedErrorBar(Xt,mean(data2,1),std(data2,[],1)/sqrt(size(data2,1)),'lineprops', '-k');hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel([xMeassure],'FontSize',12,'FontWeight','bold')
ylabel([yMeassure],'FontSize',12,'FontWeight','bold')
%%%%%%%%%%%%%%%%%%%%%%%%%% AXIS SCALE %%%%
ylim(yaxis)
xlim ([Xt(1) Xt(end)])

title(figure1,'Fontsize',12,'FontWeight','bold','interpreter', 'none');
h = title(figure1,'Fontsize',12,'FontWeight','bold','interpreter', 'none');
P = get(h,'Position');
%set(h,'Position',[P(1) P(2)+0.1 P(3)])

ax = gca; % current axes
ax.XTickMode = 'manual';
ax.TickDirMode = 'manual';
ax.TickDir = 'in';
ax.XColor = 'black';
ax.YColor = 'black';
set(gca,'LineWidth',1,'Fontsize',10,'clipping', 'on')

%% Sombreado
% pmask = pvals;
winlim = start_point:1:end_sample;

for pv=1:length(pmask)-2   
    if pmask(pv) == 1
        p = patch('Faces', [1 2 3 4], 'Vertices', [Xt(winlim(pv)) yaxis(1); Xt(winlim(pv)) yaxis(2);...
            Xt(winlim(pv+1)) yaxis(2); Xt(winlim(pv+1)) yaxis(1)]);
        set(p,'facealpha',.2,'linestyle','none');
    end  
end

%p3 = plot([x1(1) x1(1)],[yaxis(1)+0.05 yaxis(2)-0.05],'color','k','LineWidth',2,'LineStyle','-');
hold on;
box on;

legend({condition1,condition2});
allChildren = get(gca, 'Children');                % list of all objects on axes
displayNames = get(allChildren, 'DisplayName');    % list of all legend display names
% Remove object associated with "data1" in legend
delete(allChildren(strcmp(displayNames, 'data1')))
%saveas(fig,char(strcat(ROIp,'.png')));
end