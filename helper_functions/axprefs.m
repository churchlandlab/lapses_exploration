function axprefs(ax,fontSize,type)
if nargin == 1
    fontSize = 18;
    type = 'none';
elseif nargin ==2
    type = 'none';
end

% fix up axes the way we like them
set(ax,'TickDir','out','FontSize',fontSize,'Box','off','fontname','Helvetica');
set(ax,'ticklength',[0.02 0.02]);
%     offsetAxes
hglobal = gcf;
hglobal.PaperPositionMode = 'auto';
fig_pos = hglobal.PaperPosition;
hglobal.PaperSize = [1.1*fig_pos(3) 1.1*fig_pos(4)];

switch type
    case 'pmf'
        axis square
        
xlabel('Stimulus rate (Hz)');
ylabel('p(Chose high)');
set(gca,'XTick',[9,16]);
xlim([7,18])
ylim([0,1])
set(gca,'YTick',[0,1]);
% set(gcf,'Position',[1294 549 387 406]);
end
end



