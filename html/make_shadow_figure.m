function fig = make_shadow_figure(dat, dn, bw)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

dat = dat';  % transpose for a more natural image, with time of day along x axis

M = size(dat,2);  % number of times, number of days in the data
day_idx = 1:M;

fig = figure;
%     fig.PaperSize = [3.325 2.325];
%     fig.PaperPosition = [.25 .125 2.875 2.125];
%     fig.Position = [200 200 600 600];

% handle NaNs
tmp=dat;
tmp(isnan(dat))= -1;

caxis([min(tmp(:)) max(tmp(:))]);

if bw==1  
    colormap(flip(gray))
elseif bw==0
    tmpcolormap = [0.75 0.75 0.75; 0 0 0; parula(max(dat(:)))];
    colormap(tmpcolormap); 
else
    tmpcolormap = [0 0 0; parula(max(dat(:)))];
    colormap(tmpcolormap); 
end

image(tmp);

% set up labels for time and date axes
% date axis
dates = unique(floor(dn));  % datenum for each day in time list
dv = datevec(dates);
u = dv(:,3)==1;  % tick at the first day of each month
yticksat = day_idx(u);

tmp = dn(floor(dn)==min(floor(dn))) - min(floor(dn));  % get fractional day on first day
tmp = 24*tmp;
xtickhours = [3 6 9 12 15 18 21]';
xticksat = NaN(size(xtickhours));
for i=1:length(xtickhours)
    xticksat(i) = find(tmp>=xtickhours(i), 1, 'first');
end
xticklbls = dn(xticksat);

%     image(dat);

%     if any(any(isnan(dat)))
%     end

xlim([330 1130])
%ylim([0.5 183.5])

axis square
ax = gca;
ax.FontSize = 10;
ax.YDir = 'normal';
ax.YTick = yticksat;
ax.YTickLabel = datestr(dates(u),'mm/dd');
ax.XTick = xticksat;
ax.XTickLabel = datestr(xticklbls,'HH:MM');
ylabel('Date')
xlabel('Local Time')

end

