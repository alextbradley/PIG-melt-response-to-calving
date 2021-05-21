%Make figure 2 in the manuscript: schematic of the shelf geometry and hydrographic conditions
%Alex Bradley (aleey@bas.ac.uk), 21/05/2021

%
% Preliminaries
%
plot_defaults
fig = figure(1); clf; 
%fig.Position(3:4) = 

%
% schematic of ice shelf geometry
%
yy = 0:400:(319*400);
yyf = max(yy) - yy;
[~, idx] = min(abs(yy - 84*1e3));
H = [100, 150, 200];
linestyles = ["-","--", "--"];
h_profiles = zeros(3, length(yy));
for i = 1:3
h_profiles(i,:)=(310 + H(i))/2.64*atan(0.17*yy/1000 - 3) + 0.47*(H(i)+400) - 1051.3;
h_profiles(i, idx+1:end) = nan;
end

pos = [0.1, 0.12, 0.5, 0.82];
subplot('Position', pos); box on; hold on
for i = 1:3
if i == 1
fillX = [yyf(1:idx), flip(yyf(1:idx))]/1e3;
fillY = [h_profiles(i,1:idx), -300*ones(1,idx)];
fill(fillX, fillY, [173, 216, 230]/255, 'linewidth', 1.5)
else
plot((max(yy) - yy)/1e3, h_profiles(i,:), 'k', 'linestyle', linestyles(i))
end
end

%add the ridge
fillX = [yy, flip(yy)]/1e3;
latg = [1.62e6:400:1.748e6-400];
bump = 400*exp(-(latg-1.67e6).^2/(2*12000^2)) - 1095;
fillY = [-1100*ones(1,length(yy)), bump];
fill(fillX, fillY, [181, 101, 29]/255, 'Linewidth', 1.5)
%plot(yyf/1e3, bump, 'k', 'linewidth', 1.5)	


xlabel('Y (km)')
ylabel('depth(m)')
xlim([min(yy), max(yy)]/1e3)
ylim([-1100, -400])
north = text(5,-550, "North", 'FontSize', 16);
set(north, 'Rotation', 90)

south = text(123,-550, "South", 'FontSize', 16);
set(south, 'Rotation', 90)

