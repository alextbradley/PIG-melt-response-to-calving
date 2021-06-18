%Make figure 1 in the PIG calving manuscript: plots of (a) PIG bathymetry with 2009 and 2020 ice front on (b) ice topo and bathymetry along the contour, (c) gap along the contour

%Alex Bradley (aleey@bas.ac.uk) 27/05/2021. MIT license.

%
% Plot info
%
saveflag = 1;
addpath('plot_tools');
plot_defaults
positions = [0.1,0.1,0.4,0.85;
             0.58, 0.31, 0.38, 0.25;
             0.58, 0.59, 0.38, 0.25];

close all

%
% Grid
%
nx = 360;
ny = 320;
dx = 400;
dy = 400;
x = dx:dx:nx*dx;
y = dy:dy:ny*dy;
[xx,yy] = meshgrid(x,y);


%
% get grid in image space
%
load('image_to_model_points.mat', 'ximage', 'yimage', 'xmod', 'ymod');
[model2image, image2model] = genmaps_image2model(ximage, yimage, xmod, ymod); %maps from model to image and vice 
yi = zeros(size(yy));
xi = zeros(size(xx));
for i = 1:320
cc = [xx(i,:); yy(i,:)];
ci = model2image(cc);
yi(i,:) = ci(2,:);
xi(i,:) = ci(1,:);
end


% 
% Bathy and topos
%
fid = fopen("../gendata_realistic/bathy_files/bathymetry.shice");
bathy = fread(fid, 'real*8', 'b');
bathy = reshape(bathy, [nx,ny]);
bathyng = bathy; bathyng(bathy ==0) = nan;

fid = fopen("../gendata_realistic/topo_files/shelfice_topo_scn1.shice");
topo2009 = fread(fid, 'real*8', 'b');
topo2009 = reshape(topo2009, [nx,ny]);

fid = fopen("../gendata_realistic/topo_files/shelfice_topo_scn2.shice");
topo2020 = fread(fid, 'real*8', 'b');
topo2020 = reshape(topo2020, [nx,ny]);

% 
% engineer topos so that gl doesn't show up at zero contour
%
for p = 2:nx-1
    for q = 2:ny-1
        if any( bathy(p,q+1) == 0 || bathy(p, q-1) == 0 || bathy(p+1,q) == 0 || bathy(p-1,q) == 0)
            topo2009(p,q) = nan;
            topo2020(p,q) = nan;
        end
    end
end

%
% plot bathy and topos and get data
%
figure(1); clf; hold on
contourf(x,y,bathyng', 20, 'linestyle', 'none')
[cgl,~] = contour(x,y,bathy', [0,0], 'k');
[c2009,~] = contour(x,y,topo2009', [0,0], 'r');
[c2020,~] = contour(x,y,topo2020', [0,0], 'b');


%
%transect
%
yidx = 5:140;
xidx = 256*ones(1,length(yidx));
sline =  sqrt((x(xidx) - x(xidx(1))).^2 + (y(yidx) - y(yidx(1))).^2); %arclength along line
xline = x(xidx);
yline = y(yidx);
plot(xline, yline, 'k--', 'linewidth', 2);

%
%get bathy and topo along line
%
bathyline = nan(1,length(yidx));
topoline  = nan(1,length(yidx));
for i = 1:length(bathyline)
if bathy(xidx(i), yidx(i)) ~=0
bathyline(i) = bathy(xidx(i), yidx(i));
end
if topo2009(xidx(i), yidx(i)) ~=0
topoline(i)  = topo2009(xidx(i), yidx(i));
end
end
idx = ~isnan(bathyline);
sline = sline(idx);
sline = sline - sline(1);
bathyline = bathyline(idx);
topoline = topoline(idx);
%
% Open image
%
 
t = Tiff('PIG-S2-NovDec2020.tif', 'r'); %!! not in git repo!!
imageData = read(t);
figure(2); 
ax(1) = subplot('Position',positions(1,:)); imshow(imageData);

%
% add bathymetry and fronta
%
hold on
contourf(xi, yi, bathyng', 30, 'linestyle', 'none'); 
c = colorbar('Location', 'northoutside');
c.Position(2) = 0.79;
c.Position(4) = 0.02;

c.Label.String = 'sea bed depth (m)';
c.Label.VerticalAlignment = 'bottom';
c.Label.FontSize = 12;
c.Label.Interpreter = 'latex';

[cgl,~] = contour(xi,yi,bathy', [0,0], 'k', 'linewidth',2 );
[c2009,h2009] = contour(xi,yi,topo2009', [0,0], 'linecolor', plotcolor3, 'linewidth', 2);
[c2020,h2020] = contour(xi,yi,topo2020', [0,0], 'linecolor', plotcolor1, 'linewidth', 2);
95

%
% add cross section
%
cline = [xline; yline];
cline_img = model2image(cline);
plot(cline_img(1,:), cline_img(2,:), 'k--');
ptA = text(ax(1), 9700,8200, 'A', 'FontSize',12, 'FontWeight', 'bold');
ptB = text(ax(1), 5300,8900, 'B', 'FontSize',12, 'FontWeight', 'bold');

%
% Label fronts
%
f2009 = text(ax(1), 4500,13500, '2009', 'FontSize',12, 'color', plotcolor3);
f2020 = text(ax(1), 8700,10200, '2020', 'FontSize',12, 'color', plotcolor1);

%tidy
ylim([3000, 14847])
xlim([2000,12000])
plot(ax(1), [11000,11000], [5000,4000], 'k', 'linewidth', 3)
scalebar = text(ax(1),10950, 6600, '10 km', 'Interpreter', 'latex', 'FontSize', 12);


camroll(-90);

%
% Plot along the cross section
%
ax(2) = subplot('Position', positions(3,:)); grid on 
plot(sline/1e3, bathyline, 'color', plotcolor1, 'linewidth', 2);
hold on
plot(sline/1e3, topoline, 'color', plotcolor3, 'linewidth', 2);
xlim([0, 45])
ax(2).XTickLabels = cell(length(ax(2).XTick), 1);

ax(3) = subplot('Position', positions(2,:)); grid on
plot(sline/1e3,topoline - bathyline, 'color', plotcolor1, 'linewidth', 2);
xlim([0, 45])
ptAA = text(ax(2), 1,-700, 'A', 'FontSize', 12, 'FontWeight', 'bold');
ptBB = text(ax(2), max(ax(3).XLim)-3 ,-700, 'B', 'FontSize', 12, 'FontWeight', 'bold');



%
% tidy everything
%
ax(3).XLabel.String = 'Distance along transect (km)';
ax(3).XLabel.Interpreter = 'latex';
ax(3).XLabel.FontSize = 12;
ax(2).YLabel.String = 'depth (m)';
ax(2).YLabel.FontSize = 12;
ax(2).YLabel.Interpreter = 'latex';
ax(2).YLim = [-850, -200];
ax(2).YTick = [-800, -600, -400,-200];
ax(3).YTick = [0,100,200,300, 400];
ax(3).YLabel.String = 'gap (m)';
ax(3).YLabel.FontSize = 12;
ax(3).YLabel.Interpreter = 'latex';

fig = gcf; fig.Position(3:4) = [900, 600];
grid(ax(3), 'on')
grid(ax(2), 'on')

%
% Save?
%
if saveflag
saveas(gcf, 'plots/figure1.eps', 'epsc');
end
