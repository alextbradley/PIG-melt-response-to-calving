%Make figure 1 in the PIG calving manuscript: plots of (a) PIG bathymetry with 2009 and 2020 ice front on (b) ice topo and bathymetry along the contour, (c) gap along the contour

%Alex Bradley (aleey@bas.ac.uk) 27/05/2021. MIT license.

%
% Plot info
%
addpath('plot_tools');
positions = [0.1,0.1,0.4,0.85;
             0.6, 0.27, 0.38, 0.24;
             0.6, 0.54, 0.38, 0.24];

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
xidx = 245*ones(1,length(yidx));
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
subplot('Position',positions(1,:)); imshow(imageData);

%
% add bathymetry and fronta
%
hold on
contourf(xi, yi, bathyng', 30, 'linestyle', 'none');
[cgl,~] = contour(xi,yi,bathy', [0,0], 'k', 'linewidth',2 );
[c2009,~] = contour(xi,yi,topo2009', [0,0], 'linecolor', plotcolor3, 'linewidth', 2);
[c2020,~] = contour(xi,yi,topo2020', [0,0], 'linecolor', plotcolor1, 'linewidth', 2);


%
% add cross section
%
cline = [xline; yline];
cline_img = model2image(cline);
plot(cline_img(1,:), cline_img(2,:), 'k--');

ylim([3000, 14847])
xlim([2000,12000])
camroll(-90);

subplot('Position', positions(2,:)); 
plot(sline/1e3, bathyline, 'color', plotcolor1, 'linewidth', 2);
hold on
plot(sline/1e3, topoline, 'color', plotcolor3, 'linewidth', 2);
xlim([0, 46])

subplot('Position', positions(3,:));
plot(sline/1e3,topoline - bathyline, 'color', plotcolor1, 'linewidth', 2);
ax = gca; ax.XTickLabels = cell(length(ax.XTick), 1);
xlim([0, 46])

fig = gcf; fig.Position(3:4) = [900, 600];
