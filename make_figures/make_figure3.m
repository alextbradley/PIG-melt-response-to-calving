%Make figure 3 in the Idealized PIG calving manuscript: results from the default run.
%(a) Plot of melt rate with boundary current overlain
%(b) Plot of the water column thickness with barotropic stream function overlain
%(c) Plot of the bottom temperature with bottom current overlain
%
% NB: Many of the data files referred to in this script are too large to be hosted online. These files are hosted internally as BAS.
% Please email Alex Bradley (aleey@bas.ac.uk) to obtain a copy.
%Alex Bradley (aleey@bas.ac.uk) 27/05/2021. MIT license.

%
% Flags
%
gendata = 1; %specify whether to pass through the generate data loop
save_flag = 1;

%
% Preliminaries
%
addpath("plot_tools");
plot_defaults
label_size = 11;
ax_fontsize = 10;
figure(1); clf;
fig = gcf; fig.Position(3:4) = [900, 500];

%
% Data locations
%
rootdir = '/data/oceans_output/shelf/aleey/mitgcm/APIGi_';
topodir = '/data/hpcdata/users/aleey/mitgcm/matlab/interp_APIGi/topo_files';
bathy_path = '/data/hpcdata/users/aleey/mitgcm/matlab/interp_APIGi/bathy_files/bathymetry.shice';

%grid details
nx=120; % number of grid cells along longitudinal direction
ny=320; % number of grid cells along latitudinal direction
nz=110; % number of vertical grid cells
dx=400;
dy=400;
dz=10;
X = 0:dx:(nx-1)*dx;
Y = 0:dx:(ny-1)*dy;
Z = 0:dz:(nz-1)*dz;
[XX,YY] = meshgrid(X,Y);
YYt = YY';
idx = (YYt < 30e3); %inner cavity definition

%parameters
secs_per_year = 365.25*24*60*60;
density_ice = 918.0;
lambda = 7.61*1e-4;%constants in liquidus
gamma = 5.73*1e-2;
T0 = 8.32*1e-2;

%time details
ntout1 = 6;
ntout2 = 12; %define time period to average over


%
% Generate data loop
%
run_nos = "077"; 
sz = length(run_nos);
extent = 84;
H = 400; %ridge height (always 400);
W = 100; %ridge gap

%generate data loop
if gendata

%load bathy
fid = fopen(bathy_path);
bathy = fread(fid, 'real*8', 'b');
bathy = reshape(bathy, [nx, ny]);
%bathy(bathy == 0) = nan;

%draft
topo_fname=  ['shelfice_topo_H' num2str(H) '_W' num2str(W) '_extent' num2str(extent) 'km.bin'];
topo_fid = fopen(strcat(topodir, '/',topo_fname));
topo = fread(topo_fid, 'real*8', 'b');
topo = reshape(topo, [nx, ny]);

%melt rates
state2D_fname = strcat(rootdir, run_nos, '/run/state2D.nc');
melt = ncread(state2D_fname, 'SHIfwFlx', [1, 1, ntout1], [Inf, Inf, 1+ntout2- ntout1]);
melt = mean(melt, 3); %average over months ntout1 to ntout2
melt = -melt * secs_per_year / density_ice;
melt(topo == 0) = nan;

%Theta
Theta_fname = strcat(rootdir, run_nos, '/run/stateTheta.nc');
Theta = ncread(Theta_fname, 'THETA', [1,1,1,ntout1], [Inf, Inf, Inf, 1+ntout2 - ntout1]);
Theta = mean(Theta, 4);

%Salinity
Salt_fname = strcat(rootdir, run_nos, '/run/stateSalt.nc');
Salt = ncread(Salt_fname, 'SALT', [1,1,1,ntout1], [Inf, Inf, Inf, 1+ntout2 - ntout1]);
Salt = mean(Salt, 4);

%Velocities
UVEL_fname = strcat(rootdir, run_nos, '/run/stateUvel.nc');
UVEL = ncread(UVEL_fname, 'UVEL', [1,1,1,ntout1], [Inf, Inf, Inf,  1+ntout2 - ntout1]);
UVEL = mean(UVEL, 4);
VVEL_fname = strcat(rootdir, run_nos, '/run/stateVvel.nc');
VVEL = ncread(VVEL_fname, 'VVEL', [1,1,1,ntout1], [Inf, Inf, Inf,  1+ntout2 - ntout1]);
VVEL = mean(VVEL, 4);


%boundary layer quantities
Nb = 3; %number of grid pts to take mean over
Sbl = nan(nx,ny); Tbl = nan(nx,ny); Ubl = nan(nx, ny); Vbl = nan(nx,ny);
Sbot = nan(nx,ny); Tbot = nan(nx,ny); Ubot = nan(nx,ny); Vbot = nan(nx,ny);
for p = 1:nx
for q = 1:ny
	%work out the bottom and velocity
	idx = find((bathy(p,q) - (-Z) < 0), 1, 'last'); %gives you the index of first Z grid point above the bathymetry
         Tbot(p,q) = double(mean(Theta(p,q,idx-Nb+1:idx)));
         Ubot(p,q) = double(mean(UVEL(p,q,idx-Nb+1:idx)));
         Vbot(p,q) = double(mean(VVEL(p,q,idx-Nb+1:idx)));


        if topo(p, q) < 0 %if we're in the cavity
                idxtop = find((topo(p,q) - (-Z)) > 0, 1, 'first'); %gives you the index of first Z grid point above the bathymetry
                idxtop = find(Theta(p,q,:) ~= 0);
                idxtop = idxtop(1);
                Sbl(p,q) = double(mean(Salt(p,q,idxtop:idxtop+Nb-1)));
                Tbl(p,q) = double(mean(Theta(p,q,idxtop:idxtop+Nb-1)));
                Ubl(p,q) = double(mean(UVEL(p,q,idxtop:idxtop+Nb-1)));
                Vbl(p,q) = double(mean(VVEL(p,q,idxtop:idxtop+Nb-1)));

                if 1 %account for partial cell in the mean calculation
                draft = topo(p,q);
                partial_cell_frac = abs(rem(draft, dz)) / dz;
                draft_rounded = draft + abs(rem(draft, dz));
                [~,idx_top] = min(abs(-Z - draft_rounded));
                vec = [partial_cell_frac,1,1]';
                Sbl(p,q) = sum(vec.*squeeze(Salt(p,q,idxtop:idxtop+Nb-1)))/sum(vec);
                Tbl(p,q) = sum(vec.*squeeze(Theta(p,q,idxtop:idxtop+Nb-1)))/sum(vec);

                Ubl(p,q) = sum(vec.*squeeze(UVEL(p,q,idxtop:idxtop+Nb-1)))/sum(vec);
                Vbl(p,q) = sum(vec.*squeeze(VVEL(p,q,idxtop:idxtop+Nb-1)))/sum(vec);
                end

        end
end %end loop over y grid
end %end loop over x grid

%compute BSF
vvel = squeeze(sum(VVEL, 3)) * dz; %units m^2 /s
stream=zeros(size(vvel));
stream(nx,:)=vvel(nx,:)*dx;
for p=nx-1:-1:1
 stream(p,:)=stream(p+1,:) + vvel(p,:)*dx;
end
stream = stream/1e6; %convert to sv

end %end gendata loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%
width = 0.22;
gap = 0.02;
startx = (1 - 3*width + 2*gap)/2 - 0.05;
starty = 0.1;
height = 0.8;
positions = [startx, starty, width, height;
	     startx + gap + width, starty, width, height;
	     startx + 2*gap + 2*width, starty, width, height];
	    
%
% Plot 1: Melt rate and BL velocities
%
ax(1) = subplot('Position', positions(1,:)); hold on; box on
contourf(X/1e3,Y/1e3,melt', 50, 'linestyle', 'none');
cmap = lighter_blue_parula(100,0.2); 
colormap(ax(1), cmap);

%add the velocity arrows
idxX = 1:8:120;
idxY = 1:8:320;
[XX,YY] = meshgrid(X,Y);
XX = XX/1e3; YY = YY/1e3;
velscale =15;
plot(X/1e3, 50*ones(length(X), 1), 'w--', 'linewidth', 1.5)
quiver(XX(idxY, idxX),YY(idxY, idxX),velscale *Ubl(idxX, idxY)', velscale*Vbl(idxX, idxY)', 'autoscale', 'off', 'color', 'k')

%colorbar and arrow
c = colorbar;
c.Location = 'north';
c.Label.String = 'melt rate (m/yr)';
c.Position(end) = c.Position(end) - 0.02;
c.Position(2) = c.Position(2) + 0.03;
c.FontSize = 10;
plot([20, 30], [100,100], 'k', 'linewidth', 1);
text(31, 100, '0.6 m/s')
xlabel('x (km)');
ylabel('y (km)')


%
% Plot 2: 1/h and BSF contours
%
ax(2) = subplot('Position', positions(2,:)); hold on; box on;
column_thickness = topo - bathy;
contourf(X/1e3,Y/1e3,1e3* (1./column_thickness)', 20, 'linestyle', 'none');
plot(X/1e3, 50*ones(length(X), 1), 'w--', 'linewidth', 1.5)
colormap(ax(2), cmap);
c = colorbar;
c.Location = 'north';
c.Position(end) = c.Position(end) - 0.02;
c.Position(2) = c.Position(2) + 0.03;
c.Label.String = '1/h (10^{-3} m^{-1})';
c.FontSize = 10;
yticks([])
xlabel('x (km)');

%add bsf
streamsm = smooth2a(stream, 2,2);
streamsm(:,end-32:end) = nan;
axnew = axes;
axnew.Position = ax(2).Position;
[C,h] =contour(X/1e3,Y/1e3, streamsm', [-0.7, -0.5, -0.3, -0.1], 'k');
clabel(C,h);
hold on
streamsm(1:4,:) = nan; streamsm(end-3:end,:) = nan; streamsm(:,1:20) = nan; streamsm(:,end-4:end) = nan; %remove borders and near Gl where stream is messy
[C,h] =contour(X/1e3,Y/1e3, streamsm', [0,0], 'k');
clabel(C,h);

xticks([]);
yticks([]);
set(axnew, 'color', 'none')

%
% Plot 3
%
ax(3) = subplot('Position', positions(3,:)); hold on; box on
Tbotsm = smooth2a(Tbot, 2,2);
contourf(X/1e3,Y/1e3,Tbotsm', 50, 'linestyle', 'none');
cmap = lighter_blue_parula(100,0.2); 
colormap(ax(3), cmap);

plot(X/1e3, 50*ones(length(X), 1), 'w--', 'linewidth', 1.5)
c = colorbar;
c.Location = 'north';
c.Position(end) = c.Position(end) - 0.02;
c.Position(2) = c.Position(2) + 0.03;
c.Label.String = 'Bottom temp(\circC)';
c.FontSize = 10;
yticks([])
Ubot(:,end-48:end)= nan;
quiver(XX(idxY, idxX),YY(idxY, idxX),velscale *Ubot(idxX, idxY)', velscale*Vbot(idxX, idxY)', 'autoscale', 'off', 'color', 'k')
xlabel('x (km)')

%
% Save
%
if save_flag
saveas(gcf, "plots/figure3", 'epsc')
end

