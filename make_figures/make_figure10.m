%Make figure 10 of the IG calving manuscript: realistic domain melt rate (2009 geometry) and non-cumulative melt rate anomalies

% NB: Many of the data files referred to in this script are too large to be hosted online. These files are hosted internally as BAS.
% Please email Alex Bradley (aleey@bas.ac.uk) to obtain a copy.
%Alex Bradley (aleey@bas.ac.uk) 27/05/2021. MIT license.

%
% Flags
%
gendata = 1; %specify whether to pass through the generate data loop
save_flag = 0;

%
% Preliminaries
%
addpath("plot_tools");
plot_defaults
label_size = 11;
ax_fontsize = 11;
figure(1); clf;
fig = gcf; fig.Position(3:4) = [1085, 540];
%plot colours
background_color = 0.8*[1,1,1]; %ice color
background_color = [1,1,1];
ocean_color = [226,249, 255]/255;
cols = [7,43, 184;
        184, 7, 43]/255; %colours for each different bc
bathycontourcolor = [255,0,255]/255;

%
% Data locations
%
rootdir = '/data/oceans_output/shelf/aleey/mitgcm/rPIG_'; %output data NOT in github repo (contact for copy)
topodir = '../gendata_realistic/topo_files/';
bathypath = '../gendata_realistic/bathy_files/bathymetry.shice';

%grid details
nx=360; % number of grid cells along longitudinal direction
ny=320; % number of grid cells along latitudinal direction
nz=120; % number of vertical grid cells
dx=400;
dy=400;
dz=10;
X = ncread(strcat(rootdir, "078", '/run/state2D.nc'), 'LONGITUDE');
Y = ncread(strcat(rootdir, "078", '/run/state2D.nc'), 'LATITUDE');
[XX,YY] = meshgrid(X,Y);
YYt = YY';

%put onto the latlon grid
xshift = 870;
yshift = -15000;
XX = XX + xshift;
YY = YY + yshift;
[phi, lambda] = polarstereo_inv(XX,YY);
lambda = lambda - 101;
phi = phi - 0.07;
lambda = double(lambda);
phi = double(phi);


%parameters
secs_per_year = 365.25*24*60*60;
density_ice = 918.0;

%time details
ntout1 = 6;
ntout2 = 12; %define time period to average over


%
% Generate data loop
%
run_nos = ["078", "082", "083", "084", "085", "086"];
sz = length(run_nos);

%setup storage
melt_scenarios = cell(1,sz);
topo_scenarios = cell(1,sz);

%load bathy
bathyfid = fopen(bathypath);
bathy = fread(bathyfid, 'real*8', 'b');
bathy = reshape(bathy, [nx,ny]);
bathy = double(bathy);

%loop over runs
for i = 1:sz
%draft
topo_fname=  ['shelfice_topo_scn', num2str(i), '.shice'];
topo_fid = fopen(strcat(topodir, '/',topo_fname));
topo = fread(topo_fid, 'real*8', 'b');
topo = reshape(topo, [nx,ny]);
topo_scenarios{i} = topo;


%melt rates
state2D_fname = strcat(rootdir, run_nos(i), '/run/state2D.nc');
melt = ncread(state2D_fname, 'SHIfwFlx', [1, 1, ntout1], [Inf, Inf, 1+ntout2- ntout1]);
melt = mean(melt, 3); %average over months ntout1 to ntout2
melt = -melt * secs_per_year / density_ice;
melt_scenarios{i} = melt;

end


%
% Make the plot
%
for q = 1:sz
icetopo = cell2mat(topo_scenarios(q));
melt = cell2mat(melt_scenarios(q));

%if baseline, store as such. Otherwise, get the anomaly
if q == 1
melt_diff = melt; 
else
melt_diff = melt - cell2mat(melt_scenarios(q-1)); %
end
melt_diff(icetopo == 0) = nan;
if q == 1
melt_diff_sat = saturate(melt_diff, 126,0);
else
melt_diff_sat = saturate(melt_diff, 22,-20);
end


%make plot
axs(q) = subplot(2,3,q);
contourf(lambda,phi,melt_diff_sat', 20, 'linestyle', 'none')
hold on
contour(lambda,phi,icetopo', [0,0], 'k', 'linewidth', 1.5)
if q == 1
        colormap(axs(q), parula);
else
        colormap(axs(q),redblue);
end
if q == 1 || q==2
c = colorbar;
c.Location = 'west';
c.Position(1) = c.Position(1)  - 0.005; %relatvie horiztonal pos
c.Position(2) = 0.77+0.004;
c.Position(end) = 0.145;
c.Color = 0*[1,1,1]; %set color to white/blck
if q == 1
    c.Label.String = 'melt rate (m/yr)';
else
    c.Label.String = 'anomaly (m/yr)';
end
end
xlim([-102.6, -99])
ylim([-75.45,-74.75])
set(subplot(2,3,q),'Color',background_color)
        
%add the open ocean
bathynoice = bathy;
bathynoice(icetopo < 0) = 0; %remove cavity
bathynoice(bathy == 0) = 0; %remove grounded ice
bathynoice(bathynoice ~= 0) = 1; %make constant
cf = contourf(lambda, phi, bathynoice',[1,1]);
idx = find((cf(1,:) == 1));
c1 = cf(1,2:idx(2)-1);
c2 = cf(2,2:idx(2)-1);
fill(c1,c2,ocean_color, 'linewidth', 1.5, 'edgecolor', 'k');
        
%add 750m bathymetric contour
%contour(lambda, phi, -1e-2*bathy', -1e-2*[-750, -750], 'color',bathycontourcolor, 'linewidth' ,1, 'linestyle', '--'); %he weird "* -1e-2" is to get bring -750 into range of colourar, avoid skewing it
grid on
%yticks([-75.25, -75, -74.75])
end %end loop over runs
fig = gcf; fig.Position(3:4) = [1280, 687.333];
