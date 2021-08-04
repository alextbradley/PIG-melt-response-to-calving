%Make figure 11 of the PIG calving manuscript: realistic domain melt rate (2009 geometry) and non-cumulative melt rate anomalies. NB: make_figure10 must be ran first to get the definition of the inner cavity.

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
figure(2); clf;
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
lambdashift = -101;
phishift = -0.07;
XX = XX + xshift;
YY = YY + yshift;
[phi, lambda] = polarstereo_inv(XX,YY);
lambda = lambda + lambdashift;
phi = phi + phishift;
lambda = double(lambda);
phi = double(phi);


%parameters
secs_per_year = 365.25*24*60*60;
density_ice = 918.0;

%time details
ntout1 = 11;
ntout2 = 12; %define time period to average over


%
% Generate data loop
%
if gendata
run_nos = ["078", "082", "083", "084", "085", "086"];
%run_nos = ["141", "142", "143", "144", "145", "146"];
sz = length(run_nos);

%setup storage
bsf_scenarios = cell(1,sz);
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

%velocities
UVEL_fname = strcat(rootdir, run_nos(i), '/run/stateUvel.nc');
UVEL = ncread(UVEL_fname, 'UVEL', [1,1,1,ntout1], [Inf, Inf, Inf, 1+ntout2 - ntout1]);
VVEL_fname = strcat(rootdir, run_nos(i), '/run/stateVvel.nc');
VVEL = ncread(VVEL_fname, 'VVEL', [1,1,1,ntout1], [Inf, Inf, Inf, 1 + ntout2 - ntout1]);
VVEL = squeeze(mean(VVEL,4));
UVEL = squeeze(mean(UVEL,4));
nt = length(ncread(UVEL_fname, 'TIME'))


%compute bsf
vvel = squeeze(sum(VVEL, 3)) * dz; %units m^2 /s
stream=zeros(size(vvel));
stream(nx,:)=vvel(nx,:)*dx;
for p=nx-1:-1:1
 stream(p,:)=stream(p+1,:) + vvel(p,:)*dx;
end
stream = stream/1e6; %convert to sv
bsf_scenarios{i} = stream;

end
end

%
% Make the plot
%
width = 0.30;
height = 0.44;
gapx = 0.02;
gapy = 0.02;
ncols = 3;
nrows = 2;
startx = (1 -width*ncols - (ncols-1)*gapx)/2;
starty = 1-(1 -height*nrows - (nrows -1)*gapy)/2;
starty = 0.97;
positions = zeros(4, nrows* ncols);
for i = 1:nrows*ncols
q = 1 + mod(i-1,ncols); %index in x direction
p = ceil(i/ncols); %index in y directio
positions(:,i) = [startx + (q-1)*gapx + (q-1)*width, starty - p*height - (p-1)*gapy, width, height];
end
cmap = lighter_blue_parula(100,0.1);


for q = 1:sz
icetopo = cell2mat(topo_scenarios(q));
bsf = cell2mat(bsf_scenarios(q));
topo = cell2mat(topo_scenarios(q));
bsf(bathy == 0) = nan;

bsf_sat = saturate(bsf, 0.45, -0.4);
%make plot
axs(q) = subplot('Position',positions(:,q));
contourf(lambda,phi,bsf_sat', 50, 'linestyle', 'none')
hold on
contour(lambda,phi,bathy', [0,0], 'k', 'linewidth', 1.5)
contour(lambda,phi,icetopo', [0,0], 'k', 'linewidth', 1.5)
colormap(axs(q),redblue);
if q == 1
c = colorbar;
c.Location = 'west';
c.Position(1) = positions(1,q) + 0.005; %relatvie horiztonal pos
c.Position(2) = 0.78;
c.Position(end) = 0.145;
c.Color = 0*[1,1,1]; %set color to white/blck
c.Label.String = 'BSF (Sv)';
c.Label.FontSize = 10;
end
xlim([-102.6, -99])
ylim([-75.45,-74.75])
set(axs(q),'Color',background_color)
        
%add the open ocean
bathynoice = bathy;
bathynoice(icetopo < 0) = 0; %remove cavity
bathynoice(bathy == 0) = 0; %remove grounded ice
bathynoice(bathynoice ~= 0) = .01; %make constant
%cf = contourf(lambda, phi, bathynoice',[.01,.01]);
%idx = find((cf(1,:) == 1));
%c1 = cf(1,2:idx(2)-1);
%c2 = cf(2,2:idx(2)-1);
%fill(c1,c2,ocean_color, 'linewidth', 1.5, 'edgecolor', 'k');
        
%add 750m bathymetric contour
%contour(lambda, phi, -1e-2*bathy', -1e-2*[-750, -750], 'color',bathycontourcolor, 'linewidth' ,1, 'linestyle', '--'); %he weird "* -1e-2" is to get bring -750 into range of colourar, avoid skewing it
grid on
axs(q).YTick = -75.4:0.1:-74.8;
if mod(q,3) == 1
axs(q).YTickLabel = {"24'", "18'",  "12'", "6'",  "75S", "54'", "48'"};
else
axs(q).YTickLabel = cell(length(axs(q).YTick),1);
end

axs(q).XTick = [-102.5:0.5:-99];
if floor((q-1)/3) == 1
axs(q).XTickLabel = {"30'", "102W",  "30'", "101W", "30'", "100W", "30'", "99W"};
else
axs(q).XTickLabel = cell(length(axs(q).XTick),1);
end

%add the definitions of the two cavity regions
realistic_inner_cavity_definition; %bring inner cavity definition into scope (a1,b1,a2,b2)
[XXns,YYns] = meshgrid(X,Y);

in1 = inpolygon(XXns',YYns', a1,b1);
in2 = inpolygon(XXns',YYns', a2,b2);
idx1 = (icetopo < 0) & in1;
idx2 = (icetopo < 0) & in2;
if ((q == 1) || (q == 6))
bsf_copy = bsf;
bsf_copy(idx1) = 0.01;
bsf_copy(~idx1) =0;
contour(lambda,phi,bsf_copy',[.01,.01], '--', 'linewidth', 1.75, 'linecolor', 'k')
bsf_copy = bsf_sat;
bsf_copy(idx2) = 0.01;
bsf_copy(~idx2) =0;
contour(lambda,phi,bsf_copy',[.01,.01], '--', 'linewidth', 1.75, 'linecolor', 'k')

end


end %end loop over runs
fig = gcf; fig.Position(3:4) = [1280, 687.333];

%
% Save
%
if save_flag
saveas(gcf, "plots/figure13.eps", "epsc")
end
