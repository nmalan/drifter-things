%Plots and calculates shear information from ACT cruise shipboard ADCP for use in stability analysis for drifter paper
% N Malan, UCT, May 2016

close all;clear all;clc
%set drifter release point
drift_lat=-33.598;
drift_lon=27.607;

%load SADCP
load act0213_wh300.mat

%plot surface velocity section
	f=300
  	l=500
	surf_U=squeeze(adp.u(1,:));
	surf_V=squeeze(adp.v(1,:));

m_proj('mercator','longitude',[25.5 29],'latitude', [-35 -33])
m_vec(20,adp.lon(f:l),adp.lat(f:l),surf_U(f:l),surf_V(f:l),[0 1 127/255],'shaftwidth',0.01,'headlength',0.1,'headwidth',0.1 )
hold on
%plot drifter release location
m_plot(drift_lon,drift_lat,'rx','markersize',10,'linewidth',3)
m_grid('box','fancy')
m_gshhs_i('patch',[244/255 164/255 96/255])

%put in scale arrow
m_vec(1,25.7,-34.9,1,0,[0 1 127/255],'shaftwidth',0.01,'headlength',0.1,'headwidth',0.1 )
m_text(25.75,-34.85,'0.5 m/s')
m_text(25.6,-33.1,'a)','FontSize',18)
%add bathymetry
clevs=[-200];
[fat guy]=m_tbase('contour',-200:-200:-600,'linewidth',1,'linecolor',[.5 .5 .5]);
clabel(fat,guy,clevs)

%print -depsc S_ADCP_track.eps
%print -r300 -dpng S_ADCP_track.png
%
% %plot map of section location
% figure(2)
%  m_proj('mercator','longitude',[17 30],'latitude', [-46 -23])
%  m_etopo2
%   %m_scatter(adp.lon,adp.lat,adp.tday)
%  m_grid('box','fancy')
%  m_gshhs_i('patch',[.7 .7 .7])

%plot section

%get direction and distance of section
section_dir=azimuth(adp.lat(l),adp.lon(l),adp.lat(f),adp.lon(f))
section_dist=m_lldist(adp.lon(f:l),adp.lat(f:l));
section_dist=permute(section_dist,[2 1]);
section_d=section_dist
section_dist=cumsum([0 section_dist]);
%rotate velocities
[cross_u,cross_v]=uv_rotate(surf_U(f:l),surf_V(f:l),(180-section_dir));
cross_u=fixgaps(cross_u);

%plot crosstrack
figure(2)
% plot(-section_dist,-cross_u,'k','linewidth',2)
% xlim([-260 0])
% ylabel('crosstrack velocity [m/s]')
% xlabel('distance along transect [km]')

% %ax.xtick([-250 -200 -150 -100 -50 0])
% set(gca,'XTickLabel',{'0','50','100','150','200','250'})
% hold on
% %plot drifter release point
% plot(-section_dist(183),-cross_u(183),'rx','markersize',10,'linewidth',3)
%print -depsc SADCP_drifter.eps
%plot 2D section
U = adp.u(:,f:l);
V = adp.u(:,f:l);
[U_rot_2d,V_rot_2d]=uv_rotate(U,V,(180-section_dir));
cross_u=fixgaps(cross_u);

pcolor(-section_dist,adp.depth,-U_rot_2d)


%% Now calculate vorticity along the section

%create 1km resolution distance grid
grid_dist = [0:1:250];

%use geometry on angle of section - each 1km along section = dx899m and dy438m
dx = 899
dy = -438

%find nearest velocity values
ind = nearestpoint(grid_dist,section_dist);
u = surf_U(f:l);
v = surf_V(f:l);
grid_v = v(ind);
grid_u = u(ind);
grid_cross = cross_u(ind);
grid_dist=fliplr(grid_dist)

figure(3)
subplot(2,1,1)
plot(grid_dist,-(filtreNaN_lowpass(grid_cross,10)))
title('crosstrack velocity [m/s]')

%calculate vorticity dvdx-dudy
f=2*(7.2921*10^-4)*sind(-34)

du=diff(grid_u);
dv=diff(grid_v);

for i = 1:length(grid_dist)-1;
vtcty(i) = (dv(i)/dx-du(i)/dy)/f
end

subplot(2,1,2)
plot(grid_dist(1:end-1),filtreNaN_lowpass(vtcty,5))
title('crosstrack rel. vorticity/f')

%% Now, let's calculate Rayleigh's condition of instability (beta-Uyy) must change sign in the domain.

%First calculate beta:
beta0 = 2*(7.2921*10^-5/6371000)*cos(34)
%Then take the second derivative of the velocity profile. (i.e. the gradient of vorticity)
smth_vtcty = fixgaps(filtreNaN_lowpass(vtcty,5));
Uyy = gradient(smth_vtcty);

Ray = beta0-Uyy;

%get value for U_s (where Ray vanishes to 0)
Us = crossing(Ray);

%use biggest zero crossing - index 22
Us = grid_u(22);

% now calculate Fjortift's criterion:
for i = 1:length(Uyy)
	Fjor(i) = (Ray(i))*(grid_u(i)-Us)
end

%print -dpng -r300 adcp_relvorticity_raw.png

figure(4)
axes('position',[0.1 0.35  0.6 0.55])
plot(-section_dist,-cross_u,'k','linewidth',2)
xlim([-250 0])
ylabel('crosstrack velocity [m/s]')
text(-245,2.8,'b)','FontSize',18)
set(gca,'FontSize',13)

%ax.xtick([-250 -200 -150 -100 -50 0])
set(gca,'XTickLabel',{'0','50','100','150','200','250'})
hold on
%plot drifter release point
plot(-section_dist(183),-cross_u(183),'rx','markersize',10,'linewidth',3)
vline(-section_dist(183), 'r--') %183

axes('position',[0.1 0.1 0.6 0.2])
[haxes,hline1,hline2] = plotyy(grid_dist(1:end-1),fixgaps(filtreNaN_lowpass(vtcty,5)),grid_dist(1:end-1),fixgaps(Fjor))
ylabel('\zeta/f')
ylabel(haxes(2),'(\beta_0 - U_{yy})(U-U_s)');
xlabel('distance along transect [km]')
set(hline1,'LineWidth',2);
text(5,1.6,'c)','FontSize',18)
hold on
%plot(grid_dist(1:end-1),fixgaps(Ray),'Color', [.7 .7 .7],'linewidth',2)
hline(0,'k')
vline(35, 'r--')
set(gca,'FontSize',13)
ylim(haxes(1),[-2.5 2.5])
ylim(haxes(2),[-1 1])
%print -depsc SADCP_drifter.eps
