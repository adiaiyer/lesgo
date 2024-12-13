% Reads in lesgo's binary output data files to facilitate post-processing
%
% author: Joel Bretheim
%         Thanks to Richard Stevens for providing the basic scanning routine
%         which is embedded in the get*.m functions
%
% requires:  lesgo_param.out (in working directory)  
%            binary output files (also, user must specify which snapshot)

clear all; close all; clc;

% specify which files to load
avgVelocities    = true;
dns_profiles     = true; % Re_tau = 1000 channel flow
domain_snapshots = false;
x_snapshots      = false;
y_snapshots      = false;
z_snapshots      = false;
points           = false;

% specify file names (must choose a particular velocity snapshot)
snap_time = 150010;
xloc = 0.1;
yloc = 0.1;
zloc = 0.1;
ploc = [0.1,0.1,0.1];

% read in computational domain parameters from lesgo_param.out 
%% 
p = getParams('output/lesgo_param.out');
%% 

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% fetch average velocity fields
if avgVelocities
    [u_f,v_f,w_uv_f] = getAvgVelUV(p);
    w = getAvgVelW(p);
    [ wx,wy,wz ] = getAvgVort(p);
    [uu,vv,ww,uw,vw,uv] = getReyStress(p);
    [u2,v2,w2,uw2,vw2,uv2] = getAvgVel2(p);
    [txx,txy,tyy,txz,tyz,tzz] = getAvgTau(p);
    [fx,fy,fz] = getAvgForce(p);
    CS = getCSOpt2(p);
end

nx=numel(p.x);
ny=numel(p.y);
nz=numel(p.z_w);

u_w=zeros(nx,ny,nz+1);
v_w=zeros(nx,ny,nz+1);
u_w(:,:,2:nz)=0.5*(u(:,:,1:nz-1)+u(:,:,2:nz));
v_w(:,:,2:nz)=0.5*(v(:,:,1:nz-1)+v(:,:,2:nz));
% fetch instantaneous snapshot velocity fields
if domain_snapshots
    [ubig,vbig,wbig] = getSnap(p,snap_time);
    [theta] = getheta(p,snap_time);
    [Cs,LM,MM,QN,NN,nusgs] = getSGS(p,snap_time);
end

close all
X=p.x;Y=p.y;Z_uv=p.z_uv;
xoB=X(10);
xoE=X(23);
yoB=Y(39);
yoE=Y(91);

figure
set(gcf,'Color','w');
zloc=1;
%u2D(:,:)=ubig(:,:,zloc)-squeeze(mean(mean(ubig(:,:,zloc),1),2));
u2D(:,:)=u_f(:,:,zloc)-squeeze(mean(mean(u_p(:,:,zloc),1),2));

contourf(X,Y,u2D',22,'linestyle','none');hold on
title(sprintf('$z=%0.2f$ (m)',Z_uv(zloc)*1000),'interpreter','latex')
xlabel('$x$ (km)','interpreter','latex')
ylabel('$y$ (km)','interpreter','latex')

plot([xoB xoE],[yoB yoB],'k','linewidth',2);
plot([xoB xoE],[yoE yoE],'k','linewidth',2);
plot([xoB xoB],[yoB yoE],'k','linewidth',2);
plot([xoE xoE],[yoB yoE],'k','linewidth',2);

colormap('parula')
colorbar
set(gca,'DataAspectratio',[1 1 1]);

figure
yloc=32;
X=p.x;Y=p.y;Z_uv=p.z_uv;
%u2D_xz(:,:)=squeeze(mean(ubig(:,1:p.ny,:),2))-squeeze(mean(mean(ubig(:,:,:),1),2))';
u2D_xz(:,:)=squeeze(mean(u_f(:,1:p.ny,:),2))-squeeze(mean(mean(u_p(:,:,:),1),2))';
contourf(X,Z_uv,u2D_xz'.*(u2D_xz'<0),22,'linestyle','none')
xlabel('$x$ (km)','interpreter','latex')
ylabel('$z$ (km)','interpreter','latex')
title('$u-U(z)$','interpreter','latex')
xlim([0 30])
xticks([0:2:30])
caxis([-1 0])
ylim([0 1])
colormap('parula')
colorbar
grid on
set(gca,'DataAspectratio',[1 .1 1]);

close all
v = VideoWriter('farm_yz_no_veer.avi');  % Specify the filename and format
v.FrameRate = 30;                  % Set the frame rate (e.g., 30 frames per second)
v.Quality = 100;                   % Set the video quality (range: 0-100)
open(v);                           % Open the video file for writing
density = 10;
 swidth = 0.1;
for xloc=10:1:110
figure(1)
X=p.x;Y=p.y;Z_uv=p.z_uv;
u2D_yz(:,:)=squeeze(u_f(xloc,:,:))-squeeze(mean(mean(u_p(:,:,:),1),2))';
vp2D_yz(:,:)=squeeze(v_f(xloc,:,:))-squeeze(mean(mean(v_p(:,:,:),1),2))';
v2D_yz(:,:)=squeeze(v_f(xloc,:,:))-0*squeeze(mean(mean(v_p(:,:,:),1),2))';
w2D_yz(:,:)=squeeze(w_uv_f(xloc,:,:));
contourf(Y-Y(64),Z_uv,u2D_yz'.*(u2D_yz'<0),22,'linestyle','none');hold on
[verts, averts]=streamslice(Y-Y(64),Z_uv,vp2D_yz',w2D_yz',density,'k');
s=streamline([verts averts]);
set( s, 'Color', 'k','linewidth',swidth);
title(sprintf('$x=%0.2f$ km',X(xloc)),'interpreter','latex');
xlabel('$y$ km','interpreter','latex');
ylabel('$z$ km','interpreter','latex');
caxis([-1 0])
ylim([0 1])
xlim([-5 5])
colormap('parula')
colorbar
set(gca,'DataAspectratio',[1 .25 1]);
%pause(.05)
% Capture the current frame
currentFrame = getframe(gcf);  % Capture the current figure frame

% Write the frame to the video file
writeVideo(v, currentFrame);   % Write the frame to the video file
end

close(v);  


close all
v = VideoWriter('my_video_xy.avi');  % Specify the filename and format
v.FrameRate = 5;                  % Set the frame rate (e.g., 30 frames per second)
v.Quality = 100;                   % Set the video quality (range: 0-100)
open(v);                           % Open the video file for writing

for zloc=1:27
 figure(1)
u2D_xy(:,:)=squeeze(u_f(:,:,zloc))-squeeze(mean(mean(u_p(:,:,zloc),1),2))';
contourf(X,Y-Y(64),u2D_xy'.*(u2D_xy'<0),22,'linestyle','none');hold on

plot([xoB xoE],[yoB yoB]-5,'k','linewidth',2);
plot([xoB xoE],[yoE yoE]-5,'k','linewidth',2);
plot([xoB xoB],[yoB yoE]-5,'k','linewidth',2);
plot([xoE xoE],[yoB yoE]-5,'k','linewidth',2);
title(sprintf('$z=%0.0f$ m',Z_uv(zloc)*1000),'interpreter','latex');
xlabel('$x$ (km)','interpreter','latex');
ylabel('$y$ (km)','interpreter','latex');
caxis([-1 0])
xlim([0 30])
ylim([-5 5])
colormap('parula')
colorbar
set(gca,'DataAspectratio',[1 1 1]);
pause(.1)
% Capture the current frame
currentFrame = getframe(gcf);  % Capture the current figure frame

% Write the frame to the video file
writeVideo(v, currentFrame);   % Write the frame to the video file
end

close(v);  

close all
figure
% set(gcf,'Units','Centimeters','Color','w','Position',[10 10 30 30])
for k=1:numel(p.z_uv)
u3D(:,:,k)=(u_f(:,:,k))-squeeze(mean(mean(u_p(:,:,k),1),2));
end
up_mean(1:numel(p.z_uv))=squeeze(mean(mean(u_p(1:97,:,1:numel(p.z_uv)),1),2));
vp_mean(1:numel(p.z_uv))=squeeze(mean(mean(v_p(1:97,:,1:numel(p.z_uv)),1),2));


Ct_prime=1.33;
D=0.100;
zh=0.100;
vonk=0.4;
sx=7.85;
sy=5.2333;
zo=1e-6;
induction_factor = Ct_prime / (4 + Ct_prime);
Ct_noprime = 4*(induction_factor) * (1 - induction_factor);
cft=pi*Ct_noprime/(4.*sx*sy);
nu_w=sqrt(cft/2) * D * (1/(vonk^2)*(1/zh)*log(zh/zo));
zoCalaf = zh*(1+D/(2*zh))^(nu_w/(1+nu_w))*exp(-(cft/(2*vonk^2) ...
           +(log((zh/zo)*(1-D/(2*zh))^(nu_w/(1+nu_w))))^(-2))^(-0.5))
           
% ax=axes;
% ax.Units='Centimeters';
[Y3D,X3D,Z3D]=meshgrid(Y-Y(64),X,Z_uv);
isosurface(Y3D,X3D,Z3D,u3D.*(u3D<0),-.4);hold on
plot3(vp_mean,0*up_mean,Z_uv,'b','linewidth',2);hold on
set(gca,'DataAspectratio',[1 1 .1]);
% Set lighting and material properties
% lighting gouraud
% shading interp
% material dull
% Set the color of the isosurface
%color = [0.5, 0.8, 1]; % RGB color values (light blue)
%patch(iso, 'FaceColor', color, 'EdgeColor', 'none');view(30,30)
xlabel('$y$ (km)','interpreter','latex');
ylabel('$x$ (km)','interpreter','latex');
zlabel('$z$ (km)','interpreter','latex');
ylim([0 30])
zlim([0 .4])
%ax.Position=[0 5 30 15];
% Save the isosurface plot as a 3D model file (e.g., OBJ or STL)

% ax=axes;
% ax.Units='Centimeters';
% plot3(vp_mean*5,up_mean/3,Z_uv,'k','linewidth',2);hold on
% plot3(vp_mean*5,0*up_mean/3,Z_uv,'b','linewidth',2);hold on
% plot3(vp_mean*5*0,up_mean/3,Z_uv,'r','linewidth',2);hold on
% view(30,30)
% set(gca,'DataAspectratio',[1 1 .1]);
% ylim([-5 5])
% xlim([-5 5])
% zlim([0 .4])
% ax.Position=[5 5 5 5];
% box off






u_mean(1:Nz,counter)=squeeze(mean(mean(ubig(Xs:Xe,:,1:Nz),1),2));
v_mean(1:Nz,counter)=squeeze(mean(mean(vbig(Xs:Xe,:,1:Nz),1),2));
w_mean(1:Nz,counter)=squeeze(mean(mean(wbig(Xs:Xe,:,1:Nz),1),2));

uMean1=squeeze(mean(mean(ubig1(1:271,:,:),1),2));
uMean2=squeeze(mean(mean(ubig2(1:271,:,:),1),2));
uMean3=squeeze(mean(mean(ubig3(1:271,:,:),1),2));
uMean4=squeeze(mean(mean(ubig4(1:271,:,:),1),2));
uMean5=squeeze(mean(mean(ubig5(1:271,:,:),1),2));

Z=p.z_uv*1000;Nz=numel(Z);D=100;
semilogy(uMean1(1:Nz),Z/D,'r');hold on
semilogy(uMean2(1:Nz),Z/D,'b');hold on
semilogy(uMean3(1:Nz),Z/D,'k');hold on
semilogy(uMean4(1:Nz),Z/D,'m');hold on
semilogy(uMean5(1:Nz),Z/D,'g');hold on

uMean=squeeze(mean(mean(ubig(1:271,:,:),1),2));
vMean=squeeze(mean(mean(vbig(1:271,:,:),1),2));
wMean=squeeze(mean(mean(wbig(1:271,:,:),1),2));
Z=p.z_uv*1000;Nz=numel(Z);D=100;
semilogy(uMean(1:Nz),Z/D,'r');hold on
semilogy(vMean(1:Nz),Z/D,'b');hold on
semilogy(wMean(1:Nz),Z/D,'k');hold on


if x_snapshots
    [ux,vx,wx] = getSnapX(p,snap_time,xloc);
end
if y_snapshots
    [uy,vy,wy] = getSnapY(p,snap_time,yloc);
end
if z_snapshots
    [uz,vz,wz] = getSnapZ(p,snap_time,zloc);
end
if points
    [t,up,vp,wp] = getPoint(ploc(1),ploc(2),ploc(3));
end
if dns_profiles
    dns_data = importdata('dns_profiles.txt');
    dns_data = dns_data.data;
    dns_z  = dns_data(:,1)/1000; % z-plus -> z/h
    dns_u  = dns_data(:,2); % u-plus
    dns_uw = dns_data(:,3);
    dns_uu = dns_data(:,4);
    dns_ww = dns_data(:,5);
    dns_vv = dns_data(:,6);
    dns_tau= dns_data(:,8);
    dns_tot= dns_data(:,9);
end
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% averages across p.x and p.y directions (already time-averaged)
Xs=1;
Xe=281;
dz=p.z_uv(2)-p.z_uv(1);
dvMean=gradient(vMean(1:432));
dvMeandz=dvMean/dz;
for j=1:432
    wx2(1:360,1:128,j)=wx(1:360,1:128,j)+dvMeandz(j);
    u2(1:360,1:128,j)=u(1:360,1:128,j)-uMean(j);
    v2(1:360,1:128,j)=v(1:360,1:128,j)-vMean(j);
end
if avgVelocities
    uMean = squeeze(mean(mean(u(Xs:Xe,:,:))));
    vMean = squeeze(mean(mean(v(Xs:Xe,:,:))));
    wMean = squeeze(mean(mean(w_uv(Xs:Xe,:,:))));
    u_wMean = squeeze(mean(mean(u_w(Xs:Xe,:,:))));
    v_wMean = squeeze(mean(mean(v_w(Xs:Xe,:,:))));
    
    uuMean_p = squeeze(mean(mean(uu(Xs:Xe,:,:))));
    vvMean_p = squeeze(mean(mean(vv(Xs:Xe,:,:))));
    wwMean_p = squeeze(mean(mean(ww(Xs:Xe,:,:))));
    uwMean_p = squeeze(mean(mean(uw(Xs:Xe,:,:))));
    vwMean_p = squeeze(mean(mean(vw(Xs:Xe,:,:))));
    uvMean_p = squeeze(mean(mean(uv(Xs:Xe,:,:))));
    
    txxMean_p = squeeze(mean(mean(txx(Xs:Xe,:,:))));
    tyyMean_p = squeeze(mean(mean(tyy(Xs:Xe,:,:))));
    tzzMean_p = squeeze(mean(mean(tzz(Xs:Xe,:,:))));
    txzMean_p = squeeze(mean(mean(txz(Xs:Xe,:,:))));
    tyzMean_p = squeeze(mean(mean(tyz(Xs:Xe,:,:))));
    txyMean_p = squeeze(mean(mean(txy(Xs:Xe,:,:))));
    
    uu2Mean = squeeze(mean(mean(u2(Xs:Xe,:,:))))-uMean.^2;
    vv2Mean = squeeze(mean(mean(v2(Xs:Xe,:,:))))-vMean.^2;
    ww2Mean = squeeze(mean(mean(w2(Xs:Xe,:,:))))-wMean.^2;
    uw2Mean = squeeze(mean(mean(uw2(Xs:Xe,:,:))))-u_wMean.*wMean;
    vw2Mean = squeeze(mean(mean(vw2(Xs:Xe,:,:))))-v_wMean.*wMean;
    uv2Mean = squeeze(mean(mean(uv2(Xs:Xe,:,:))))-u_wMean.*v_wMean;

    
    
    fxMean = squeeze(mean(mean(fx(Xs:Xe,:,:))));
    txzMean = squeeze(mean(mean(txz(Xs:Xe,:,:))));    
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% basic plots
if avgVelocities
    figure
    plot(p.z_uv,uMean,'b')
    if dns_profiles
        hold on
        plot(dns_z,dns_u,'r')
    end
    ylim([0,30])
    xlabel('z','interpreter','tex')
    ylabel('<u>','interpreter','tex')

    figure
    kappa = 0.41;  z0 = .0001186;zh=0.1;
    loglaw = 1/kappa * log(p.z_uv ./ z0);    % rough wall
    semilogx(p.z_uv/zh,loglaw,'k')
    hold on
    semilogx(p.z_uv,uMean,'ob')
    if dns_profiles
      semilogx(dns_z,dns_u,'r')
    end
    hold off
    xlim([0.01,1])
    ylim([0,25])
    xlabel('z','interpreter','tex')
    ylabel('<u>','interpreter','tex')
    if dns_profiles
        legend('Log Law','LES','DNS','Location','best')
    else
        legend('Log Law','LES','Location','best')
    end

    
    fxMeanw=zeros(nz,1);
    fxMeanw(2:nz,1)=0.5*(fxMean(2:nz)+fxMean(1:nz-1));
    fxw=trapz(p.z_w,fxMeanw)-cumtrapz(p.z_w,fxMeanw);
    
    figure
    plot( -uw2Mean(1:nz),p.z_w,'g','linewidth',1)
    hold on
    plot( -txzMean(1:nz),p.z_w,'b','linewidth',1)
    plot( -txzMean(1:nz) - uw2Mean(1:nz) + fxw,p.z_w,'k','linewidth',1)
    plot(- fxw,p.z_w,'m','linewidth',1)
    
    plot( -uwMean(1:nz),p.z_w,'k--','linewidth',1)
    plot( -(uw2Mean(1:nz)-uwMean(1:nz)),p.z_w,'r--','linewidth',1)
    
    
    
    if dns_profiles
        plot(dns_z,-dns_uw,'b')
        plot(dns_z,dns_tau,'c')
        plot(dns_z,dns_tau-dns_uw,'k')
    end
    plot(p.z_w, (1-p.z_w),'k')
    hold off
    xlabel('z','interpreter','tex')
    ylabel('<u''w''>','interpreter','tex')
    
    figure
    plot(p.z_uv, uuMean, 'ob')
    if dns_profiles
        hold on
        plot(dns_z,dns_uu,'b')
    end
    ylim([0,8])
    xlabel('z','interpreter','tex')
    ylabel('<u''u''>','interpreter','tex')
    
    figure
    plot(p.z_uv, vvMean, 'ob')
    if dns_profiles
        hold on
        plot(dns_z,dns_vv, 'b')
    end
    ylim([0,4])
    xlabel('z','interpreter','tex')
    ylabel('<v''v''>','interpreter','tex')
    
    figure
    plot(p.z_w, wwMean, 'ob')
    if dns_profiles
        hold on
        plot(dns_z,dns_ww,'b')
    end
    ylim([0,2])
    xlabel('z','interpreter','tex')
    ylabel('<w''w''>','interpreter','tex')
end

if domain_snapshots
    figure
    [X,Y] = meshgrid(p.x,p.y);
    pcolor(X,Y,ubig(:,:,4)')
    xlabel('x')
    ylabel('y')
    title(['Streamwise velocity at z = ',num2str(p.z_uv(4))])
    shading interp; colorbar;
end
    
if x_snapshots
    figure
    [Y,Z] = meshgrid(p.y,p.z_uv);
    pcolor(Y,Z,ux')
    xlabel('y')
    ylabel('z')
    title(['Streamwise velocity at x = ',num2str(xloc)])
    shading interp; colorbar;
end

if y_snapshots
    figure
    [X,Z] = meshgrid(p.x,p.z_uv);
    pcolor(X,Z,uy')
    xlabel('x')
    ylabel('z')
    title(['Streamwise velocity at y = ',num2str(yloc)])
    shading interp; colorbar;
end

if z_snapshots
    figure
    [X,Y] = meshgrid(p.x,p.y);
    pcolor(X,Y,uz')
    xlabel('x')
    ylabel('y')
    title(['Streamwise velocity at z = ',num2str(zloc)])
    shading interp; colorbar;
end

if points
    figure
    plot(t,[up,vp,wp])
    xlabel('t')
    ylabel('Velocity')
    legend('u','v','w')
end


%%
y1=[0 1.5];
x1=[1*3.75-.125*3.75 1*3.75-.125*3.75];
plot(x1,y1,'k','linewidth',.5);hold on
x1=[.75*3.75-.125*3.75*0 .75*3.75-0*.125*3.75];
plot(x1,y1,'k','linewidth',.5);hold on
x1=[.75*3.75-.125*3.75*1 .75*3.75-1*.125*3.75];
plot(x1,y1,'k','linewidth',.5);hold on

%turbine line coordinate
ptop=[0,.05];
pbot=[0,-.05];

%flow angle
ptop=[0-3.75/2,0.75-1.5/2];
pbot=[3.75-3.75/2,.75-1.5/2];

%local pert lines
ptop=[0,5];
pbot=[0,-5];


theta=20;
A=[cosd(theta) -sind(theta);sind(theta) cosd(theta)];
ptopr=mtimes(A,ptop');
pbotr=mtimes(A,pbot');

%turbine
xt2=[pbotr(1) ptopr(1)]+.5*0;
yt2=[pbotr(2) ptopr(2)]+.75*0;

%flow angle
xt=[pbotr(1) ptopr(1)]+3.75/2;
yt=[pbotr(2) ptopr(2)]+1.5/2;

%local pert lines
xt=[pbotr(1) ptopr(1)]+24;
yt=[pbotr(2) ptopr(2)]+12;

plot(xt2,yt2,'k','linewidth',.5)

%SGS
Cs_opt2_4d=max(0,QN./NN);
Cs_opt2_2d=max(0,LM./MM);
tf1 = 2;
tf2 = 4;
Beta=(Cs_opt2_4d./Cs_opt2_2d).^(log(tf1)/(log(tf2)-log(tf1)));
Betaclip=max(Beta,1/tf1/tf2);
Cs_opt2=Cs_opt2_2d./Betaclip;

nut_f_smag(:)=mean(mean(nusgs_f_smag(:,:,:),1),2);
cs_f_smag(:)=mean(mean(Cs_f_smag(:,:,:),1),2);
cs2(:)=mean(mean(Cs_opt2(:,:,:),1),2);

lm(:)=mean(mean(LM(:,:,:),1),2);
mm(:)=mean(mean(MM(:,:,:),1),2);
qn(:)=mean(mean(QN(:,:,:),1),2);
nn(:)=mean(mean(NN(:,:,:),1),2);


hold on
plot(p.z_uv,nut(1:432))
plot(p.z_uv,cs(1:432))
plot(p.z_uv,lm(1:432))
plot(p.z_uv,mm(1:432))
plot(p.z_uv,qn(1:432))
plot(p.z_uv,nn(1:432))

num=lm./mm;
den=max((qn./nn)./(lm./mm),.125);
plot(p.z_uv,num(1:432)./den(1:432))

dx=p.x(2)-p.x(1);
dy=p.y(2)-p.y(1);
dz=p.z_uv(2)-p.z_uv(1);
delta=(dx*dy*dz)^(1/3);
Co=0.16;
kappa=0.41;
l=(Co^2*(kappa*p.z_uv).^(-2)+delta^(-2)).^(-1/2);
Csmag=Co^2*l.^2/delta^2;

%Vorticity contour
figure
D=0.1;%km
xloc=68;
wx2D(:,:)=wx(xloc,52:78,1:66);
contourf(p.y(52:78)/D-7.5,p.z_uv(1:66)/D-1,wx2D',22,'linestyle','none')
set(gca,'DataAspectratio',[1 1 1]);
set(gca,'fontsize',16,'fontname','times');
set(gcf,'Color','w','units','centimeters');
xlabel('$y/D$','fontsize',16,'interpreter','latex')
%ylabel('$z/D$','fontsize',16,'interpreter','latex')
yticks([])
title(['$x/D$ = ',num2str(round((p.x(xloc)-0.5)/D))],'interpreter','latex')
colormap(gca,redblue)
caxis([-30 30])
colorbar
ax=gca;
ax.Units='centimeters';
circle(0,0,.05/D)
ax.Position=[2+8.5*(0) 4 6 6];


%% Otsu segmentation
filename='wx_otsu_2.gif';
for i=[49,59,69,79,89,99,109,119,129,139,149,159,168,178,188]
    xloc=i;
    wx2D=zeros(27,66);
    wx2D_otsu(:,:)=wx_otsu(i,:,:);
    %wx2D(:,5:66)=wx(i,52:78,5:66)*(D/ub);%.*abs(wx_otsu(i,:,:));
    wx2D(:,5:66)=wx(i,52:78,5:66)*(D/ub).*abs(wx_otsu(i,:,:));

    contourf(p.y(52:78)/D-7.5,p.z_uv(1:66)/D-1,wx2D',50,'linestyle','none');hold on
    contour(p.y(52:78)/D-7.5,p.z_uv(5:66)/D-1,abs(wx2D_otsu'),[1],'linestyle','-','Color','k','linewidth',1);
    caxis([min(wx2D(:)) max(wx2D(:))])
    %caxis([-10 10]);
    xlim([-1.5 1.5]);
    ylim([-1 2]);
    
    %caxis([-30 30])
    set(gca,'DataAspectratio',[1 1 1]);
    set(gca,'fontsize',16,'fontname','times');
    set(gcf,'Color','w','units','centimeters');
    xlabel('$y/D$','fontsize',16,'interpreter','latex')
    ylabel('$z/D$','fontsize',16,'interpreter','latex')
    %yticks([0:1:2])
    title(sprintf('$x/D$ =%d ',round((p.x(xloc)-0.5)/D )),'interpreter','latex')
    colormap(gca,bluewhitered(50))
    cba=colorbar;
    t=get(cba,'Limits');
    set(cba,'Ticks',[t(1),t(1)/2,0,t(2)/2,t(2)]);
    cba.Ruler.TickLabelFormat = '%.2f';
    ax1=gca;
    ax1.Units='centimeters';
    circle(0,0,.05/D)
    set(gcf,'Units','Centimeters','pos',[5 5 11 10]);
    ax1.Position=[2 2 6 6];
    if i==49
       gif(filename,'DelayTime',1.5,'LoopCount',5,'frame',gcf)
    else
       gif('frame',gcf)
    end
    pause(1)
end

%% wx_max and circulation
counter=1;
for i=49:193
 
    wx_p(:,:)=wxp(i,:,1:432).*wx_otsu_p(i,:,1:432);
    wx_n(:,:)=wxp(i,:,1:432).*wx_otsu_n(i,:,1:432);

    wxmax_p(counter)=max(max(wx_p));
    wxmax_n(counter)=max(max(-wx_n));
    
    alpha_p(counter) = thresh_p(i)/wxmax_p(counter);
    alpha_n(counter) = thresh_n(i)/wxmax_n(counter);
    
    Gamma_p(counter)=sum(sum(wx_p))*dy*dz/(1-alpha_p(counter));
    Gamma_n(counter)=sum(sum(wx_n))*dy*dz/(1-alpha_n(counter));
    
    xmax(counter)=(X(i)-0.5)/D;
    counter=counter+1;
    
end

plot(xmax,wxmax_p*D/ub,'r','linewidth',1);hold on
%plot(xmax,wmax_p*D/ub,'r--');
plot(x_wxmax_analytical, wxmax_Shapiro_analytical,'r','linewidth',1);
plot(xmax,-wxmax_n*D/ub,'b','linewidth',1);
%plot(xmax,wmax_n*D/ub,'b--');
plot(xn_Carl, wxn_Carl,'bv');
h=legend('$\omega^+_{x,max}$','$\omega^+_{x,max}$ \ (Shapiro et al, 2020)');
set(h,'interpreter','latex','fontsize',13,'location','Northeast');

plot(xmax,Gamma_p/(D*ub),'r','linewidth',1);hold on
plot(x_circ_analytical, circ_analytical,'r','linewidth',2)
plot(xmax,Gamma_n/(D*ub),'b','linewidth',1);hold on
plot(x_circ_analytical, -circ_analytical,'b','linewidth',2)

%%
 h=legend('$\Gamma^+$ Shapiro et al (2020)','$\Gamma^-$ Shapiro et al (2020)',...
     '$\Gamma^+$ (Const Smag)','$\Gamma^-$ (Const Smag)',...
     '$\Gamma^{\prime +}$ (no wind veer, Lag scale-dep)','$\Gamma^{\prime -}$ (no wind veer, Lag scale-dep)',...
     '$\Gamma^{ +}$ (Lag scale-dep)','$\Gamma^{-}$ (Lag scale-dep)');
set(h,'interpreter','latex','fontsize',13,'location','Northeast Outside');

for i=1:359
    wxmax(i)=max(max(sqrt(wx(i,52:78,1:66).^2+wy(i,52:78,1:66).^2+wz(i,52:78,1:66).^2)));
end
loglog(p.x,wxmax)
xlim([0.5 3]);
set(gca,'fontsize',16,'fontname','times');
set(gcf,'Color','w','units','centimeters');
ax=gca;
ax.Units='centimeters';
ax.Position=[2 4 8 4];

%% mean plot
ax1 = axes('XAxisLocation','bot',...
    'YAxisLocation','left',...
    'Color','none','Units','Centimeters');
line(uMean_p(1:432)/ub,(Z),'Color','r','linewidth',1);
xlim([0 1.5]);
ylim([0 2]);
xticks([0:.25:2]);
ax1.Position=[2 2 6 6];
ax1.XColor = 'r';

set(ax1,'fontsize',16,'fontname','times');
xlabel('$\langle U\rangle/U_h$','fontsize',16,'interpreter','latex');
ylabel('$z$','fontsize',16,'interpreter','latex');
set(gcf,'Color','w','units','centimeters');

ax2 = axes('XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none','Units','Centimeters','YTick',[]);
ax2.Position=[2 2 6 6];
ax2.XColor='b';
line(vMean_p(1:432)/ub,(Z),'Color','b','linewidth',1);
xlim([-.4 .1]);
xticks([-0.4:.1:.1]);
ylim([0 2]);
set(ax2,'fontsize',16,'fontname','times');
xlabel('$\langle V\rangle/U_h$','fontsize',16,'interpreter','latex');
ylabel('$z$','fontsize',16,'interpreter','latex');


ax3=axes;
line(uMean_p(1:432)/ub,(Z-zt)/D,'Color','r','linewidth',1);
ax3.XColor = 'r';
ax3.YColor = 'k';
ax3.Units='centimeters';
ax3.Position=[10 2 6 6];
xlabel('$\langle U\rangle/U_h$','fontsize',16,'interpreter','latex');
ylabel('$(z-z_h)/D$','fontsize',16,'interpreter','latex');
set(ax3,'fontsize',16,'fontname','times');
ylim([-1 2]);
xlim([0.2 1.25]);
xticks([0:.25:1.5]);

ax4 = axes('XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none','Units','Centimeters');
ax4.Position=ax3.Position;
ax4.YTick=[];
ax4.XColor = 'b';
line(vMean_p(1:432)/ub,(Z-zt)/D,'Parent',ax4,'Color','b','linewidth',1);
ax4.YLim=[-1 2];
ax4.XLim=[-.125 0.05];
set(ax4,'fontsize',16,'fontname','times');
xlabel('$\langle V\rangle/U_h$','fontsize',16,'interpreter','latex');
ax4.YTick=[];
ax4.XTick=[-.1:.05:0.05];