%LES parameters
Lx  = 2000;%m
Ly  = 2000;%m
Lz  = 1000;%m
Hbar = Lz;
z_i = 1000;%m
z0  = 2e-4;%m
Nx  = 256;
Ny  = 256;
Nz  = 192;
ak  = 0.1;
ustar_LESGO = 1;%m/s
c_by_ustar = 7;
wave_angle = 0;%rad
Ug  = 8;%m/s
Vg  = 0;
G   = 8;%m/s
fc  = 1e-4;%1/s
Cr  = -1/3600;%K/s
nu  = 0;%Molecular viscosity


g_accl         = 9.81;%m/s^2
c_phase        = c_by_ustar*ustar_LESGO;
cx_phase       = c_phase*cos(wave_angle);
cy_phase       = c_phase*sin(wave_angle);
k_wavno        = g_accl/c_phase^2;
kx_wavno       = k_wavno*cos(wave_angle);
ky_wavno       = k_wavno*sin(wave_angle);
a_amp          = ak/k_wavno;
omega_freq     = c_phase*k_wavno;

dx = Lx/Nx;
dy = Ly/Ny;
dz = Lz/Nz;
x=zeros(Nx+1,1);
y=zeros(Ny+1,1);

for j = 1:Ny+1
    y(j) = (j-1)*dy;
end

for i = 1: Nx+1
    x(i) = (i - 1)*dx;
end

z=zeros(1,Nz+1);
zw=zeros(1,Nz+1);
for k=1:Nz+1
    z(1,k)=(k-0.5)*dz;
end
zw(1,:) = z-dz/2;


kx  = [0:Nx/2, -Nx/2+1:-1] * (2 * pi / Lx);
ky  = [0:Ny/2, -Ny/2+1:-1] * (2 * pi / Ly);
kx_3D=zeros(Nx,Ny,Nz);
for j=1:Ny
    for k=1:Nz
        kx_3D(:,j,k)=kx;
    end
end
ky_3D=zeros(Nx,Ny,Nz);
for i=1:Nx
    for k=1:Nz
        ky_3D(i,:,k)=ky;
    end
end


%Initialize phase average quantities
  u_phase_avg    = zeros(Nx, Nz);
u_w_phase_avg    = zeros(Nx, Nz);
  w_phase_avg    = zeros(Nx, Nz);
  p_phase_avg    = zeros(Nx, Nz);
p_w_phase_avg    = zeros(Nx, Nz);
  W_phase_avg    = zeros(Nx, Nz);
  U_phase_avg    = zeros(Nx, Nz);

 uW_phase_avg    = zeros(Nx, Nz);
 uU_phase_avg    = zeros(Nx, Nz);
 wU_phase_avg    = zeros(Nx, Nz);
 wW_phase_avg    = zeros(Nx, Nz);


eta_phase_avg    = zeros(Nx, 1 ); 
taup_13_phase_avg= zeros(Nx, Nz);
taup_11_phase_avg= zeros(Nx, Nz);

  tau13_phase_avg= zeros(Nx, Nz);
  tau31_phase_avg= zeros(Nx, Nz);
  tau11_phase_avg= zeros(Nx, Nz);
  tau33_phase_avg= zeros(Nx, Nz);

tau13_SGS_phase_avg= zeros(Nx, Nz);
tau31_SGS_phase_avg= zeros(Nx, Nz);
tau11_SGS_phase_avg= zeros(Nx, Nz);
tau33_SGS_phase_avg= zeros(Nx, Nz);

 tau13_nu_phase_avg= zeros(Nx, Nz);
 tau31_nu_phase_avg= zeros(Nx, Nz);
 tau11_nu_phase_avg= zeros(Nx, Nz);
 tau33_nu_phase_avg= zeros(Nx, Nz);

 tau13_wave_phase_avg= zeros(Nx, Nz);
 tau31_wave_phase_avg= zeros(Nx, Nz);
 tau11_wave_phase_avg= zeros(Nx, Nz);
 tau33_wave_phase_avg= zeros(Nx, Nz);

n_array = 725000%:500:925000;

Nt = numel(n_array)
addpath('/Users/ghaneshnarasimhan/Library/Mobile Documents/com~apple~CloudDocs/Desktop/JHU/LESGO/lesgo/post-processing/matlab-scripts/');

ke_data=importdata('./output/check_ke.dat');
time_array=ke_data(:,1);
p = getParams('./output/lesgo_param.out');

counter = 1;
for n = n_array
   tic
   counter

   dzw_3D=dz*ones(Nx,Ny,Nz);
   dz_3D=dzw_3D;

   time = time_array(n);
   
   eta  = zeros(Nx,Ny);
   for i=1:Nx
       for j=1:Ny
        eta(i,j)  = a_amp*cos( kx_wavno*x(i) + ky_wavno*y(j) - omega_freq*time );
       end
   end
   
   p_uv         = get_pres_snap(p,n);
   [u,v,w_uv] = getSnap(p,n);
  
   z_phys = zeros(Nx,Nz);
   x_phys = zeros(Nx,Nz);
   for i=1:Nx
        z_phys(i,:) = eta(i,1)+(zw(1,1:Nz)/Hbar).*(Hbar-eta(i,1));%zw(:).*(Hbar-eta_phase_avg(i))+eta_phase_avg(i);
        x_phys(i,:) = x(i);
   end

   %zeta_x, zeta_z - Hao & Shen (2022)
   zetax_w = zeros(Nx,Ny,Nz);
   zetay_w = zeros(Nx,Ny,Nz);
   zetaz_w = zeros(Nx,Ny,Nz);
  
   eta_x = zeros(Nx,Ny);
   for j=1:Ny
     eta_hat_1       = fft(eta(1:Nx,j))/(Nx);
     eta_x(1:Nx,j)   = real(ifft(sqrt(-1)*kx'.*eta_hat_1))*(Nx);
   end
  
   eta_y = zeros(Nx,Ny);
   for i=1:Nx
     eta_hat_2         = fft(eta(i,1:Ny))/(Ny);
     eta_y(i,1:Ny)   = real(ifft(sqrt(-1)*ky.*eta_hat_2))*(Ny);
   end
   eta_y(:,Ny) = eta_y(:,1);

   for k=1:Nz
     zetax_w(:,:,k) = eta_x(:,:).*(zw(k)-Hbar)./(Hbar-eta(:,:)) ;%( eta_x(:,:) - zw(k) * eta_x(:,:) )./(eta(:,:) + Hbar);
     zetay_w(:,:,k) = eta_y(:,:).*(zw(k)-Hbar)./(Hbar-eta(:,:)) ;%( eta_y(:,:) - zw(k) * eta_y(:,:) )./(eta(:,:) + Hbar);
     zetaz_w(:,:,k) = Hbar./(Hbar-eta(:,:));%1./(eta(:,:) + Hbar);
   end
 
   %interpolate u & p to w-grid
   u_w  = zeros(Nx,Ny,Nz);
   v_w  = zeros(Nx,Ny,Nz);
   p_w  = zeros(Nx,Ny,Nz);
   w  = zeros(Nx,Ny,Nz);
   
   %check the values at the bottom boundary
   u_w(:,:,  1   )  = u(:,:,1);
   u_w(:,:, Nz-1 )  = u(:,:,Nz);
   u_w(:,:, Nz )    = u(:,:,Nz);
   u_w(:,:,2:Nz-2)  = 0.5*(u(:,:,2:Nz-2) + u(:,:,3:Nz-1));

   v_w(:,:,  1   )  = v(:,:,1);
   v_w(:,:, Nz-1 )  = v(:,:,Nz);
   v_w(:,:, Nz   )  = v(:,:,Nz);
   v_w(:,:,2:Nz-2)  = 0.5*(v(:,:,2:Nz-2) + v(:,:,3:Nz-1));

   w(:,:,  1   )  = w_uv(:,:,1);
   w(:,:, Nz-1 )  = w_uv(:,:,Nz);
   w(:,:, Nz )    = w_uv(:,:,Nz);
   w(:,:,2:Nz-2)  = 0.5*(w_uv(:,:,2:Nz-2) + w_uv(:,:,3:Nz-1));

   p_w(:,:,  1   )  = p_uv(:,:,1);
   p_w(:,:, Nz-1 )  = p_uv(:,:,Nz);
   p_w(:,:, Nz   )  = p_uv(:,:,Nz);
   p_w(:,:,2:Nz-2)  = 0.5*(p_uv(:,:,2:Nz-2) + p_uv(:,:,3:Nz-1));
   
   
   %Calc Sij
   Sij=zeros(Nx,Ny,Nz,9);
   %Sij=calc_Sij(x,y,z_phys,u_w,v_w,w,Lx,Ly,Nx,Ny,Nz);
   Sij=calc_Sij_3D(kx_3D,ky_3D,dzw_3D,zetaz_w,u_w(1:Nx,1:Ny,1:Nz),v_w(1:Nx,1:Ny,1:Nz),w(1:Nx,1:Ny,1:Nz),Nx,Ny,Nz,zw(1,1:Nz));
   
  
   %tau_SGS
   [ txx,tyy,tzz,txy,txz,tyz ] = get_tau_SGS(p,n);
   tau13_SGS=zeros(Nx,Ny,Nz);
   tau31_SGS=zeros(Nx,Ny,Nz);
   tau11_SGS=zeros(Nx,Ny,Nz);
   tau33_SGS=zeros(Nx,Ny,Nz);

   
   tau13_SGS(:,:,:)=txx(1:Nx,1:Ny,1:Nz).*zetax_w(:,:,:)./zetaz_w(:,:,:)  + txz(1:Nx,1:Ny,1:Nz);
   tau31_SGS(:,:,:)=txz(1:Nx,1:Ny,1:Nz)./zetaz_w(:,:,:);
   tau11_SGS(:,:,:)=txx(1:Nx,1:Ny,1:Nz)./zetaz_w(:,:,:);
   tau33_SGS(:,:,:)=txz(1:Nx,1:Ny,1:Nz).*zetax_w(:,:,:)./zetaz_w(:,:,:)  + tzz(1:Nx,1:Ny,1:Nz);
   

   tau13_nu=zeros(Nx,Ny,Nz);
   tau31_nu=zeros(Nx,Ny,Nz);
   tau11_nu=zeros(Nx,Ny,Nz);
   tau33_nu=zeros(Nx,Ny,Nz);


   tau13_nu(:,:,:)=-2*nu*Sij(:,:,:,1).*zetax_w(:,:,:)./zetaz_w(:,:,:)  - 2*nu*Sij(:,:,:,3);
   tau31_nu(:,:,:)=-2*nu*Sij(:,:,:,7)./zetaz_w(:,:,:);
   tau11_nu(:,:,:)=-2*nu*Sij(:,:,:,1)./zetaz_w(:,:,:);
   tau33_nu(:,:,:)=-2*nu*Sij(:,:,:,7).*zetax_w(:,:,:)./zetaz_w(:,:,:)  - 2*nu*Sij(:,:,:,9);
   

   
   %Calculate contravariant velocity
          W = zeros(Nx,Ny,Nz);
          U = zeros(Nx,Ny,Nz);
    taup_13 = zeros(Nx,Ny,Nz);
    taup_11 = zeros(Nx,Ny,Nz);
    %taup_31 = 0;
    %taup_33 = same as pressure p;
   U(:,:,:) = (u_w(:,:,:) - c_phase)./zetaz_w(:,:,:);
   W(:,:,:) = (u_w(:,:,:) - c_phase).*(zetax_w(:,:,:)./zetaz_w(:,:,:)) + w(:,:,:);
   taup_13(:,:,:) =  -p_w(:,:,:).*zetax_w(:,:,:)./zetaz_w(:,:,:);
   taup_11(:,:,:) =  -p_w(:,:,:)./zetaz_w(:,:,:);
   %disp([size(W),max(W(:)), min(W(:))])



   %average in y-direction
   u_xz   = zeros(Nx,Nz);
   u_w_xz = zeros(Nx,Nz);
   w_xz   = zeros(Nx,Nz);
   p_xz   = zeros(Nx,Nz);
   p_w_xz = zeros(Nx,Nz);  
   W_xz   = zeros(Nx,Nz);
   U_xz   = zeros(Nx,Nz);

   uW_xz   = zeros(Nx,Nz);
   uU_xz   = zeros(Nx,Nz);
   wU_xz   = zeros(Nx,Nz);
   wW_xz   = zeros(Nx,Nz);

   taup_13_xz = zeros(Nx,Nz);
   taup_11_xz = zeros(Nx,Nz);

   tau13_SGS_xz = zeros(Nx,Nz);
   tau31_SGS_xz = zeros(Nx,Nz);
   tau11_SGS_xz = zeros(Nx,Nz);
   tau33_SGS_xz = zeros(Nx,Nz);

    tau13_nu_xz = zeros(Nx,Nz);
    tau31_nu_xz = zeros(Nx,Nz);
    tau11_nu_xz = zeros(Nx,Nz);
    tau33_nu_xz = zeros(Nx,Nz);

   u_xz(:,:)   = squeeze(mean(  u(1:Nx,1:Ny,1:Nz),2));
   u_w_xz(:,:) = squeeze(mean(u_w(1:Nx,1:Ny,1:Nz),2));

   w_xz(:,:)   = squeeze(mean(  w(1:Nx,1:Ny,1:Nz),2));

   p_xz(:,:)   = squeeze(mean( p_uv(1:Nx,1:Ny,1:Nz),2));
   p_w_xz(:,:) = squeeze(mean(p_w(1:Nx,1:Ny,1:Nz),2));  
   
   W_xz(:,:)   = squeeze(mean(W(1:Nx,1:Ny,1:Nz),2));
   U_xz(:,:)   = squeeze(mean(U(1:Nx,1:Ny,1:Nz),2)); 
   taup_13_xz(:,:) = squeeze(mean(taup_13(1:Nx,1:Ny,1:Nz),2));
   taup_11_xz(:,:) = squeeze(mean(taup_11(1:Nx,1:Ny,1:Nz),2));

        uW_xz(:,:) = squeeze(mean(u_w(1:Nx,1:Ny,1:Nz).*W(1:Nx,1:Ny,1:Nz),2));
        uU_xz(:,:) = squeeze(mean(u_w(1:Nx,1:Ny,1:Nz).*U(1:Nx,1:Ny,1:Nz),2));
        wW_xz(:,:) = squeeze(mean(w(1:Nx,1:Ny,1:Nz).*W(1:Nx,1:Ny,1:Nz),2));
        wU_xz(:,:) = squeeze(mean(w(1:Nx,1:Ny,1:Nz).*U(1:Nx,1:Ny,1:Nz),2));
        

   tau13_SGS_xz(:,:) = squeeze(mean(tau13_SGS(1:Nx,1:Ny,1:Nz),2));
   tau31_SGS_xz(:,:) = squeeze(mean(tau31_SGS(1:Nx,1:Ny,1:Nz),2));
   tau11_SGS_xz(:,:) = squeeze(mean(tau11_SGS(1:Nx,1:Ny,1:Nz),2));
   tau33_SGS_xz(:,:) = squeeze(mean(tau33_SGS(1:Nx,1:Ny,1:Nz),2));

    tau13_nu_xz(:,:) = squeeze(mean(tau13_nu(1:Nx,1:Ny,1:Nz),2));
    tau31_nu_xz(:,:) = squeeze(mean(tau31_nu(1:Nx,1:Ny,1:Nz),2));
    tau11_nu_xz(:,:) = squeeze(mean(tau11_nu(1:Nx,1:Ny,1:Nz),2));
    tau33_nu_xz(:,:) = squeeze(mean(tau33_nu(1:Nx,1:Ny,1:Nz),2));

   %disp('reached')
   %perform phase shift in x-direction
     u_shifted = zeros(Nz, Nx);
   u_w_shifted = zeros(Nz, Nx);
     w_shifted = zeros(Nz, Nx);
     p_shifted = zeros(Nz, Nx);
   p_w_shifted = zeros(Nz, Nx);

     W_shifted = zeros(Nz, Nx);
     U_shifted = zeros(Nz, Nx);
   eta_shifted = zeros( 1, Nx);
   taup_13_shifted = zeros(Nz, Nx);
   taup_11_shifted = zeros(Nz, Nx);
   
     
   uW_shifted = zeros(Nz, Nx);
   uU_shifted = zeros(Nz, Nx);
   wU_shifted = zeros(Nz, Nx);
   wW_shifted = zeros(Nz, Nx);

   tau13_SGS_shifted = zeros(Nz, Nx);
   tau31_SGS_shifted = zeros(Nz, Nx);
   tau11_SGS_shifted = zeros(Nz, Nx);
   tau33_SGS_shifted = zeros(Nz, Nx);

    tau13_nu_shifted = zeros(Nz, Nx);
    tau31_nu_shifted = zeros(Nz, Nx);
    tau11_nu_shifted = zeros(Nz, Nx);
    tau33_nu_shifted = zeros(Nz, Nx);

   for k = 1:Nz
            u_shifted(k, :) = phase_avg_shift(  u_xz(1:Nx, k)', kx, Nx, time, c_phase);
          u_w_shifted(k, :) = phase_avg_shift(u_w_xz(1:Nx, k)', kx, Nx, time, c_phase);
            w_shifted(k, :) = phase_avg_shift(  w_xz(1:Nx, k)', kx, Nx, time, c_phase);
            p_shifted(k, :) = phase_avg_shift(  p_xz(1:Nx, k)', kx, Nx, time, c_phase);
          p_w_shifted(k, :) = phase_avg_shift(p_w_xz(1:Nx, k)', kx, Nx, time, c_phase);
            W_shifted(k, :) = phase_avg_shift(  W_xz(1:Nx, k)', kx, Nx, time, c_phase); 
            U_shifted(k, :) = phase_avg_shift(  U_xz(1:Nx, k)', kx, Nx, time, c_phase);

      taup_13_shifted(k, :) = phase_avg_shift( taup_13_xz(1:Nx, k)', kx, Nx, time, c_phase);
      taup_11_shifted(k, :) = phase_avg_shift( taup_11_xz(1:Nx, k)', kx, Nx, time, c_phase);

        uW_shifted(k, :) = phase_avg_shift( uW_xz(1:Nx, k)', kx, Nx, time, c_phase);
        uU_shifted(k, :) = phase_avg_shift( uU_xz(1:Nx, k)', kx, Nx, time, c_phase);
        wW_shifted(k, :) = phase_avg_shift( wW_xz(1:Nx, k)', kx, Nx, time, c_phase);
        wU_shifted(k, :) = phase_avg_shift( wU_xz(1:Nx, k)', kx, Nx, time, c_phase);

      tau13_SGS_shifted(k, :) = phase_avg_shift( tau13_SGS_xz(1:Nx, k)', kx, Nx, time, c_phase);
      tau31_SGS_shifted(k, :) = phase_avg_shift( tau31_SGS_xz(1:Nx, k)', kx, Nx, time, c_phase);
      tau11_SGS_shifted(k, :) = phase_avg_shift( tau11_SGS_xz(1:Nx, k)', kx, Nx, time, c_phase);
      tau33_SGS_shifted(k, :) = phase_avg_shift( tau33_SGS_xz(1:Nx, k)', kx, Nx, time, c_phase);


       tau13_nu_shifted(k, :) = phase_avg_shift(  tau13_nu_xz(1:Nx, k)', kx, Nx, time, c_phase);
       tau31_nu_shifted(k, :) = phase_avg_shift( tau31_nu_xz(1:Nx, k)', kx, Nx, time, c_phase);
       tau11_nu_shifted(k, :) = phase_avg_shift( tau11_nu_xz(1:Nx, k)', kx, Nx, time, c_phase);
       tau33_nu_shifted(k, :) = phase_avg_shift( tau33_nu_xz(1:Nx, k)', kx, Nx, time, c_phase);
       
   end 
   %storing surface pressure as a function of time
   p_shifted_surf = zeros(Nx,1);
   p_w_shifted_surf = zeros(Nx,1);
   p_shifted_surf(:,1) = p_shifted(1,:)-mean(p_shifted(1,:));
   p_w_shifted_surf(:,1) = p_w_shifted(1,:)-mean(p_w_shifted(1,:));

   %phase average
                 eta_shifted = phase_avg_shift(squeeze(mean(eta(1:Nx, :),2))', kx, Nx, time, c_phase);
       eta_phase_avg(1:Nx,:) = eta_phase_avg(1:Nx,:) + (1 / Nt) * eta_shifted';
         u_phase_avg(1:Nx,:) =   u_phase_avg(1:Nx,:) + (1 / Nt) *   u_shifted';
       u_w_phase_avg(1:Nx,:) = u_w_phase_avg(1:Nx,:) + (1 / Nt) * u_w_shifted';
         w_phase_avg(1:Nx,:) =   w_phase_avg(1:Nx,:) + (1 / Nt) *   w_shifted';
         p_phase_avg(1:Nx,:) =   p_phase_avg(1:Nx,:) + (1 / Nt) *   p_shifted';
       p_w_phase_avg(1:Nx,:) = p_w_phase_avg(1:Nx,:) + (1 / Nt) * p_w_shifted'; 
         W_phase_avg(1:Nx,:) =   W_phase_avg(1:Nx,:) + (1 / Nt) *   W_shifted';
         U_phase_avg(1:Nx,:) =   U_phase_avg(1:Nx,:) + (1 / Nt) *   U_shifted';
   taup_13_phase_avg(1:Nx,:) = taup_13_phase_avg(1:Nx,:) + (1 / Nt) * taup_13_shifted';
   taup_11_phase_avg(1:Nx,:) = taup_11_phase_avg(1:Nx,:) + (1 / Nt) * taup_11_shifted';
        
   uW_phase_avg(1:Nx,:) = uW_phase_avg(1:Nx,:) + (1 / Nt) * uW_shifted';
   uU_phase_avg(1:Nx,:) = uU_phase_avg(1:Nx,:) + (1 / Nt) * uU_shifted';
   wW_phase_avg(1:Nx,:) = wW_phase_avg(1:Nx,:) + (1 / Nt) * wW_shifted';
   wU_phase_avg(1:Nx,:) = wU_phase_avg(1:Nx,:) + (1 / Nt) * wU_shifted';

   tau13_SGS_phase_avg(1:Nx,:) = tau13_SGS_phase_avg(1:Nx,:) + (1 / Nt) * tau13_SGS_shifted';
   tau31_SGS_phase_avg(1:Nx,:) = tau31_SGS_phase_avg(1:Nx,:) + (1 / Nt) * tau31_SGS_shifted';
   tau11_SGS_phase_avg(1:Nx,:) = tau11_SGS_phase_avg(1:Nx,:) + (1 / Nt) * tau11_SGS_shifted';
   tau33_SGS_phase_avg(1:Nx,:) = tau33_SGS_phase_avg(1:Nx,:) + (1 / Nt) * tau33_SGS_shifted';


   tau13_nu_phase_avg(1:Nx,:) = tau13_nu_phase_avg(1:Nx,:) + (1 / Nt) * tau13_nu_shifted';
   tau31_nu_phase_avg(1:Nx,:) = tau31_nu_phase_avg(1:Nx,:) + (1 / Nt) * tau31_nu_shifted';
   tau11_nu_phase_avg(1:Nx,:) = tau11_nu_phase_avg(1:Nx,:) + (1 / Nt) * tau11_nu_shifted';
   tau33_nu_phase_avg(1:Nx,:) = tau33_nu_phase_avg(1:Nx,:) + (1 / Nt) * tau33_nu_shifted';
   %disp([counter,min(u_w_phase_avg(:)),max(u_w_phase_avg(:))]);
   counter = counter + 1;


   tau13_mean_array=zeros(Nz,5);
   %tau13_tot = tau13_visc + tau13_p + tau13_SGS + tau13_turb  + tau13_wave;
   tau13_mean_array(1:Nz,1)=squeeze(mean(tau13_nu_shifted(1:Nz,1:Nx),2));%tau_visc
   tau13_mean_array(1:Nz,2)=squeeze(mean(taup_13_shifted(1:Nz,1:Nx),2));%tau_p
   tau13_mean_array(1:Nz,3)=squeeze(mean(tau13_SGS_shifted(1:Nz,1:Nx),2));%tau_SGS
   tau13_mean_array(1:Nz,4)=-squeeze(mean(uW_shifted - u_w_shifted.*W_shifted,2));%tau_turb

   for k=1:Nz
       tau13_mean_array(k,5)=-squeeze(mean((u_w_shifted(k,:) - mean(u_w_shifted(k,:),2)) ...
                                             .*(  W_shifted(k,:) - mean(W_shifted(k,:),2)),2));%tau_wave
   end

   save(sprintf('inst_data_%0.0d.mat',n), ...
         'zw','zz','x','time', ...
         'tau13_mean_array', ...
         'p_shifted_surf', ...
         'p_w_shifted_surf');

   toc
end
%tau_turb
tau13_phase_avg(:,:) = -( uW_phase_avg(1:Nx,:) - u_w_phase_avg(1:Nx,:).*W_phase_avg(:,:) );
tau31_phase_avg(:,:) = -( wU_phase_avg(1:Nx,:) -   w_phase_avg(1:Nx,:).*U_phase_avg(:,:) );
tau11_phase_avg(:,:) = -( uU_phase_avg(1:Nx,:) - u_w_phase_avg(1:Nx,:).*U_phase_avg(:,:) );
tau33_phase_avg(:,:) = -( wW_phase_avg(1:Nx,:) -   w_phase_avg(1:Nx,:).*W_phase_avg(:,:) );
%tau_wave
for k=1:Nz
    tau13_wave_phase_avg(:,k) = -(u_w_phase_avg(:,k)-mean(u_w_phase_avg(:,k),1)) ...
                               .*(W_phase_avg(:,k)-mean(W_phase_avg(:,k),1));
    tau31_wave_phase_avg(:,k) = -(w_phase_avg(:,k)-mean(w_phase_avg(:,k),1)) ...
                               .*(U_phase_avg(:,k)-mean(U_phase_avg(:,k),1));
    tau11_wave_phase_avg(:,k) = -(u_w_phase_avg(:,k)-mean(u_w_phase_avg(:,k),1)) ...
                               .*(U_phase_avg(:,k)-mean(U_phase_avg(:,k),1));
    tau33_wave_phase_avg(:,k) = -(w_phase_avg(:,k)-mean(w_phase_avg(:,k),1)) ...
                               .*(W_phase_avg(:,k)-mean(W_phase_avg(:,k),1));
    
    
end


save('phase_averaged_data_SBL_Cr1_CPS.mat','zw','z', ...
     'u_phase_avg', ...
     'u_w_phase_avg', ...
     'w_phase_avg', ...
     'p_phase_avg', ...
     'p_w_phase_avg', ...
     'eta_phase_avg', ...
     'U_phase_avg', ...
     'W_phase_avg', ...
     'taup_13_phase_avg', ...%tau_p
     'taup_11_phase_avg', ...
     'tau13_SGS_phase_avg', ...%tau_SGS
     'tau31_SGS_phase_avg',...
     'tau11_SGS_phase_avg',...
     'tau33_SGS_phase_avg',...
     'tau13_nu_phase_avg', ...%tau_visc
     'tau31_nu_phase_avg',...
     'tau11_nu_phase_avg',...
     'tau33_nu_phase_avg',...
     'tau13_phase_avg',... %tau_turb
     'tau31_phase_avg',...
     'tau11_phase_avg',...
     'tau33_phase_avg',...
     'tau13_wave_phase_avg',...%tau_wave 
     'tau31_wave_phase_avg',...
     'tau11_wave_phase_avg',...
     'tau33_wave_phase_avg');



