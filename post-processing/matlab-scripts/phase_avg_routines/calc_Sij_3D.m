function [Sij]=calc_Sij_3D(kx_3D,ky_3D,dz,zetaz,u,v,w,Nx,Ny,Nz,zw)
   
   dudx=zeros(Nx,Ny,Nz);
   dvdx=zeros(Nx,Ny,Nz);
   dwdx=zeros(Nx,Ny,Nz);

   
   n=1;
   dudx(1:Nx,:,:)=ddx_i_Fourier_3D(u(1:Nx,:,:),kx_3D,Nx,n);
   dvdx(1:Nx,:,:)=ddx_i_Fourier_3D(v(1:Nx,:,:),kx_3D,Nx,n);
   dwdx(1:Nx,:,:)=ddx_i_Fourier_3D(w(1:Nx,:,:),kx_3D,Nx,n);

   dudy=zeros(Nx,Ny,Nz);
   dvdy=zeros(Nx,Ny,Nz);
   dwdy=zeros(Nx,Ny,Nz);

   n=2;
   dudy(:,1:Ny,:)=ddx_i_Fourier_3D(u(:,1:Ny,:),ky_3D,Ny,n);
   dvdy(:,1:Ny,:)=ddx_i_Fourier_3D(v(:,1:Ny,:),ky_3D,Ny,n);
   dwdy(:,1:Ny,:)=ddx_i_Fourier_3D(w(:,1:Ny,:),ky_3D,Ny,n);
   
   dudzeta=zeros(Nx,Ny,Nz);
   dvdzeta=zeros(Nx,Ny,Nz);
   dwdzeta=zeros(Nx,Ny,Nz);
   dudz=zeros(Nx,Ny,Nz);
   dvdz=zeros(Nx,Ny,Nz);
   dwdz=zeros(Nx,Ny,Nz);
   % 
   % dudzeta(:,:,2:Nz-1)=(u(:,:,3:Nz)-u(:,:,1:Nz-2))./(dz(:,:,3:Nz)+dz(:,:,2:Nz-1));
   %      dudzeta(:,:,1)=(u(:,:,2)-u(:,:,1))./dz(:,:,1);
   %     dudzeta(:,:,Nz)=(u(:,:,Nz)-u(:,:,Nz-1))./dz(:,:,Nz);
   % 
   % dvdzeta(:,:,2:Nz-1)=(v(:,:,3:Nz)-v(:,:,1:Nz-2))./(dz(:,:,3:Nz)+dz(:,:,2:Nz-1));
   %      dvdzeta(:,:,1)=(v(:,:,2)-v(:,:,1))./dz(:,:,1);
   %     dvdzeta(:,:,Nz)=(v(:,:,Nz)-v(:,:,Nz-1))./dz(:,:,Nz);
   % 
   % dwdzeta(:,:,2:Nz-1)=(w(:,:,3:Nz)-w(:,:,1:Nz-2))./(dz(:,:,3:Nz)+dz(:,:,2:Nz-1));
   %      dwdzeta(:,:,1)=(w(:,:,2)-w(:,:,1))./dz(:,:,1);
   %     dwdzeta(:,:,Nz)=(w(:,:,Nz)-w(:,:,Nz-1))./dz(:,:,Nz);


   % for i=1:Nx
   %     for j=1:Ny
   %         dudzeta(i,j,:)=gradient(squeeze(u(i,j,:)))./gradient(zw);
   %         dvdzeta(i,j,:)=gradient(squeeze(v(i,j,:)))./gradient(zw);
   %         dwdzeta(i,j,:)=gradient(squeeze(w(i,j,:)))./gradient(zw);           
   %     end
   % end

   zw_3D = repmat(reshape(zw, 1, 1, []), Nx, Ny, 1);
   dudzeta(:,:,2:end-1) = (u(:,:,3:end) - u(:,:,1:end-2)) ./ (zw_3D(:,:,3:end) - zw_3D(:,:,1:end-2));
   dudzeta(:,:,1) = (u(:,:,2) - u(:,:,1)) ./ (zw_3D(:,:,2) - zw_3D(:,:,1));

   dudz(:,:,:) = zetaz.*dudzeta;
   dvdz(:,:,:) = zetaz.*dvdzeta;
   dwdz(:,:,:) = zetaz.*dwdzeta;

   Sij=zeros(Nx,Ny,Nz,9);
   
   Sij(:,:,:,1)=0.5*(dudx+dudx);%S11
   Sij(:,:,:,2)=0.5*(dudy+dvdx);%S12
   Sij(:,:,:,3)=0.5*(dwdx+dudz);%S13


   Sij(:,:,:,4)=0.5*(dvdx+dudy);%S21
   Sij(:,:,:,5)=0.5*(dvdy+dvdy);%S22
   Sij(:,:,:,6)=0.5*(dwdy+dvdz);%S23

   Sij(:,:,:,7)=0.5*(dwdx+dudz);%S31
   Sij(:,:,:,8)=0.5*(dwdy+dvdz);%S32
   Sij(:,:,:,9)=0.5*(dwdz+dwdz);%S33
   
   
end