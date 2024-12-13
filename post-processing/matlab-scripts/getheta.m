function [ theta, pi_x, pi_y, pi_z ] = getheta(p,step)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for i=1:p.nproc
    
    % Open the file
    fname = ['./output/theta.',num2str(step),'.c',num2str(i-1),'.bin'];
    fid=fopen(fname,'r');
    if (fid < 0) 
        error('getSnap:fname',['Could not open file ',fname]);
    end

    % Determine the interval of the matrix where the data should be stored
    zmin=p.zmin_buf(i);
    zmax=p.zmax_buf(i);

    % Scan the data
    N = p.nx*p.ny*p.nz2;
    dummy=fread(fid,N,'double',p.fmt); 
    theta(1:p.nx,1:p.ny,zmin:zmax)=reshape(dummy,p.nx,p.ny,p.nz2);
    dummy=fread(fid,N,'double',p.fmt); 
    pi_x(1:p.nx,1:p.ny,zmin:zmax)=reshape(dummy,p.nx,p.ny,p.nz2);
    dummy=fread(fid,N,'double',p.fmt); 
    pi_y(1:p.nx,1:p.ny,zmin:zmax)=reshape(dummy,p.nx,p.ny,p.nz2);
    dummy=fread(fid,N,'double',p.fmt); 
    pi_z(1:p.nx,1:p.ny,zmin:zmax)=reshape(dummy,p.nx,p.ny,p.nz2);
    
    
    fclose(fid);

end
end

