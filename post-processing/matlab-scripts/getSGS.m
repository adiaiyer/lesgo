function [ Cs,LM,MM,QN,NN,nusgs ] = getSGS(p,step)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for i=1:p.nproc
    
    % Open the file
    fname = ['./output/sgs.',num2str(step),'.c',num2str(i-1),'.bin'];
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
    Cs(1:p.nx,1:p.ny,zmin:zmax)=reshape(dummy,p.nx,p.ny,p.nz2);
    dummy=fread(fid,N,'double',p.fmt); 
    LM(1:p.nx,1:p.ny,zmin:zmax)=reshape(dummy,p.nx,p.ny,p.nz2);
    dummy=fread(fid,N,'double',p.fmt); 
    MM(1:p.nx,1:p.ny,zmin:zmax)=reshape(dummy,p.nx,p.ny,p.nz2);
    dummy=fread(fid,N,'double',p.fmt); 
    QN(1:p.nx,1:p.ny,zmin:zmax)=reshape(dummy,p.nx,p.ny,p.nz2);
    dummy=fread(fid,N,'double',p.fmt); 
    NN(1:p.nx,1:p.ny,zmin:zmax)=reshape(dummy,p.nx,p.ny,p.nz2);
    dummy=fread(fid,N,'double',p.fmt); 
    nusgs(1:p.nx,1:p.ny,zmin:zmax)=reshape(dummy,p.nx,p.ny,p.nz2);
    
    fclose(fid);

end
end

