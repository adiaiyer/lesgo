function f_ph_avg = phase_avg_shift(f,kx,Nx,time,c) 
   f_hat = fft(f)/Nx;
   f_ph_avg = real(ifft(f_hat.*exp(1j*kx*c*time))*Nx);
end
