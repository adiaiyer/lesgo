function [dui_dxi] = ddx_i_Fourier_3D(ui, ki, N, n)
    % Calculate Fourier transform of ui
    uhati = fft(ui, [], n) / N;
    
    % Calculate derivative in Fourier space
    dui_dxi = real(ifft(sqrt(-1) * ki .* uhati, [], n)) * N;
end