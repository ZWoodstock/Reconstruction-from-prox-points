function out = bandlim(x,band)
%bandlimits the real signal x by keeping the DFT indices where
%band==1.
%band: must be a binary 0-1 vector, where '1' corresponds to the
%kept Fourier indices. Note that for real signals, 'band' must be
%constructed properly (e.g. it must satisfy certain symmetry
%constraints).
out = real(ifft(band.*fft( x)));
end
