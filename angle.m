x = segment(1).data;
figure(21), plot(x);
y = fft(x);
z = fftshift(y);
fs = 2000;
t = 0:1/fs:1-1/fs;
ly = length(y);
f = (-ly/2:ly/2-1)/ly*fs;

stem(f,abs(z))
xlabel 'Frequency (Hz)'
ylabel '|y|'
grid

tol = 1e-6;
z(abs(z) < tol) = 0;

theta = angle(z);

stem(f,theta/pi)
xlabel 'Frequency (Hz)'
ylabel 'Phase / \pi'
grid