% authors: nadav goldin, amichai horvitz
% made for SSP course, huji 2015
% Example of simple sin recovery sum using OMP with Gaussian matrix
% versus traditonal nyquist rate sampling

% parameters
amplitude = 100;
freqHz1 = 100;
freqHz2 = 900;
maxFreq = max (freqHz1, freqHz2);
OMP_samples = 185;
fsHz = freqHz2*100; 
fsHz_nyq=1801;


% plot original DFT way above nyquist
dt = 1/fsHz;
sine = amplitude*sin(2*pi*freqHz1*(0:dt:1-dt)) +  amplitude*sin(2*pi*freqHz2*(0:dt:1-dt));
N = length(sine);
transform = fft(sine,N)/N;
magTransform = abs(transform);
faxis = linspace(-fsHz/2,fsHz/2,N);
subplot(2,2,2)
plot(faxis,fftshift(magTransform));
xlabel('Frequency (Hz)')
axis([ floor(-maxFreq*1.2) floor(maxFreq*1.2) 0 amplitude*1.2 ])
title('Original signal DFT(way above nyquist)')
% plot time domain
subplot(2,2,1)
plot( (0:dt:1-dt), sine)
grid on;
xlabel('time[S]');
ylabel('Amplitude');
header=sprintf('Time Domain of sum of 2 sine wave %dHz, %dHz', freqHz1, freqHz2);
axis([ 0 dt*maxFreq -max(sine)*1.2 max(sine)*1.2 ])
title(header);

% plot DFT at nyquist rate
dt_below= 1/fsHz_nyq;
sine_below = amplitude*sin(2*pi*freqHz1*(0:dt_below:1-dt_below)) +  amplitude*sin(2*pi*freqHz2*(0:dt_below:1-dt_below));
N_below = length(sine_below);
transform_below = fft(sine_below,N_below)/N_below;
magTransform_below = abs(transform_below);
faxis_below = linspace(-fsHz_nyq/2,fsHz_nyq/2,N_below);
subplot(2,2,3)
plot(faxis_below,fftshift(magTransform_below));
axis([ floor(-maxFreq*1.2) floor(maxFreq*1.2) 0 amplitude*1.2 ])
xlabel('Frequency (Hz)')
header = sprintf('DFT at nyquist rate(%dHz), requires at least %d samples', fsHz_nyq, max(freqHz1, freqHz2)*2 );
title(header);


% Create Gaussian Matrix and run OMP
sparsity_level=length ( find ( magTransform > 1e-5 ) );
matrix = randn( OMP_samples, length(magTransform) );
b = matrix*magTransform';
[ x_hat , ~, ~, ~ ] = OMP( matrix, b, sparsity_level );
x_hat = x_hat';

% Plot OMP result
subplot(2,2,4)
plot(faxis,fftshift(x_hat));
xlabel('Frequency (Hz)')
error = norm ( magTransform - x_hat ) / norm (magTransform);
header = sprintf('Recovery using OMP with %d samples and random gaussian matrix, error: %f', length(b), error );
axis([ floor(-maxFreq*1.2) floor(maxFreq*1.2) 0 amplitude*1.2 ])
title(header);
