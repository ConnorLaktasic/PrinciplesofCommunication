%%Homework 1. 7.5-7.13;7.15
%% 7.5
f = 100; Ts = 1/1000; time = 5.0;
t = Ts:Ts:time;
w=sin(2*pi*f*t);
N=2^8;
ssf = (-N/2:N/2-1)/(Ts*N);
fw = fft(w(1:N));
fws = fftshift(fw);
plot(ssf,abs(fws))

%a.) As f becomes to large it begins to alias to a lower frequency.
% f and Ts are inversely related
%b.) There is an inverse relationship between f and Ts.
%c.) As N becomes to smaller the width of the peak becomes larger and
% and as N becomes to large the width of the peak becimes narrower.
%% 7.6
f = 100; Ts = 1/1000; time = 5.0;
t = Ts:Ts:time;
w=sin(2*pi*f*t).^4;
N=2^10;
ssf = (-N/2:N/2-1)/(Ts*N);
fw = fft(w(1:N));
fws = fftshift(fw);
plot(ssf,abs(fws))

%The spectrum of sin^2 has 3 peaks, one at f = 0 with a magnitude larger
%than sin and two at f =  plus and minus 200 with a magnitude of half the
%the peak at f = 0. There is also an additional peak in at f = 0. 
%For sin^3 there are 4 peaks. The spectrum of sin^k would 5 peaks if it was
%even and have 4 peaks if it were odd. The largest k that would make sense
%is k = 4.

%% 7.7
f = 100; Ts = 1/1000; time = 5.0;
t = Ts:Ts:time;
w=sinc(2*pi*f*t).^2;
N=2^10;
ssf = (-N/2:N/2-1)/(Ts*N);
fw = fft(w(1:N));
fws = fftshift(fw);
plot(ssf,abs(fws))

%The spectrum of the sinc funtion has two wide peaks. The spectrum of
%sinc.^2 is one large peak with two local max peaks, one on each side of
%the main peak.

%% 7.8
f = 100; Ts = 1/1000; time = 5.0;
t = Ts:Ts:time;
w=sinc(2*pi*f*t)+j*exp(1).^-t;
N=2^10;
ssf = (-N/2:N/2-1)/(Ts*N);
fw = fft(w(1:N));
fws = fftshift(fw);
plot(ssf,abs(fws))

%Use specsin2 method.
%% 7.9
f = 100; Ts = 1/1000; time = 5.0;
t = Ts:Ts:time;
phi = 0;
w=sin(2*pi*f*t+phi).^2;
N=2^10;
ssf = linspace(0,180,1024);
fw = fft(w(1:N));
fws = fftshift(fw);
plot(ssf,angle(fws)*(180/pi))

%a.)The phase of stays pretty constant at low values of phi, then becomes very
%different for larger values such as 1.5 and 3.14. At 3.14 it inverts. 
%b.) The phase acts similar to before but this time there are two phase
%jumps instead of one. 

%% 7.10
filename = 'gong.wav';
[x,sr] = wavread(filename);
Ts = 1/Sr;
N=2^15;
x=x(1:N);
sound(x, 1/Ts);
time = Ts*(0:length(x)-1);
subplot(2,1,1); plot(time,x)
magx = abs(fft(x));
ssf = (0:N/2-1)/(Ts*N);
subplot(2,1,2), plot(ssf, magx(1:N/2))

% The gong sound have more lower frequency components during the first 0.1
% seconds than it does in a 0.1 second section in the middle.

%% 7.11
filename = 'gong.wav';
[x,sr] = wavread(filename);
Ts = 1/Sr;
N=2^15;
x=x(1:N);
sound(x, 1/Ts);
time = Ts*(0:length(x)-1);
subplot(2,1,1); semilogy(time,x)
magx = abs(fft(x));
ssf = (0:N/2-1)/(Ts*N);
subplot(2,1,2), semilogy(ssf, magx(1:N/2))

%You can better see the distribution of the points use the semilogy
%function than you can using the plot function. 

%% 7.12
filename = 'gong2.wav';
[x,sr] = wavread(filename);
Ts = 1/Sr;
N=2^15;
x=x(1:N);
sound(x, 1/Ts);
time = Ts*(0:length(x)-1);
subplot(2,1,1); semilogy(time,x)
magx = abs(fft(x));
ssf = (0:N/2-1)/(Ts*N);
subplot(2,1,2); semilogy(ssf, magx(1:N/2))

% The magnitude of the sound is greater for all frequencies and times
% throught the sound.

%% 7.13
filename = 'rumble.wav';
[x,sr] = wavread(filename);
Ts = 1/Sr;
N=2^15;
x=x(1:N);
sound(x, 1/Ts);
time = Ts*(0:length(x)-1);
subplot(2,1,1); semilogy(time,x)
magx = abs(fft(x));
ssf = (0:N/2-1)/(Ts*N);
subplot(2,1,2); semilogy(ssf, magx(1:N/2))

%The sounds has all of the same frequency components throught out the sound
%but the magnitude of the components changes. They start small and grow
%throughout the sound.

%% 7.15
a = [0.9];
b = 2;
lena = length(a) - 1;
lenb = length(b);
d = randn(1,20);
if lena >= lenb,
    h = impz(b,a);
    yfilt = filter(h,1,d)
end
yfilt2 = filter(b,a,d)
y = zeros(lena,1); 
x = zeros(lenb,1);
for k = 1:length(d)-lenb
    x = [d(k);x(1:lenb-1)];
    ytim(k) = -a(2:lena+1)*y+b*x;
    y=[ytim(k);y(1:lena-1)];
end
