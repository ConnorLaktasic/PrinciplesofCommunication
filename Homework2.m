%% Homework 2. 9.1-9.5
%% 9.1
str = '01234 I wish I were an Oscar Meyer wiener 56789'
m = letters2pam(str);
N=length(m);
M=100;
mup = zeros(1,N*M);
mup(1:M:N*M)=m;
p=hamming(M);
x=filter(p,1,mup);
figure(1), plotspec(x,1/M)
t=1/M:1/M:length(x)/M;
fc = 50;
c = cos(2*pi*fc*t);
r=c.*x;
c2 = cos(2*pi*fc*t);
x2 = r.*c2;
f1 = 50;
fbe = [0 0.1 0.2 1];
damps = [1 1 0 0];
b = firpm(f1, fbe, damps);
x3 = 2*filter(b,1,x2);
y = filter(fliplr(p)/(pow2(p)*M),1,x3);
z = y(0.5*f1+M:M:N*M);
figure(2), plot([1:length(z)],z,'.');
mprime = quantalph(z,[-3,-1,1,3]);
cvar = (mprime-z)*(mprime-z)/length(mprime);
lmp = length(mprime);
pererr = 100*sum(abs(sign(mprime-m(1:lmp))))/lmp
reconstructed_message = pam2letters(mprime)
% Changing fc changes the carrier frequency. The frequency which the
% message is sent at changes as fc changes. When fc increases the message
% is carried at a higher frequency and when fc deceases the message is
% carried at a lower frequency. The limiting factor is how close the
% signals are together.

%% 9.2
str = '01234 I wish I were an Oscar Meyer wiener 56789'
m = letters2pam(str);
N=length(m);
M=10;
mup = zeros(1,N*M);
mup(1:M:N*M)=m;
p=hamming(M);
x=filter(p,1,mup);
figure(1), plotspec(x,1/M)
t=1/M:1/M:length(x)/M;
fc = 50;
c = cos(2*pi*fc*t);
r=c.*x;
c2 = cos(2*pi*fc*t);
x2 = r.*c2;
f1 = 50;
fbe = [0 0.1 0.2 1];
damps = [1 1 0 0];
b = firpm(f1, fbe, damps);
x3 = 2*filter(b,1,x2);
y = filter(fliplr(p)/(pow2(p)*M),1,x3);
z = y(0.5*f1+M:M:N*M);
figure(2), plot([1:length(z)],z,'.');
mprime = quantalph(z,[-3,-1,1,3]);
cvar = (mprime-z)*(mprime-z)/length(mprime);
lmp = length(mprime);
pererr = 100*sum(abs(sign(mprime-m(1:lmp))))/lmp
reconstructed_message = pam2letters(mprime)

%The nyquist rate is the limiting factor. The sampling frequencies must be
%sampled at a rate of at least 2 times the frequency of the signal to
%reconstruct the signal exactly.
%% 9.3
str = '01234 I wish I were an Oscar Meyer wiener 56789';
fc = 30;
m = letters2pam(str);
N=length(m);
M=100;
mup = zeros(1,N*M);
mup(1:M:N*M)=m;
p=hamming(M);
c2 = cos(2*pi*fc*t);
x2 = r.*c2;
f1 = 50;
fbe = [0 0.1 0.2 1];
damps = [1 1 0 0];
b = firpm(f1, fbe, damps);
x3 = 2*filter(b,1,x2);
y = filter(fliplr(p)/(pow2(p)*M),1,x3);
z = y(0.5*f1+M:M:N*M);
figure(2), plot([1:length(z)],z,'.');
mprime = quantalph(z,[-3,-1,1,3]);
cvar = (mprime-z)*(mprime-z)/length(mprime);
lmp = length(mprime);
pererr = 100*sum(abs(sign(mprime-m(1:lmp))))/lmp
reconstructed_message = pam2letters(mprime)

%If the frequency at the beginning is removed there are frequencies and
%noise outside of the area of interest that will corrupt the signal. If
%there are other user present there transmission will not be removed and
%will corrupt the signal when it is trying to be received.
%% 9.4
str = '01234 I wish I were an Oscar Meyer wiener 56789'
m = letters2pam(str);
N=length(m);
M=10;
mup = zeros(1,N*M);
mup(1:M:N*M)=m;
p=hamming(M);
x=filter(p,1,mup);
figure(1), plotspec(x,1/M)
t=1/M:1/M:length(x)/M;
fc = 50;
c = cos(2*pi*fc*t);
r=c.*x;
c2 = cos(2*pi*fc*t);
x2 = r.*c2;
f1 = 50;
fbe = [0 0.1 0.2 1];
damps = [1 1 0 0];
b = firpm(f1, fbe, damps);
x3 = 2*filter(b,1,x2);
y = filter(fliplr(p)/(pow2(p)*M),1,x3);
z = y(0.5*f1+M:M:N*M);
figure(2), plot([1:length(z)],z,'.');
mprime = quantalph(z,[-3,-1,1,3]);
cvar = (mprime-z)*(mprime-z)/length(mprime);
lmp = length(mprime);
pererr = 100*sum(abs(sign(mprime-m(1:lmp))))/lmp
reconstructed_message = pam2letters(mprime)

%The limits of the LFP design at the beginning is that it doesnt completly
%get rid of other frequencies near by, there is some overlap that makes it
%through the filter. The lowest cut off frequency is 0 and the highest
%cut off frequency is 1.

%% 9.5
str = '01234 I wish I were an Oscar Meyer wiener 56789'
m = letters2pam(str);
N=length(m);
M=10;
mup = zeros(1,N*M);
mup(1:M:N*M)=m;
p=hamming(M);
x=filter(p,1,mup);
figure(1), plotspec(x,1/M)
t=1/M:1/M:length(x)/M;
fc = 50;
c = cos(2*pi*fc*t);
r=c.*x;
c2 = cos(2*pi*fc*t);
x2 = r.*c2;
f1 = 50;
fbe = [0 0.1 0.2 1];
damps = [1 1 0 0];
b = firpm(f1, fbe, damps);
x3 = 2*filter(b,1,x2);
y = filter(fliplr(p)/(pow2(p)*M),1,x3);
z = y(0.5*f1+M:M:N*M);
figure(2), plot([1:length(z)],z,'.');
mprime = quantalph(z,[-3,-1,1,3]);
cvar = (mprime-z)*(mprime-z)/length(mprime);
lmp = length(mprime);
pererr = 100*sum(abs(sign(mprime-m(1:lmp))))/lmp
reconstructed_message = pam2letters(mprime)

%Using the same specifications as before you can make the LFP last 0.10 or
%ten percent bandwidth.
