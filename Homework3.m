%% Homework 3. 11.3,11.4,12.1,12.2,12.3,13.1,13.2
%% 11.3
N = 1000;
m = pam(N,2,1);
M = 20;
mup = zeros(1, N*M);
mup(1:M:N*M) = m;
ph = hamming(M);
xh = filter(ph, 1, mup);
pr = hamming(M);
xr = filter(pr, 1, mup);
ps = hamming(M);
xs = filter(ps, 1, mup);
neye = 5;
ch = floor(length(xh)/(neye*M));
cr = floor(length(xr)/(neye*M));
cs = floor(length(xs)/(neye*M));
xp = xh(N*M-eye*M*c+1:N*M);
xo = xr(N*M-eye*M*c+1:N*M);
xi = xs(N*M-eye*M*c+1:N*M);
plot(reshape(xp,neye*M,c))
plot(reshape(xo,neye*M,c))
plot(reshape(xi,neye*M,c))

%The plots show the eye diagrams fro rectangular, hamming, and sinc pulse
%shapes.
%% 11.4
N = 1000;
m = pam(N,2,1);
M = 20;
v = 3;
mup = zeros(1, N*M);
mup(1:M:N*M) = m;
ph = hamming(M);
xh = filter(ph, 1, mup);
xh = v*randn + xh;
pr = hamming(M);
xr = filter(pr, 1, mup);
xr = v*randn + xr;
ps = hamming(M);
xs = filter(ps, 1, mup);
xs = v*randn + xs;
neye = 5;
ch = floor(length(xh)/(neye*M));
cr = floor(length(xr)/(neye*M));
cs = floor(length(xs)/(neye*M));
xp = xh(N*M-eye*M*c+1:N*M);
xo = xr(N*M-eye*M*c+1:N*M);
xi = xs(N*M-eye*M*c+1:N*M);
plot(reshape(xp,neye*M,c))
plot(reshape(xo,neye*M,c))
plot(reshape(xi,neye*M,c))

%For the rectangular diagram v can be 4 and still have an open eye. For the
%hamming diagram the v can be 7 and still have an open eye. For the sinc
%diagram the v can be  4 and still have an open eye.

%% 12.1
n = 10000;
m = 2;
beta = 0.3;
l = 50;
chan = [1];
toffset = -0.3;
pulshap = srrc(1,beta,m,toffset);
s = pam(n,4,5);
sup = zeros(1,n*m);
sup(1:m:n*m) = s;
hh = conv(pulshap,chan);
r = conv(hh,sup);
matchfilt = srrc(1,beta,m,0);
x = conv(r, matchfilt);
tnow = 1*m+1;
tau = 0;
xs = zeros(1,n);
tausave = zeros(1,n);
tausave(1) = tau;
i = 0;
mu = 0.01;
delta = 0.1;
while tnow < lenght(x)-2*1*m
    i=i+1;
    xs(i) = interpsinc(x,tnow+tau,1);
    x_deltap=interpsinc(x,tnow+tau+delta, 1);
    x_deltam=interpsinc(x,tnow+tau+delta, 1);
    dx = x_deltap-x_deltam;
    qx = quantalph(xs(i),[-3,-1,1,3]);
    tau = tau+mu*dx*(qx-xs(i));
    tnow = tnow+m;
    tausave(i) = tau;
end

%a.) As mu becomes larger the convergence rate becomes slower, meaning that it
%converges over a longer area. When mu becomes smaller the convergence rate
%becomes faster meaning that it converges over a smaller ares. a step size
%of 0.05 work.
%b.) The signal constellation of the input affects the convergent value of
%tau as a relationship of the size to the convergence rate. The more
%quantized states in the constelation, the faster the convergence rate.

%% 12.2
n = 10000;
m = 2;
beta = 0.3;
l = 50;
chan = [1];
toffset = -0.3;
pulshap = srrc(1,beta,m,toffset);
s = pam(n,4,5);
sup = zeros(1,n*m);
sup(1:m:n*m) = s;
hh = conv(pulshap,chan);
r = conv(hh,sup);
matchfilt = srrc(1,beta,m,0);
x = conv(r, matchfilt);
tnow = 1*m+1;
tau = 0;
xs = zeros(1,n);
tausave = zeros(1,n);
tausave(1) = tau;
i = 0;
mu = 0.01;
delta = 0.1;
while tnow < lenght(x)-2*1*m
    i=i+1;
    xs(i) = interpsinc(x,tnow+tau,1);
    x_deltap=interpsinc(x,tnow+tau+delta, 1);
    x_deltam=interpsinc(x,tnow+tau+delta, 1);
    dx = x_deltap-x_deltam;
    qx = quantalph(xs(i),[-3,-1,1,3]);
    tau = tau+mu*dx*(qx-xs(i));
    tnow = tnow+m;
    tausave(i) = tau;
end

%Implementing a rectangular shape is worse than the raised cosine shape of
%SRRC. This is because the SRRC does not spread out as much in the
%frequency spectrum as the rectangular shape does.

%% 12.3
n = 10000;
m = 2;
beta = 0.3;
l = 50;
chan = [1];
toffset = -0.3;
pulshap = srrc(1,beta,m,toffset);
pulshap = pulshap + randn;
s = pam(n,4,5);
sup = zeros(1,n*m);
sup(1:m:n*m) = s;
hh = conv(pulshap,chan);
r = conv(hh,sup);
matchfilt = srrc(1,beta,m,0);
x = conv(r, matchfilt);
tnow = 1*m+1;
tau = 0;
xs = zeros(1,n);
tausave = zeros(1,n);
tausave(1) = tau;
i = 0;
mu = 0.01;
delta = 0.1;
while tnow < lenght(x)-2*1*m
    i=i+1;
    xs(i) = interpsinc(x,tnow+tau,1);
    x_deltap=interpsinc(x,tnow+tau+delta, 1);
    x_deltam=interpsinc(x,tnow+tau+delta, 1);
    dx = x_deltap-x_deltam;
    qx = quantalph(xs(i),[-3,-1,1,3]);
    tau = tau+mu*dx*(qx-xs(i));
    tnow = tnow+m;
    tausave(i) = tau;
end

%Adding noise slows down the convergence of the timing-offset parameter
%tau, meaning that tau becomes larger when noise is added. It does not
%change the final convergence value.

%% 13.1
b = [0.5 1 -0.6];
m = 1000;
x = linspace(0,1000,512);
s = sign(randn(1,m));
r = filter(b,1,s);
n = 3;
for i = 0:3
    delta = i;
    p = length(r)-delta;
    R = toeplitz(r(n+1:p)),r(n+1:-1:1);
    S=s(n+1-delta:p-delta)';
    f = inv(R'*R)*R'*S;
    Jmin = S'*S-S'*R*inv(R'*R)*R'*S;
    y = filter(f,1,r);
    g = freqz(y);
    figure(i+1);
    plot(x,freqz(y))
    dec = sign(y);
    err = 0.5*sum(abs(dec(delta+1:m)-s(1:m-delta)))
end

%The products are very close to unity. When the four plots are placed side
%by side you can visually see how similar the frequency response it. The
%delay doesn't affect the frequency response significantly.

%% 13.2
b = [0.5 1 -0.6];
m = 1000;
x = linspace(0,1000,512);
s = sign(randn(1,m));
r = filter(b,1,s);
n = 3;
delta = 1;
sd = 4;
p = length(r)-delta;
R = toeplitz(r(n+1:p)),r(n+1:-1:1);
S=s(n+1-delta:p-delta)';
f = inv(R'*R)*R'*S;
Jmin = S'*S-S'*R*inv(R'*R)*R'*S;
r = filter(b,1,s)+sd*randn(size(s));
y = filter(f,1,r);
y = y + r;
g = freqz(y);
figure(i+1);
plot(x,freqz(y))
dec = sign(y);
err = 0.5*sum(abs(dec(delta+1:m)-s(1:m-delta)))

%a.) sd can be as large as 4 and not get any errors for the equalizer with
%delay of 2.
%c.)sd can be as larger as 7 and not get any erros for the equalizer with a
% delay of 1.
%d.) The equalizer with a delay of 1 is a better equalizer becausee sd can
% can be larger and not get errors than the equalizer of delay 2.