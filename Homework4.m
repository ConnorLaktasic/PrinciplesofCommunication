%% Homework 4 14.2, 14.13, 14.14, 14.16

%% 14.2
%a.)
N = 6;
%xi = binary representation
px1 = 1/6;
px2 = 1/6;
px3 = 1/6;
px4 = 1/6;
px5 = 1/6;
px6 = 1/6;

Ix1 = log2(1/px1)
Ix2 = log2(1/px2)
Ix3 = log2(1/px3)
Ix4 = log2(1/px4)
Ix5 = log2(1/px5)
Ix6 = log2(1/px6)

%b.)
px1 = 1/4;
px2 = 1/4;
px3 = 1/8;
px4 = 1/8;
px5 = 1/4;
px6 = 1/4;

Ix1 = log2(1/px1)
Ix2 = log2(1/px2)
Ix3 = log2(1/px3)
Ix4 = log2(1/px4)
Ix5 = log2(1/px5)
Ix6 = log2(1/px6)

%% 14.13
m = 1000;
p = 1/15;
s= 1.0;
x = pam(m,4,s);
L=sqrt(1/5);
n=sqrt(p)*randn(1,m);
y=x+n;
qy=quantalph(y,[-3*L,-L,L,3*L]);
err = sum(abs(sign(qy'-x)))/m;

%The noise for two-level transmission is the most. The noise for a four
%level transmission is the second best and the noise for the a 6-level
%transmission is the best. As the levels increase the noise decreases. They
%are inversly related.

%% 14.14
m = 1000;
p = 0.01;
s= 1.0;
x = pam(m,4,s);
L=sqrt(1/5);
n=sqrt(p)*randn(1,m);
y=x+n;
qy=quantalph(y,[-3*L,-L,L,3*L]);
err = sum(abs(sign(qy'-x)))/m

%The error probability is 0.892. The power S that is required to make the
%two-level and six-level tranmission is 2.34. I cannot think of a way to
%calculate this value.

%% 14.16
N = 5;
px1 = 1/16;
px2 = 1/8;
px3 = 1/4;
px4 = 1/16;
px5 = 1/2;

%a.) What is the entropy of this source?
H = (1/16)*log2(1/px1) + (1/8)*log2(1/px2) + (1/4)*log2(1/px3) + (1/16)*log2(1/px4) + (1/2)*log2(1/px5)
%H = 1.875

%b._ Build the Huffman chart.
%                   1        |
%x5 ------------------------|
%                   0      | 
%                1        |
%x3 ---------------------| 0.5
%                0      |
%           1          |
%x2 ------------------| 0.25
%           0        |
%      1           |
%x1 -------------| 0.125
%      0       |  
%            |
%x4 -------|
%
%

%c.) Show that the huffman code is x1<>0001, x2<>001, x3<>01, x4<>0000,
%x5<>1.
%Going down the list we can see that x5 is the first branch therefore
%x5<>1. x3 is the second branch therefore x3<>01. x2 is the third branch
%therfore x2<> 001. x1 is the fourth branch therefore x1<>0001. x4 is the
%last branch therefore x4<>0000;

%d.) What is the efficiency of this code?
E = H/(14/5)
% E = 0.6696

%e.) If this source were encoded naively, how many bits per symbol would be
%needed? What is the efficiency?
% If the code was coded naively there would be four bits per symbol. The
% efficiency of this would be .
E = H/(20/5)
%E = 0.4688




