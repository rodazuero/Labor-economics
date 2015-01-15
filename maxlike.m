
%% Emax Estimation

%Loading the data
a=load('PS3_data.csv');
L=a(:,3);
W=a(:,4);
Y=a(:,5);
S=a(:,7);
A=a(:,6);
B=a(:,8);



%Given a set of parameters of:


gamma0=0;
gamma1=0.1;
gamma2=0.01;
gamma3=0;
gamma4=0;
alpha0=1;
alpha1=0.1;
stheta=0.2;
sxi=0.2;


options = optimset('Display','iter','LargeScale','off','MaxFunEvals',5000);



x=[gamma0, gamma1, gamma2, gamma3, gamma4, alpha0, alpha1, stheta, sxi];
xo=[1.0,    0.01,    0,     0,    0,     2.5,     0.03,    0.2  , 0.3];

f = @(x) likelihood(x,Y,S,A,B,W,L);
x0 = fminunc(f,xo,options)
tic
f = @(x) likelihood(x,Y,S,A,B,W,L);
f(x)
toc

