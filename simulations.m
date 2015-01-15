%% Problem Set III

%% For simulation (Store all the values of cutoffs). 


%The parameters obtained are:

xo=x0;

%Loading the data
a=load('PS3_data.csv');
AGE=a(:,2);
L=a(:,3);
W=a(:,4);
Y=a(:,5);
S=a(:,7);
A=a(:,6);
B=a(:,8);

%% Loading parameters


gamma0=xo(1);
gamma1=xo(2);
gamma2=xo(3);
gamma3=xo(4);
gamma4=xo(5);
alpha0=xo(6);
alpha1=xo(7);
stheta=xo(8);
sxi=xo(9);

%% We will find the three cutoff values and the 2 emax functions
%  for every period and for every individual:

%Also, note that I will have to generate the shock and the wages

%Setting the seed
rng(1);

%Generate the shocks: (This should be done 50 times, think about how to
%better accomplish this.

SHOCKS=normrnd(0,sxi,1000,15,50);

%In this vector I will store the history of labor force participation
%For each individual, time and for each simulation. 
LABOR=zeros(1000,15,50);



%These two functions will help us for the truncated log-normal
%expectations when double truncation or just one:

%Double truncation
TREXP = @(xiR,xiT,sxi)((normcdf(sxi-xiR/sxi)-normcdf(sxi-xiT/sxi))/(normcdf(xiT/sxi)-normcdf(xiR/sxi)));


%Single truncation
TREXP2=@(ximax,sxi)(normcdf(sxi-ximax/sxi)/(1-normcdf(ximax/sxi)));

%Wages total (random and non-random components)

XIMAT=zeros(1000,3,15,15);
%% Begin with the loop:

%The first loop will be for individuals
    for i=1:1000
        
        %Clean the EMAXT and ZIMAT matrices. Previously I included an
        %additional dimension for individuals but somehow this increased
        %computational burden by a lot. 
        
        %Storing Emax values for different levels of experience
        EMAXT=zeros(15,15);
        
        %Storing tresholds for different levels of experience
        
        %Time period for each individual
        for t=15:-1:1
            
                %Obtain the index of location for tXi
                ii=15*(i-1)+t;
                %Actual realized experience
                EXP=A(ii);
                %Non random component of wages as function of experience
                WNRT =@(A) ((gamma0+gamma1*S(ii)+gamma2*A-gamma3*A^2-gamma4*B(ii)));
                
                %Last period to be treated different as no dynamics into
                %consideration
                
                if t==15
                    %Experience at the begining of time:
                    b=ii-14; %Time indicator for beginning

                    %Initial experience:
                    AI=A(b);

                    %I loop over all possible experience levels
                    for uu=1:1:15; %This should be a 14 as experience of 15 is not important for end period
                         j=AI+uu-1; %Index of experience
                         %% Shocks
                         %  Reservation shock in last period:
                         xire=log(max(alpha1*Y(ii)+alpha0,1e-270))*(Y(ii)*alpha1+alpha0<3)+log(max(2*(alpha1*Y(ii)+alpha0-1.5),1e-270))*((Y(ii)*alpha1+alpha0>=3))-WNRT(j);


                        %  Shock tax in last period
                        xitax=log(3)-WNRT(j);

                        %  Maximum shock
                        ximax=max(xire,xitax);

                        %Previously I was storing the true values
                        %exclusively. Now I will store for every potential
                        %level of experience:
                        XIMAT(i,1,uu,t)=xire;
                        XIMAT(i,2,uu,t)=xitax;
                        XIMAT(i,3,uu,t)=ximax;



                        %% Emax functions
                        %1.  From not working:
                        UNWT=(Y(ii)*(1+alpha1)+alpha0)*normcdf(xire/sxi);

                        %2. From working and not being taxed:
                        if xire<xitax
                            UWNT=(Y(ii)+exp(WNRT(j))*exp(sxi^2/2)*TREXP(xire,xitax,sxi))*(normcdf(xitax/sxi)-normcdf(xire/sxi));
                        else
                            UWNT=0;
                        end;
                       
                        %3. From working and being taxed:
                        UWT=(Y(ii)+1.5+0.5*exp(WNRT(j)+(sxi^2)/2)*min(TREXP2(ximax,sxi),1.0e+20))*(1-normcdf(ximax/sxi));
                        EMAXT(t,uu)=UWT+UWNT+UNWT;
                        
                        
                    end;
                else  
                    %Now for periods below 15:
                    %Getting experience in initial time:
                    b=ii-(t-1); %Getting initial observation of individual i
                    AI=A(b); %Getting initial experience
                    for uu=1:1:t
                        j=AI+uu-1;

                        %% Shocks cutoff levels:

                        %1. Reservation wage shock

                        %1.1 If below taxing region
                        if alpha1*Y(ii)+alpha0+0.95*(EMAXT(t+1,uu)-EMAXT(t+1,uu+1))<3
                            xire=log(max(alpha1*Y(ii)+alpha0+0.95*(EMAXT(t+1,uu)-EMAXT(t+1,uu+1)),1e-27))-WNRT(j);

                        %1.2 If  above taxing region
                        else 
                            xire=log(max(2*(alpha1*Y(ii)+alpha0+0.95*(EMAXT(t+1,uu)-EMAXT(t+1,uu+1))-1.5),1e-27))-WNRT(j);
                        end


                        %2. Taxing wage shock
                        xitax=log(3)-WNRT(j);

                        ximax=max(xire,xitax);
                        %Previously I was storing the true values
                        %exclusively. Now I will store for every potential
                        %level of experience:
                        %
                        
                        XIMAT(i,1,uu,t)=xire;
                        XIMAT(i,2,uu,t)=xitax;
                        XIMAT(i,3,uu,t)=ximax;
                        

                        %% Emax functions

                        %1. From not working
                        UNWT=(Y(ii)*(1+alpha1)+alpha0+0.95*EMAXT(t+1,uu))*normcdf(xire/sxi);

                        %2. From working and not being taxed:
                        if xire<xitax
                            UWNT=(Y(ii)+exp(WNRT(j)+sxi^2/2)*TREXP(xire,xitax,sxi)+0.95*EMAXT(t+1,uu+1))*(normcdf(xitax/sxi)-normcdf(xire/sxi));
                        else
                            UWNT=0;
                        end;

                        %3. From working and being taxed:
                        UWT=(Y(ii)+1.5+0.5*exp(WNRT(j)+sxi^2/2)*TREXP2(ximax,sxi)+0.95*EMAXT(t+1,uu+1))*(1-normcdf(ximax/sxi));
                        EMAXT(t,uu)=UWT+UWNT+UNWT;
                        
                    end
                end 


                

        end
        
                        %Previously, this section was devoted to estimating the
                %likelihood function. Now, I am not interested in it but on
                %getting the optimal path. 
                

                
                %Simulation 
                 %For each possible simulation
                
                for rr=1:1:50
                    %For each time period
                    
                    %Initialize with initial experience
                    
                    for tt=1:1:15
                        ii=15*(i-1)+tt;
                         b=ii-(tt-1);
                         %Difference between initial and actual experience
                         uu=A(ii)-A(b)+1;
                        if SHOCKS(i,tt,rr)>XIMAT(i,1,uu,tt)
                            LABOR(i,tt,rr)=1;
                        end
                    end
                end
        
    end
%Before estimating any average for the simulations, I need to collapse
%across the dimension of the 50 replications:

L2=mean(LABOR,3);
% Rearranging:

v= [ L2(:) ];

%% Once obtained the simulations we report accordingly to 
%what we are told:

%I. Simulated number of periods working


%1. For all individuals
sum(v)/1000 
sum(L)/1000

%2. S<12
%Count how many individuals have less than 12:
sum((S<12))/15
%We get that a total of 102 individuals have less than 12

%Averaging for all individuals who have less than 12:
sum(v(S<12))/102
sum(L(S<12))/102

%3. S=12

sum((S==12))/15
%We get that a total of 281 individuals have less than 12
sum(v(S==12))/281
sum(L(S==12))/281

%4. S \in [13,15]

sum((S<=15 & S>=13))/15
%We get that a total of 269 individuals have less than 12
sum(v(S<=15 & S>=13))/269
sum(L(S<=15 & S>=13))/269

%5. S >=16

sum((S>=16))/15
%We get that a total of 348 individuals have less than 12
sum(v(S>=16))/348
sum(L( S>=16))/348


%% II. Fraction of women working at each age

SIMAGE=zeros(15,2);
for aa=1:1:15
    SIMAGE(aa,1)=mean(v(AGE==aa));
    SIMAGE(aa,2)=mean(L(AGE==aa));
end


%% III. Fraction of women working by work experience
SIMWORKE=zeros(45,2)

%10 years or less

mean(v(A<10))
mean(L(A<10))

%11-20

mean(v(A>=10 & A<=20))
mean(L(A>=10 & A<=20))

%21-more

mean(v(A>=21))
mean(L(A>=21))



%% Predicted wages by education level:


%PREDICTED WAGES FOR WORKERS:
%Periods working for each individual:

%Number of periods worked in the 50 simulations
WO=sum(LABOR,3);

%Predicted wage in each observation
PREDW=zeros(15000,1);
for i=1:1:1000
    for t=1:1:15
        %For each period-individual a vector of predicted wage
        PW=zeros(50,1);
        
        %Indicator of location for time-individual
        ii=15*(i-1)+t;
        
        %Defining non random component of wages
        WNRT = (gamma0+gamma1*S(ii)+gamma2*A(ii)-gamma3*A(ii)^2-gamma4*B(ii));
        
        for rr=1:1:50
            %Predicted wage including shocks-Taking into account only
            %workers
            PW(rr)=exp(WNRT+SHOCKS(i,t,rr))*(LABOR(i,t,rr)==1);
        end
        PREDW(ii)=sum(PW)/WO(ii);
    end
end

%Less than 12
mean(PREDW(S<12))
mean(W(S<12 & L==1))

%12
mean(PREDW(S==12))
mean(W(S==12 & L==1))


%13-15
mean(PREDW(S>=13 & S<=15))
mean(W(S>=13 & S<=15 & L==1))

%16+
mean(PREDW(S>=16 & L==1))
mean(W(S>=16 & L==1))

%Predicted wages by age:
PREDWB=zeros(15,2);

for i=1:1:15
   PREDWB(i,1)=mean(PREDW(AGE==i));
   PREDWB(i,2)=mean(W(L==1 & AGE==i));
end

L11=v;
