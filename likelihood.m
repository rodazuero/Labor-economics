function l=likelihood(x,Y,S,A,B,W,L)

%%Variable names:

%Note that gamma3 and gamma4 are re-scaled

gamma0=x(1);
gamma1=x(2);
gamma2=x(3);
gamma3=x(4);
gamma4=x(5);
alpha0=x(6);
alpha1=x(7);
stheta=x(8);
sxi=x(9);

%% We will find the three cutoff values and the 2 emax functions
%  for every period and for every individual:

XIMAT=zeros(3,1000,15);
EMAXT=zeros(15,15);


%These two functions will help us for the truncated log-normal
%expectations when double truncation or just one:


%Double truncation
TREXP = @(xiR,xiT,sxi)((normcdf(sxi-xiR/sxi)-normcdf(sxi-xiT/sxi))/(normcdf(xiT/sxi)-normcdf(xiR/sxi)));


%Single truncation
TREXP2=@(ximax,sxi)(normcdf(sxi-ximax/sxi)/(1-normcdf(ximax/sxi)));

%Non random component of wages:
WNRT =@(S,A,B) (gamma0+gamma1*S+gamma2*A-gamma3*A^2-gamma4*B);

lik=0;
%% Begin with the loop:

for i=1:1000
    for t=15:-1:1
            ii=15*(i-1)+t;
            %Actual realized experience
            EXP=A(ii);
            %Non random component of wages as function of experience
            WNRT =@(A) ((gamma0+gamma1*S(ii)+gamma2*A-gamma3*A^2-gamma4*B(ii)));
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
                    
                    %I should only store the true values
                    if A(ii)==j
                        XIMAT(1,i,t)=xire;
                        XIMAT(2,i,t)=xitax;
                        XIMAT(3,i,t)=ximax;
                    end
                    

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
                    %I should only store the true values
                    if A(ii)==j
                        XIMAT(1,i,t)=xire;
                        XIMAT(2,i,t)=xitax;
                        XIMAT(3,i,t)=ximax;
                    end
                    
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
            
            
         %After estimating all the likelihood and cutoff rules, I proceed
         %with getting the likelihood function.
         
         
         %I estimated cutoff rules for different potential levels of 
         %experience but I have to call the actual real ones:
         xire=XIMAT(1,i,t);

         
         %Likelihood contribution for individual i time t:
         
         %Term going inside of cdf
         if L(ii)==1
            incdf=(xire-sxi^2/(sxi^2+stheta^2)*(log(W(ii))-WNRT(EXP)))/(sxi*sqrt(1-(sxi^2/(sxi^2+stheta^2))));
         
         %Term going inside of pdf:
            inpdf=(log(W(ii))-WNRT(EXP))/(sqrt(sxi^2+stheta^2));
         
         
         %Adding contribution to the likelihood if person works
         lik=lik+(L(ii)==1)*log(max((1-normcdf(incdf))*normpdf(inpdf)/(sqrt(sxi^2+stheta^2)),1e-100));
         end
         %Adding contribution to the likelihood if person doesn't work:
         lik=lik+(L(ii)==0)*log(max(normcdf(xire/sxi),1e-100));
         
         %Delete this part, only to veryfy when the loop gets NaN:
         
    end
end

l=-lik;
end

