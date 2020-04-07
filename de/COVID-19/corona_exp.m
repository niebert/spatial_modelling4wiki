% -----------------------------------------
% Exponential Growth Modelling for COVID-19
% Wikiversity Learning Resource: https://...
% Author: Prof.Dr. Anna Hundermark (University Koblenz-Landau)
% License: MIT oder GPL ???
% Description:  This script was developed for the students of the University Koblenz-Landau
% in course Spatial Modelling:
% https://de.wikiversity.org/wiki/Kurs:R%C3%A4umliche_Modellbildung
% -----------------------------------------

clear all
close all

%data WHO
%days (from the 1st of march)
%(not used)
 times=[3  4 5 6    9 10 11];
 inf_falle=[165 226 400 639     1156 1281 1500];
 
 
 % Data WHO or John Hopkins
timesWHO=[ 1 2  3  4 5 6  7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22  24 25 26 28 29 30 ];
inf_falleWHO=[  57 129 157 196 262 534  639  706   1112 1139 1296 1567 2369 3062  4800 6012 8000 11000 14000 19848 22000 23974  31370 37323 43646 53000 61120 66125 ];

 
 n=length(timesWHO);
 times_long=[timesWHO timesWHO(n).+timesWHO];
 
 figure(1)
 plot (times, inf_falle,'*', "linewidth",2 , timesWHO,inf_falleWHO, 's', "linewidth",2 ) 
 h=get(gcf, "currentaxes");
set(h, "fontsize", 12, "linewidth", 2); 
set(gca,'xtick',[0,5,10,11,12,13,14,15,16,17,18,20, 22, 24,26,28,30])
set(gca,'ytick',[0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,15000,20000,25000 3000 40000 50000 60000 70000])
 grid on
 %legend ("further data", "WHO data")
 title ('German corona infection increase: exponential regression')
  hold on 
 
 %initial value
  a=inf_falleWHO(1);
  %approximate growth rate,(estimated)
  k=0.26;
 expon=@(t,k) a*exp(k.*t)
 
  
 %regression for growth rate k:
 MaxN=100;
 for i=-MaxN:MaxN/5
   kn=k+i*k/MaxN;
   res(i+MaxN+1)=norm (expon(timesWHO,kn)-inf_falleWHO,2);
  
   if res(i+MaxN+1)==min(res)
     Ind=i;
     kopt= kn;
   endif
   
 endfor

 %optimal growth rate (minimal residue) ist:
 Optim_Wachsrate=kopt
 
 plot ( timesWHO, expon(timesWHO, kopt), '-',"linewidth",3 );

resi=norm (expon(timesWHO,kopt)-inf_falleWHO,2);
text (17,17000, ['Res-opt=' num2str(resi,3)], "fontsize", 18)
xlabel('days, march 2020')
  
  
  
%-----------------------linear regression fit  for growth rate -----------------------
 for i=1:n-2
wachsrate(i)=(inf_falleWHO(i+2)-inf_falleWHO(i+1))/inf_falleWHO(i+1);
endfor
k_mittelwert=sum (wachsrate)/(n-2); 
t_mittel=sum(timesWHO(2:n-1))/(n-2);

%linear function
linear=@(t,c) (t-t_mittel)*c+k_mittelwert;

 timesWHO2=timesWHO(2:n-1);


 
 c=-0.007;
 MaxN=10;
 for i=-MaxN:MaxN
   cn=c+i*1*abs(c)/MaxN;
   res1(i+MaxN+1)=norm (linear(timesWHO2,cn)-wachsrate,2);
 
   if res1(i+MaxN+1)==min(res1)
     Ind=i;
     copt= cn;
   endif
   
 endfor
 %copt
 %min(res1)
   
 
  %---------------- figures ------------------------
 figure (2)
  plot(res, '*-')
  legend ("Residuum, exponential")
   h=get(gcf, "currentaxes");
  set(h, "fontsize", 14, "linewidth", 2);
 
 
 figure (4)
  plot(res1, '*-')
  legend ("Residuum, wachstumsrate")
   h=get(gcf, "currentaxes");
  set(h, "fontsize", 14, "linewidth", 2);
  
 figure (3)
 hold on
 plot (timesWHO2, wachsrate,'*', "linewidth",3 )
 
   h=get(gcf, "currentaxes");
  set(h, "fontsize", 14, "linewidth", 2);
  
  plot (timesWHO2,linear(timesWHO2, copt), '-', "linewidth",3,t_mittel,k_mittelwert, '*', "linewidth",3  )
   set(gca,'ytick',[0,0.1, 0.15,0.17,0.2, 0.3,0.4,0.5])
   grid on
  legend ("Wachstumsrate k","Lineare Regression für k" )
  hold off;
 
 