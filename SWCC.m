clc 
clear
nsp=51; %Every peak and trough of the theoretical and actual magnetic anomalies are interpolated into nsp equal data points
A=load('C9Window.dat'); % Theoretical window
B=load('Magnetic_Profile.dat'); % Actual marine magnetic anomaly profile
dt=A(:,2);
da=B(:,2);
xt=A(:,1);
xa=B(:,1);
nt=length(dt);
na=length(da);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Search for positions of zero magnitude of the theoretical magnetic anomalies
kt=0;
skt=0; 
for i=1:nt-1    
  if (dt(i)==0)
      kt=kt+1;
  elseif(dt(i+1)/dt(i)<0)
      kt=kt+1; 
      skt=skt+1;
  end
end
ntz=kt;
dt0x=zeros(ntz,1);
dt0d=zeros(ntz,1);
kt=0;
for i=1:nt-1
  if (dt(i)==0)   
      kt=kt+1;
      dt0x(kt)=xt(i);
      dt0d(kt)=0.0;
  elseif(dt(i+1)/dt(i)<0)
      kt=kt+1;
      slope=(dt(i+1)-dt(i))/(xt(i+1)-xt(i));
      dt0x(kt)=-1.0*(dt(i)/slope)+xt(i);
      dt0d(kt)=0.0;
  end
end
xts=zeros(nt+skt,1);
dts=zeros(nt+skt,1); 
 kt=0;
 for i=1:nt-1
    dts(i+kt)=dt(i);
    xts(i+kt)=xt(i);
   if (dt(i)==0) 
       kt=kt+1;
    dts(i+kt)=dt(i);
    xts(i+kt)=xt(i);      
  elseif(dt(i+1)/dt(i)<0)
      kt=kt+1;
    dts(i+kt)=dt0d(kt);
    xts(i+kt)=dt0x(kt);
   end
 end
 nts=nt+kt;
 xts(nts)=xt(nt);
 dts(nts)=dt(nt); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Search for positions of zero magnitude of the actual magnetic anomalies
kt2=0;
skt2=0;
for i=1:na-1    
  if (da(i)==0)then
      kt2=kt2+1;
  elseif(da(i+1)/da(i)<0)
      kt2=kt2+1;
      skt2=skt2+1;
  end
end
naz=kt2;
da0x=zeros(naz,1);
da0d=zeros(naz,1);
kt2=0;
for i=1:na-1
  if (da(i)==0)   
      kt2=kt2+1;
      da0x(kt2)=xa(i);
      da0d(kt2)=0.0;
  elseif(da(i+1)/da(i)<0)
      kt2=kt2+1;
      slope=(da(i+1)-da(i))/(xa(i+1)-xa(i));
      da0x(kt2)=-1.0*(da(i)/slope)+xa(i);
      da0d(kt2)=0.0;
  end
end
xas=zeros(na+skt2,1);
das=zeros(na+skt2,1);
 kt2=0;
 for i=1:na-1
    das(i+kt2)=da(i); 
    xas(i+kt2)=xa(i);
   if (da(i)==0)
       kt2=kt2+1;
    xas(i+kt2)=xa(i); 
    das(i+kt2)=da(i);         
  elseif(da(i+1)/da(i)<0)
      kt2=kt2+1;
    das(i+kt2)=da0d(kt2);
    xas(i+kt2)=da0x(kt2);
   end
 end
 nas=na+skt2;
 xas(nas)=xa(na); 
 das(nas)=da(na);
 data=[xas,das]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%interpolate every peak and trough of the theoretical magnetic anomalies
%into nsp data points
magtx=zeros((nsp*(ntz-1)-(ntz-2)),1);
magtd=zeros((nsp*(ntz-1)-(ntz-2)),1);
ntp=nsp*(ntz-1)-(ntz-2);
for i=1:ntz-1 
    ks=0;    
    for j=1:nts
      if((xts(j)>=dt0x(i))&&(xts(j)<=dt0x(i+1)))  
         ks=ks+1;
         sx(ks)=xts(j);
         sy(ks)=dts(j);       
      end     
    end
    xmin=min(sx);
    xmax=max(sx);  
    interx=(xmax-xmin)/(nsp-1);
    cx=xmin:interx:xmax;
    cy=interp1(sx,sy,cx,'spline');
    magtx(1+(i-1)*(nsp-1):(i*nsp-i+1))=cx;
    magtd(1+(i-1)*(nsp-1):(i*nsp-i+1))=cy;    
     sx=[];
     sy=[];
end
 figure (1) 
 plot(magtx,magtd)
 axis([min(magtx),max(magtx),min(magtd),max(magtd)]) 
 title('Theoretical window')
 xlabel('Distance')
 ylabel('Magnetic anomaly')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%interpolate every peak and trough of the actual magnetic anomalies
%into nsp data points
magax=zeros((nsp*(naz-1)-(naz-2)),1);
magad=zeros((nsp*(naz-1)-(naz-2)),1); 
for i=1:naz-1 
    ks2=0;    
    for j=1:nas
      if((xas(j)>=da0x(i))&&(xas(j)<=da0x(i+1)))   
         ks2=ks2+1;
         sx(ks2)=xas(j);
         sy(ks2)=das(j);         
      end     
    end
    xmin=min(sx);
    xmax=max(sx);
    interx2=(xmax-xmin)/(nsp-1);
    cx=xmin:interx2:xmax;
    cy=interp1(sx,sy,cx,'spline');
    magax(1+(i-1)*(nsp-1):(i*nsp-i+1))=cx; 
    magad(1+(i-1)*(nsp-1):(i*nsp-i+1))=cy;    
     sx=[]; 
     sy=[];
end
figure(2)
plot(magax,magad)
axis([min(magax),max(magax),min(magad),max(magad)])
title('Actual magnetic anomaly profile')
xlabel('Distance')
ylabel('Magnetic anomaly')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Identify marine magnetic anomalies with PCC, SRCC and KRCC by the sliding
%window technique
choosex=zeros(ntp,1);
choosed=zeros(ntp,1);
rp=zeros(naz+1-ntz,1);
rk=zeros(naz+1-ntz,1);
rs=zeros(naz+1-ntz,1);

for i=1:(naz+1-ntz)
    for k=1:ntp
      choosex(k)=magax(k+(i-1)*nsp-(i-1));
      choosed(k)=magad(k+(i-1)*nsp-(i-1));
    end    
   rp(i)=corr(choosed,magtd, 'type' ,'Pearson'); 
   rs(i)=corr(choosed,magtd, 'type' ,'Spearman'); 
   rk(i)=corr(choosed,magtd, 'type' ,'Kendall');  
   minx=min(choosex);
   maxx=max(choosex);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output identification results: the first column is the sliding step, the second
%column is the identification result of the PCC, the third column is the identification
%result of the SRCC and the third column is the identification result of the KRCC

n=length(rp);
step=1:n;
data=[step',rp,rs,rk];
save('SWCC_Results.dat','data','-ascii')
[m1,index1]=max(rp);
[m2,index2]=max(rs);
[m3,index3]=max(rk);
if(index1==index2&&index2==index3)
    Sliding_step=index1;
   'The position of the identified corresponding magnetic anomalies is in the', Sliding_step
end


    

