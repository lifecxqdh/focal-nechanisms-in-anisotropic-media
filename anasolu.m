function [U,Up,Us]=anasolu(isotype,fre,loca,vel,dtt,lasttt,MM,split)
if(isotype==1)
    M=[1,0,0;0,1,0;0,0,1];
end
if(isotype==2)
    M=[0,0,1;0,0,0;1,0,0];
end
if(isotype==3)
    M=[-1,0,0;0,-1,0;0,0,2];
end
if(isotype==0)
    M=MM;
end
x=loca(1);
y=loca(2);
z=loca(3);
rho=vel(3);%密度
beta=vel(2);%横波速度
alpha=vel(1);%纵波速度
r = sqrt(x^2 + y^2 + z^2);
tt=0:dtt:lasttt;
tt=tt-1.0/fre;
t1=tt-r/alpha;
t2=tt-r/beta;
M11=zeros(length(tt),1);
M22=zeros(length(tt),1);
M1=zeros(length(tt),1);
M2=zeros(length(tt),1);
MM=zeros(length(tt),1);
gamma=[x/r,y/r,z/r];
delta=[1,0,0;0,1,0;0,0,1];
u=zeros(3,length(tt));
up=zeros(3,length(tt));
us=zeros(3,length(tt));
F=@(x,t)(x.*(1-2*(pi*fre*(t-x)).^2).*exp(-(pi*fre*(t-x)).^2));
% F=@(x,t)(x.*((t-x).*exp(-(pi*fre*(t-x)).^2)));
AN=1/(4*pi*rho*r^4);
AIP=1/(4*pi*rho*alpha^2*r^2);
AIS=1/(4*pi*rho*beta^2*r^2);
AFP=1/(4*pi*rho*alpha^3*r);
AFS=1/(4*pi*rho*beta^3*r);
for i=1:length(tt)
    MM(i)=integral(@(x)F(x,tt(i)),r/alpha,r/beta);
end
for i=1:length(tt)
    M1(i)=(1-2*(pi*fre*t1(i))^2)*exp(-(pi*fre*t1(i))^2);
    %M1(i)=(t1(i))*exp(-(pi*fre*t1(i))^2);
    M11(i)=(4*(pi*fre*t1(i))^3-6*pi*fre*t1(i))*exp(-(pi*fre*t1(i))^2);
    %M11(i)=(1-2*(pi*fre*t1(i))^2)*exp(-(pi*fre*t1(i))^2);
end
for i=1:length(tt)
    M2(i)=(1-2*(pi*fre*t2(i))^2)*exp(-(pi*fre*t2(i))^2);
    %M2(i)=(t2(i))*exp(-(pi*fre*t2(i))^2);
    M22(i)=(4*(pi*fre*t2(i))^3-6*pi*fre*t2(i))*exp(-(pi*fre*t2(i))^2);
    %M22(i)=(1-2*(pi*fre*t2(i))^2)*exp(-(pi*fre*t2(i))^2);
end
for i=1:3
    for p=1:3
        for q=1:3
            for t=1:length(tt)
                if(split==0)
                     u(i,t)=u(i,t)+AFP*gamma(i)*gamma(p)*gamma(q)*M(p,q)*M11(t)+AFS*(delta(i,p)-gamma(i)*gamma(p))*gamma(q)*M(p,q)*M22(t);
%                         +AN*(15*gamma(i)*gamma(p)*gamma(q)-3*gamma(i)*delta(p,q)-3*gamma(p)*delta(i,q)-3*gamma(q)*delta(i,p))*M(p,q)*MM(t)+...
%                         +AIP*(6*gamma(i)*gamma(p)*gamma(q)-gamma(i)*delta(p,q)-gamma(p)*delta(i,q)-gamma(q)*delta(i,p))*M(p,q)*M1(t)-AIS*(6*gamma(i)*gamma(p)*gamma(q)-gamma(i)*delta(p,q)-gamma(p)*delta(i,q)-2*gamma(q)*delta(i,p))*M(p,q)*M2(t);           
                else                    
                     up(i,t)=AFP*gamma(i)*gamma(p)*gamma(q)*M(p,q)*M11(t);      
                     us(i,t)=AFS*(delta(i,p)-gamma(i)*gamma(p))*gamma(q)*M(p,q)*M22(t);   
                end
            end
        end
    end
end
U=u;
Up=up;
Us=us;
end

