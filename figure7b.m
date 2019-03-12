clear
close all
vel=[3000,1731,1800];% velocity of P and S wave, density
thomsen=[0.346,0.28,-0.05];% Thomesn parameters;
dismodel=[325,25,-75];% strike, slip and rake
slope=-90:1:90;
theta1=2/3*pi;
theta2=1*pi;
N=100000;
kesi=rand(N,1);
yita=rand(N,1);
t=(1-kesi)*cos(theta1)+kesi*cos(theta2);
theta=acos(t);
phi=2*pi*yita;
delta=[1,0,0;0,1,0;0,0,1];
kk=length(slope);
ratio=zeros(4,kk);
Rp=zeros(kk,1);
Rs=zeros(kk,1);
Rsh=zeros(kk,1);
Rsv=zeros(kk,1);
aaaa=[0,0,45,90];
bbbb=[0,1,1,1];
for llll=1:4
typeofmedia=bbbb(llll);% isotropic=0, anisotropic=other;
bond=[aaaa(llll),0]; %  Polarization angle and azimuth
for ij=1:kk
    p=zeros(N,1);
    p1=zeros(N,1);
    p2=zeros(N,1);
    p3=zeros(N,1);
    s=zeros(N,1);
    s1=zeros(N,1);
    s2=zeros(N,1);
    s3=zeros(N,1);
    sh=zeros(N,1);
    sv=zeros(N,1);
    sump=0;
    sums=0;
    sumsh=0;
    sumsv=0;
    param=[dismodel(1),dismodel(2),dismodel(3),slope(ij),bond(1),bond(2)];
    [MMM,~]=computemoment(param,thomsen,typeofmedia,vel);
    M0=sqrt(sum(sum(MMM.^2))/2);
    for i=1:N
        gamma1=sin(theta(i))*cos(phi(i));
        gamma2=sin(theta(i))*sin(phi(i));
        gamma3=cos(theta(i));
        gamma=[gamma1,gamma2,gamma3];
        for m=1:3
            for n=1:3
                p1(i)=p1(i)+MMM(m,n)*gamma(1)*gamma(m)*gamma(n);
                p2(i)=p2(i)+MMM(m,n)*gamma(2)*gamma(m)*gamma(n);
                p3(i)=p3(i)+MMM(m,n)*gamma(3)*gamma(m)*gamma(n);
                s1(i)=s1(i)+MMM(m,n)*(delta(1,m)-gamma(1)*gamma(m))*gamma(n);
                s2(i)=s2(i)+MMM(m,n)*(delta(2,m)-gamma(2)*gamma(m))*gamma(n);
                s3(i)=s3(i)+MMM(m,n)*(delta(3,m)-gamma(3)*gamma(m))*gamma(n);
            end
        end
        sv(i)=cos(theta(i))*cos(phi(i))*s1(i)+cos(theta(i))*sin(phi(i))*s2(i)-sin(theta(i))*s3(i);
        sh(i)=-sin(phi(i))*s1(i)+cos(phi(i))*s2(i);
        p(i)=p1(i)^2+p2(i)^2+p3(i)^2;
        sump=sump+p(i);
        s(i)=s1(i)^2+s2(i)^2+s3(i)^2;
        sums=sums+s(i);
        sumsh=sumsh+sh(i)^2;
        sumsv=sumsv+sv(i)^2;
    end
    sump=sump/(N*M0^2);
    sums=sums/(N*M0^2);
    sumsh=sumsh/(N*M0^2);
    sumsv=sumsv/(N*M0^2);
    Rp(ij)=sqrt(sump);
    Rs(ij)=sqrt(sums);
    Rsh(ij)=sqrt(sumsh);
    Rsv(ij)=sqrt(sumsv);
    ratio(llll,ij)=(sums/sump)*3;
end
end
plot(slope,ratio(1,:),'-','LineWidth',1.5)
hold on
plot(slope,ratio(2,:),'--','LineWidth',1.5)
hold on
plot(slope,ratio(3,:),':','LineWidth',1.5)
hold on
plot(slope,ratio(4,:),':.','LineWidth',1.5)
xticks([-90,-60,-30,0,30,60,90])
set(gca,'XTickLabel',{'-90^{\circ}','-60^{\circ}','-30^{\circ}','0^{\circ}','30^{\circ}','60^{\circ}','90^{\circ}'})
xlim([-90 90])
grid on
h1=legend('ISO','VTI','TTI','HTI');
set(h1,'box','off')
xlabel('Slope','FontSize',14,'Fontname','Times New Roman')
ylabel('Es/Ep','FontSize',14,'Fontname','Times New Roman')