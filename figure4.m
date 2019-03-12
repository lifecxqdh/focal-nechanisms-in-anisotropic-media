clear
close all
vel=[3000,1731,1800];% velocity of P and S wave, density
typeofmedia=0;% isotropic=0, anisotropic=other;
thomsen=[0.346,0.28,-0.05];% Thomesn parameters;
bond=[0,0]; %  Polarization angle and azimuth
dismodel=[30,45,60];% strike, slip and rake
slope=[-90,-45,-22.5,-10,0,10,22.5,45,90];
% slope=[-10,0,10];
kk=length(slope);
for ii=1:kk
    param=[dismodel(1),dismodel(2),dismodel(3),slope(ii),bond(1),bond(2)];
    [MMM,~]=computemoment(param,thomsen,typeofmedia,vel);
    x=linspace(0,2*pi,256);
    y=linspace(0,pi,256);
    delta=[1,0,0;0,1,0;0,0,1];
    k=length(x);
    s=length(y);
    p=zeros(k,s);
    p1=zeros(k,s);
    p2=zeros(k,s);
    p3=zeros(k,s);
    s1=zeros(k,s);
    s2=zeros(k,s);
    s3=zeros(k,s);
    ss=zeros(k,s);
    sv=zeros(k,s);
    sh=zeros(k,s);
    m=0;
    [xx,yy]=meshgrid(x,y);
    xx=xx';
    yy=yy';
    r=1;
    txx=r.*sin(xx).*cos(yy);
    tyy=r.*sin(xx).*sin(yy);
    tzz=r.*cos(xx);
    for a=x
        m=m+1;
        n=0;
        for b=y
            n=n+1;
            theta1=a;
            phi1=b;
            gamma1=sin(theta1)*cos(phi1);
            gamma2=sin(theta1)*sin(phi1);
            gamma3=cos(theta1);
            gamma=[gamma1,gamma2,gamma3];
            for i=1:3
                for j=1:3
                    p(m,n)=p(m,n)+MMM(i,j)*gamma(i)*gamma(j)/(4*pi*vel(1)^3*vel(3));
                    p1(m,n)=p1(m,n)+MMM(i,j)*gamma(1)*gamma(i)*gamma(j)/(4*pi*vel(1)^3*vel(3));
                    p2(m,n)=p2(m,n)+MMM(i,j)*gamma(2)*gamma(i)*gamma(j)/(4*pi*vel(1)^3*vel(3));
                    p3(m,n)=p3(m,n)+MMM(i,j)*gamma(3)*gamma(i)*gamma(j)/(4*pi*vel(1)^3*vel(3));
                    s1(m,n)=s1(m,n)+MMM(i,j)*(delta(1,i)-gamma(1)*gamma(i))*gamma(j)/(4*pi*vel(2)^3*vel(3));
                    s2(m,n)=s2(m,n)+MMM(i,j)*(delta(2,i)-gamma(2)*gamma(i))*gamma(j)/(4*pi*vel(2)^3*vel(3));
                    s3(m,n)=s3(m,n)+MMM(i,j)*(delta(3,i)-gamma(3)*gamma(i))*gamma(j)/(4*pi*vel(2)^3*vel(3));
                    sv(m,n)=cos(theta1)*cos(phi1)*s1(m,n)+cos(theta1)*sin(phi1)*s2(m,n)-sin(theta1)*s3(m,n);
                    sh(m,n)=-sin(phi1)*s1(m,n)+cos(phi1)*s2(m,n);
                    ss(m,n)=sqrt(sv(m,n)^2+sh(m,n)^2);
                end
            end
        end
    end
    subplot(3,3,ii)
    surf(p1,p2,p3,p);
    shading interp
    axis image
    axis off
    caxis([-5e-5,5e-5])
    colormap([winter;flipud(autumn)])
    view(-45,0)
%     export_fig(gcf,'-dtif',['ISO',num2str(slope(ii)),'.tif'],'-r300','-transparent')
%   title([num2str(slope(ii)),'^{\circ}'])
end
% export_fig(gcf,'-dtif','ISO2.tif','-r300','-transparent')
% annotation(gcf,'textbox',...
%     [0.069 0.907 0.0839 0.0762],...
%     'String',{'(b)'},...
%     'LineStyle','none',...
%     'FontSize',14,...
%     'FontName','Times New Roman');