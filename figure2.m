clear
close all
vel=[3000,1731,1800];% velocity of P and S wave, density
typeofmedia=1;% isotropic=0, anisotropic=other;
thomsen=[0.346,0.28,-0.05];% Thomesn parameters;
bond=[0,0]; %  Polarization angle and azimuth
dismodel=[30,45,60];% strike, slip and rake
slope=-90:1:90;
kk=length(slope);
I=eye(3);
HISO=zeros(kk,1);
HDC=zeros(kk,1);
HCLVD=zeros(kk,1);
aaaa=[0,0,45,90];
bbbb=[0,1,1,1];
for llll=1:4
    typeofmedia=bbbb(llll);% isotropic=0, anisotropic=other;
    bond=[aaaa(llll),0]; %  Polarization angle and azimuth
for ii=1:kk
    param=[dismodel(1),dismodel(2),dismodel(3),slope(ii),bond(1),bond(2)];
    [MMM,~]=computemoment(param,thomsen,typeofmedia,vel);
    [~,D]=eig(MMM);
    MISO=1/3*trace(D)*I;
    Mstar=D-MISO;
    d=eig(abs(Mstar));
    mstarmax=max(d);
    dd=eig(Mstar);
    [~,index]=min(abs(dd));
    bmin=dd(index);
    kesi=-bmin/(mstarmax);
    HISO(ii)=1/3*(trace(D)/max(max(abs(D))));
    HCLVD(ii)=(2*kesi*(1-abs(HISO(ii))));
    HDC(ii)=1-abs(HISO(ii))-abs(HCLVD(ii));
end 
subplot(2,2,llll)
plot(slope,HISO,'r','LineWidth',2)
hold on
plot(slope,HDC,'b:','LineWidth',2)
hold on
plot(slope,HCLVD,'k--','LineWidth',2)
xticks([-90,-60,-30,0,30,60,90])
set(gca,'XTickLabel',{'-90^{\circ}','-60^{\circ}','-30^{\circ}','0^{\circ}','30^{\circ}','60^{\circ}','90^{\circ}'})
xlabel('Slope','FontSize',14,'Fontname','Times New Roman')
ylabel('Percentage','FontSize',14,'Fontname','Times New Roman')
h1=legend('ISO','DC','CLVD');
set(h1,'Location','SouthEast','box','off')
ylim([-1 1])
xlim([-90 90])
grid on
% text(-85,0.85,'','FontSize',18,'Fontname','Times New Roman')
end
annotation(gcf,'textbox',...
    [0.065 0.941 0.035 0.037],...
    'String',{'(a)'},...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
annotation(gcf,'textbox',...
    [0.524 0.941 0.035 0.037],'String',{'(b)'},...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
annotation(gcf,'textbox',...
    [0.065 0.459 0.035 0.037],...
    'String',{'(c)'},...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
annotation(gcf,'textbox',...
    [0.518 0.459 0.035 0.037],...
    'String',{'(d)'},...
    'LineStyle','none',...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
set(gcf,'position',[607   175   922   602])
% export_fig(gcf,'-dtif','momentdec.tif','-r300','-transparent')
% export_fig(gcf,'-dfig','momentdec.fig','-r300','-transparent')