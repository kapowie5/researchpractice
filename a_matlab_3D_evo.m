% This program simulates the heat conduction condition in a thin
% film-substrate 2-layer isotropic structure in Cartesian Coordinate
% System.
% A modulated CW laser with a Gaussian distribution on intensity is used.
% The laser has a sine waveform with a specified frequency.
% FROM ZILONG

clear
% close all
clc
tic

P=1;        %Laser power
a=10e-6*sqrt(2);     %Size of the laser which has a Gaussian distribution intensity and round shape %% a=8e-6 for 10X objective, Ti Bulk; a=11e-6 for Cr film
             %%%a=1e-6 for 50X objective, Ti bulk; a=1.2e-6 for 50X TiF,
             %%%with kf=half of bulk value
a1=a;
a2=a;

R_f=0.1;     %Optical reflectivity 
R_s=0.1;
alpha_f=6e8; %and optical absorption coefficient 6e7 for Ti @ 532nm; 8e7 for Cr
alpha_s=0;

h=2e-7;     %Film thickness, for Al it's 120 nm, for other materials it's 200 nm

Q_f0=P*alpha_f*(1-R_f)/(pi*a1*a2);
Q_s0=P*alpha_s*(1-R_f)*(1-R_s)*exp(-alpha_f*h)/(pi*a1*a2);

%Thermal properties of materials
% 
% mat='Cr';
% k_f=15;     %Film(600nm) - 20; Film(200nm) - 15; - or something similar.
% Rho_f=7140;
% c_f=448;

% % % % % % % % % mat='Ti';
% % % % % % % % % Rho_f=4506;
% % % % % % % % % c_f=523.5;
% % % % % % % % % k_f=21.9;

mat='UO2';
Rho_f=10960;
c_f=237;
k_f=8;
% 
% for k_f=10:30;

% diff_f=k_f/Rho_f/c_f;
% effu_f=sqrt(k_f*Rho_f*c_f);

% %Si
% k_s=148;
% Rho_s=2330;
% c_s=700;
% 
% mat='Al';
% k_f=237/5;
% Rho_f=2700;
% c_f=896;

% mat='Au';
% k_f=318/5;
% Rho_f=19300;
% c_f=129;

% mat='W';
% k_f=173;
% Rho_f=19250;
% c_f=132;

% mat='Zr';
% k_f=22.6;
% Rho_f=6520;
% c_f=278;

% % % % % % % % % SiO2
% % % % % % % % % k_s=8;
% % % % % % % % % Rho_s=10960;
% % % % % % % % % c_s=237;

% UO2 bulk
k_s=8;
Rho_s=10960;
c_s=237;

% 
% %Ti bulk
% k_s=k_f;
% Rho_s=Rho_f;
% c_s=c_f;

% % W bulk
% k_s=173;
% Rho_s=19250;
% c_s=132;

%Zr bulk
% k_s=22.6;
% Rho_s=6520;
% c_s=278;

Rth=0;   %Thermal resistance on the interface

% D_f=k_f/Rho_f/c_f
% D_s=k_s/Rho_s/c_s

%%
load('Spatial Alpha.mat')
% yIndex=10;
% xIndex=40;
D_f = alphaSpace(yIndex,xIndex)*1e-6
D_s = alphaSpace(yIndex,xIndex)*1e-6;
k_f = D_f*Rho_f*c_f
k_s = D_s*Rho_s*c_s;


e_f=sqrt(Rho_f*c_f*k_f);
e_s=sqrt(Rho_s*c_s*k_s);

x3=0;
% 
% for ff=5:6

f=logspace(3,6,16);


% dx=Lth/12;
% if ff==2 dx=Lth/20; end
% xmax=dx*400;
% x_reso=xmax/dx;
% x=0:dx:2*xmax;
% x_number=round(length(x)/2);
% 
% Xi_max=.5/dx*2*pi;
% Xi_reso=x_reso;
% Xi_step=Xi_max/Xi_reso;
% Xi_1=-Xi_max:Xi_step:Xi_max;
% Xi_2=-Xi_max:Xi_step:Xi_max;

% dx=Lth/10;
% xmax=Lth*100;
% % x_reso=xmax/dx;
% % x=0:dx:2*xmax;
% x=0:dx:xmax;
% % x_number=round(length(x)/2);
% % 
% Xi_max=.5*1/dx*2*pi;
% % % Xi_max=1/dx*2*pi;
% % Xi_reso=x_reso;
% % Xi_step=Xi_max/Xi_reso;
% Xi_step=1/(xmax)*2*pi;
% Xi_1=-Xi_max:Xi_step:Xi_max;
% Xi_2=-Xi_max:Xi_step:Xi_max;

for cnt=1:length(f)
clear T_ff
    
fre=f(cnt);
omega=2*pi.*fre;  %Modulated laser frequency, vary from 1Hz to 1kHz exponentially

cnt

Lth=sqrt(D_f/pi/fre);
if fre<1e4 dx=Lth/50; xmax=Lth*20;
elseif fre<1e5 dx=Lth/25; xmax=Lth*20;
else dx=Lth/10;xmax=Lth*20;
end
% x_reso=xmax/dx;
% x=0:dx:2*xmax;
x=0:dx:xmax;
% x_number=round(length(x)/2);
% 
%Spatial Wavenumber
Xi_max=.5*1/dx*2*pi;
% % Xi_max=1/dx*2*pi;
% Xi_reso=x_reso;
% Xi_step=Xi_max/Xi_reso;
Xi_step=1/(xmax)*2*pi;
Xi_1=-Xi_max:Xi_step:Xi_max;
Xi_2=-Xi_max:Xi_step:Xi_max;

    for aa=1:length(Xi_1)

%     aa
%     Xi_reso

        Eta_fa=k_f;
        Eta_fb=0;
        Eta_fc=-Xi_1.^2*k_f-Xi_2(aa).^2*k_f-i*omega*Rho_f*c_f;

        Eta_sa=k_s;
        Eta_sb=0;
        Eta_sc=-Xi_1.^2*k_s-Xi_2(aa).^2*k_s-i*omega*Rho_s*c_s;
        
        Eta_f=(-Eta_fb+sqrt(Eta_fb.^2-4.*Eta_fa.*Eta_fc))./(2*Eta_fa);
        Eta_s=(-Eta_sb+sqrt(Eta_sb.^2-4.*Eta_sa.*Eta_sc))./(2*Eta_sa);

        E1=P*alpha_f*(1-R_f)/(2*pi)*exp((-Xi_1.^2*a1^2/4 - Xi_2(aa).^2*a2^2/4));
        F1=P*alpha_s*(1-R_f)*(1-R_s)*exp(-alpha_f*h)/(2*pi)*exp(-Xi_1.^2*a1^2/4 - Xi_2(aa).^2*a2^2/4);

        E2=k_f*alpha_f^2 - Xi_1.^2*k_f - Xi_2(aa).^2*k_f - i*omega*Rho_f*c_f;
        F2=k_s*alpha_s^2 - Xi_1.^2*k_s - Xi_2(aa).^2*k_s - i*omega*Rho_s*c_s;
        
        E=E1./E2;
        F=F1./F2;

        N11= - Eta_f*k_f;
        N12= + Eta_f*k_f;
        N13=0;
        N21=(-Eta_f*k_f).*exp(-Eta_f*h);
        N22=(+Eta_f*k_f).*exp(Eta_f*h);
        N23= + Eta_s*k_s;
        N31=Rth*((+ Eta_f*k_f).*exp(-Eta_f*h))-exp(-Eta_f*h);
        N32=Rth*((- Eta_f*k_f).*exp(Eta_f*h))-exp(Eta_f*h);
        N33=1;

        R1=E.*(+alpha_f*k_f);
        R2=E.*(+alpha_f*k_f)*exp(-alpha_f*h) - F.*(alpha_s*k_s);
        R3=E.*exp(-alpha_f*h).*(1-Rth.*(+alpha_f*k_f))- F;

        DETN=N11.*N22.*N33+N12.*N23.*N31+N13.*N21.*N32-N11.*N23.*N32-N12.*N21.*N33-N13.*N22.*N31;
        DETA=R1.*N22.*N33+N12.*N23.*R3+N13.*R2.*N32-R1.*N23.*N32-N12.*R2.*N33-N13.*N22.*R3;
        DETB=N11.*R2.*N33+R1.*N23.*N31+N13.*N21.*R3-N11.*N23.*R3-R1.*N21.*N33-N13.*R2.*N31;
%         DETC=N11.*N22.*R3+N12.*R2.*N31+R1.*N21.*N32-N11.*R2.*N32-N12.*N21.*R3-R1.*N22.*N31;
        
        A=DETA./DETN;
        B=DETB./DETN;
%         C=DETC./DETN;
      
%         T_ff(fe,1:length(Xi_1),aa)=A.*exp(-Eta_f.*x3)+B.*exp(Eta_f.*x3)+E.*exp(-alpha_f.*x3);
        T_ff(1:length(Xi_1),aa)=A+B+E;  %Complex Temp
%         T_sf(fe,1:length(Xi_1),aa)=C.*exp(-Eta_s.*(x3-h))+F.*exp(-alpha_s.*(x3-h));
    end

T_frsimp=sum(sum(T_ff))/Xi_step/Xi_step;
phasim(cnt)=atan2(imag(T_frsimp),real(T_frsimp))*180/pi;
Temp(cnt)=T_frsimp;
T_re(cnt)=real(T_frsimp);
T_im(cnt)=imag(T_frsimp);
% phasim

[s1,s2]=size(T_ff);
T_mid=abs(T_ff(:,round(s2/2)));
% T_mid_phase(cnt)=angle(T_ff(101,round(s2/2)));

% q0=k_s/(h*k_f);
% q=sqrt(i*omega/D_f+q0^2/2*(1-sqrt(1-4*i*omega/(D_s*q0^2)*(1-D_s/D_f))));
% 
% x_minus=-xmax:dx:-dx;                         % going to "build up" the minus half of the phase/amplitude symmetrically
% 
% for pc=1:length(x_minus)
%     phase_minus(pc)=phase(length(phase)-pc+1);
%     Tp_minus(pc)=Temp(1,length(Temp)-pc+1);
% end
% 
% x_total=horzcat(x_minus,x);
% phase_total=horzcat(phase_minus,phase');
% 
eval(['figure(',num2str(cnt+1e3),')'])
% figure(18);
% subplot(2,1,1);
plot(T_mid,'r.');
% subplot(2,1,2);
% plot(x_total(length(x)-100:length(x)+100),phase_total(length(x)-100:length(x)+100),'b.');
% hold on
% plot(x_total(length(x)-100:length(x)+100),real(q).*x_total(length(x)-100:length(x)+100)*180/pi,'r');
% plot(x_total(length(x)-100:length(x)+100),-sqrt(omega/(2*D_s)).*x_total(length(x)-100:length(x)+100)*180/pi,'g')
% hold off
% 
title(['frequency = ',num2str(fre)]);
% 
% pause(2)

% clear T_ff

end

phasim = phasim - 180;

% freq=1/4:1/4:9;
% 

figure(107);
% % axesset = axes('Xscale','Log','Yscale','lin');
% % hold(axesset,'all');
% 
semilogx(f,abs(Temp),f,T_re,'r',f,T_im,'g','linewidth',1.5);
title('Amplitude VS Frequency');
xlabel('Log-frequency');
ylabel('Temperature Amplitude');
grid on
% hold off

figure(208);
semilogx(f,phasim,'linewidth',1.5);
title('Phase Lag VS Frequency');
xlabel('Log-frequency');
ylabel('Phase Lag');
grid on


toc

ML3D=[f' phasim' abs(Temp)' (real(Temp))' (imag(Temp))'];
T_f = abs(Temp);
% filename=['Matlab3D_',mat,num2str(ff),'-',num2str(ff+1),'.txt'];
% filename=['Matlab3D_',mat,'.txt'];
% filename=['a=30_k=',num2str(k_f),'.txt'];
filename=['ML3D@h=',num2str(h*1e9),'nm.txt'];
save(filename,'ML3D','-ascii');
save('ML3D.mat','phasim','Temp','T_f','f','D_s','D_f')

% % %% 1D
% n_f=sqrt(2*pi*i.*f/D_f); 
% n_s=sqrt(2*pi*i.*f/D_s);
% Q1D=1e-3./pi/a^2/2;
% A1D=Q1D/4.*(1 + k_s.*n_s*Rth - k_s.*n_s./n_f/k_f).*exp(-n_f*h)./(k_s.*n_s.*cosh(h.*n_f) + k_f.*n_f.*sinh(h.*n_f).*(1 + Rth.*k_s.*n_s));
% T1D=Q1D/2/k_f./n_f+2.*A1D;
% 
% ang1=angle(T1D)/pi*180;
% 
% savedata=[f' ang1' abs(T1D)' real(T1D)' imag(T1D)'];
% save('Matlab1D.txt','savedata','-ascii');

% end
% figure(208)
% hold on
% semilogx(f,ang1,'b.')
% hold off

%%
% end
a_Gaussian_Spectra_Freq