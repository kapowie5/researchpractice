%%  Used for simulation data
% V7 is meant for simulating spectra with typical uncertainties, and a
%       phase and amplitude decay based on a frequency scan
%V7_1 fixes the time spacing issue
% close all
clear all
clc
uncertainty_upload  %Current directory
% uncertainty_upload_laptop

%% Checking which source has 1% error
[val,index] = max(uVectorNums);
errorStr = uVectorNames{index};
if val < sum(uVectorNums)
    errorStr = 'uAll'
end



%% This section defines linear relationships as a function of temperature for the follow spectral features
%       -Peak Intensity (PI), Peak Wavelength (PWL), Full Width Half
%        Maximum (FWHM), and Baseline (0).
%These terms are given a letter to relate to the Gaussian distribution
%parameters
%       -Peak Intensity is A, Peak Wavelength is B, Full Width Half
%        Maximum is C, and Baseline is D.
%These linear relationships take the form A(T) = A_0 +S_A*T, where A is the
%       feature of interest as a function of temperature (T), A_0 is the
%       intercept, S_A is the slope or sensitivity of that feature to
%       temperature, and T is temperature.

%Spectral Features Intercepts - Based on Experimental Data
A_PI_0  = 68800;    %Intercept for Peak Intensity
B_PWL_0 = 602;      %Intercept for Peak Wavelength
C_FWHM_0= 9;        %Intercept for FWHM
D_0     = 1300;        %Intercept for Baseline

%Spectral Features Slopes - Loosely Based on Experimental Data
S_A     = -194*0.8;     %Slope for Peak Intensity
S_B     = 0.1;          %Slope for Peak Wavelength
S_C     = 0.1;          %Slope for FWHM
S_D     = 0;            %Slope for Baseline

%% This section defines the linear relationship of the amplitude and phase of the cosine wave as a function 
%       of axial position, z, based on the thermal diffusivity of the
%       fiber.  Again, these are linear and are of the form: 
%       Mag(z) = m_mag*|z| + b_mag, where m_mag is the slope and b_mag is
%           the intercept of the magnitude of the wave.
%       Phase(z) = m_p*|z| + b_p, where m_p is the slope and b_p is
%           the intercept of the phase of the wave.


% % Load frequency Decay Data
% %This was because the f_vector had lots of 0s in it
% load('21-May-2015 MultiFreq Results, f=0.01-20.01.mat')
% % load('C:\Users\Troy\Dropbox\3 - Belgium Spider Silk\Results\Simulation Results\MultiFrequency Results\21-May-2015 MultiFreq Results, f=0.01-20.01.mat')
% f_vector(f_vector==0)=[];
% % f_vector = f_vector(18:3:51);
% f_vector = f_vector(21:101);
% T_f = reshape(Ts_f(76,3,:),[],1);
% T_f(T_f==0) = [];
% % T_f = T_f(18:3:51,1);
% T_f = T_f(21:101,1);
% disp('Loading data')
% figure;semilogy(f_vector,abs(T_f))

% Load frequency Decay Data
load('ML3D.mat')
% load('C:\Users\Troy\Dropbox\3 - Belgium Spider Silk\Results\Simulation Results\MultiFrequency Results\21-May-2015 MultiFreq Results, f=0.01-20.01.mat')
f_vector=f;
T_f = Temp;
disp('Loading data')
figure;semilogy(f_vector,abs(T_f))

%% This section defines various variables for the simulation
WL_raw  = load('wl-3648.txt');    %Vector of length nWL, of the wavelengths that will be simulated
WL      = [WL_raw(1,850:850+1399) WL_raw(1,2600:2600+99)];
WL_index= [1:length(WL)];
nWL     = length(WL);
nz      = 101;  %Number of axial locations (this number should be odd, so that there is a z=0 value
z_left  = -2e-3;    %Beginning point of the axial distance z, in m
z_right = 2e-3;     %Ending point of the axial distance z, in m
z       = linspace(z_left,z_right,nz);  %Vector of size nz of all axial locations in the simulation
Index_z0= ceil(nz/2);                %Vector index of centerline of z
z_0     = z(Index_z0);                %The centerline value of z (should be 0)
freq    = 10;      %Frequency of the modulation, given in Hz
nF      = length(f_vector);     %number of frequencies from simulation
% nt      = 200;      %Number of time steps that simulation will run
% n_f     = 100;        %Number of periods at f_Hz that the simulation experiences
% t       = linspace(0,f_Hz*n_f,nt);      %Vector of size nt of all times in the simulation, given in seconds
dt      = 0.040;
n_f     = 100;        %Number of periods at f_Hz that the simulation experiences
t       = 0:dt:n_f;
for it=1:length(t)
    t(it) = t(it) + ((randn/16)*ut);
end
t(1)    = 0;
t0      = 0;        %Initial time of heating (usually 0)
T0_0    = 300;      %Intial temperature of the simulation before heating
n       = 0.5;      %Polynomial exponent to simulate DC heating rise
%These next variables are created to make the system look like what we
%       expect it to
S_cos   = 1;        %Sensitivity of the modulation  (increase the size of the temperature variation)
T_Rise  = [26,24,22,20,18,16,14,12,10,8,6,4,2];        %Temperature rise


%%
% This loop will go through all the axial locations.  The reason for doing
%       this is that the temperature signal should decay based on this
%     for iz = 1:nz
%         %This "if statement" is to make sure that the data looks like an
%         %   inverted absolute value function.
%         if z(iz) <= 0
%             A_mag = exp(b_mag+z(iz)*m_mag);
%             phase = m_p*z(iz) + b_p;
%         else
%             A_mag = exp(b_mag+z(iz)*-m_mag);
%             phase = -m_p*z(iz) + b_p;
%         end
%         z_mag(iz)= A_mag;       %Meant to save the data so it can be graphed and checked
%         z_phase(iz)= phase;     %Meant to save the data so it can be graphed and checked
%     end




%% Get z position
clear z
z(1)    = 0.0e-3;% + ((randn/16)*uZ);  %Vector of size nz of all axial locations in the simulation
Index_z0= 1;
nz      = 1;

    
%% Add in frequency effect
for iiF = 1:nF
    f_mag(iiF) = abs(T_f(iiF));
    f_phase(iiF) = degtorad(phasim(iiF));
end
% f_phase = unwrap(f_phase);
f_Amax = max(f_mag);


%% Create Files
% [ Folder ] = SharedSpectraFolders;
Folder.Current = pwd;
% Folder.SpectraSaveParent=[Folder.Doc,Folder.mainDir,Folder.fspec];
Folder.SpectraSaveParent=Folder.Current;
str2 = '2017-05-24';
% str2 = input('Today (YYYY-MM-DD)','s');
str = [str2,', Freq Calibrated Simulated Spectras, 1% error at ',errorStr];
str3 = 'May 24';
% str3 = input('Today (May27)','s');
zstr = num2str(z*1e6);
tstr = num2str(dt*1e3);
cd(Folder.SpectraSaveParent)
mkdir(str);


%% Simulation section
D_0 = 1300;
for iT = 3:3
cd([Folder.SpectraSaveParent,'\',str])
dummy=[errorStr,'.txt']; dlmwrite(dummy,val);
    T0 = (iT-1)*6+290;    
%     Tinf    = T0 + T_Rise;
    Tinf    = T0 + T_Rise(iT);      %Final DC temperature of simulation after heating
                
            nF = length(f_vector);
%             nF=3;
            zstr = num2str(z(1)*1e6);
            for iiF = 1:nF                   
                freq = f_vector(iiF);
                fstr = num2str(freq);  %*1e3
                
                clear t
                dt      = (1/freq)/20;
                dt_vector(iiF) = dt;
                n_f     = 100;        %Number of periods at f_Hz that the simulation experiences
                t_end   = dt*100*20;
                t       = 0:dt:t_end;
                t_vector(iiF,:) = t;
                for it=1:length(t)
                    t(it) = t(it) + ((randn/16)*ut);
                end
                t(1)    = 0;
                t0      = 0;        %Initial time of heating (usually 0)
                nt = length(t);  
                
                    for jt = 1:nt
                        %Variable T is meant to basically superimpose the AC on the DC
                        %component
                        %The DC signal is probably following something like a Gaussian (z-z0) with a width increasing like sqrt(w0+sqrt(alpha/pi/f))
%                         T_DC(iiF,jt) = (Tinf + ((randn/8)*uTfit)-T0+ ((randn/8)*uT)) * sqrt((t(jt)-t0)/t(end)) + (T0+((randn/8)*uT));
                        T_DC(iiF,jt) = T0;
                        T_mod(iiF,jt) = 20*(f_mag(iiF)/f_Amax) * cos(2*pi*(freq+((randn/128)*uFreq))*(t(jt)-t0) - f_phase(iiF));
                        T(iiF,jt) =  T_mod(iiF,jt) + T_DC(iiF,jt);   
                        Ref(iiF,jt) = cos(2*pi*freq*(t(jt)-t0)); 
                        
                        
                        A_PI    = (A_PI_0+((randn/8)*uPIb)) + (S_A+((randn/8)*uPIm))*(T(iiF,jt)+ ((randn/8)*uT));
                        B_PWL   = (B_PWL_0+((randn/8)*uPWLb)) + (S_B+((randn/8)*uPWLm))*(T(iiF,jt)+ ((randn/8)*uT));
                        C_FWHM  = (C_FWHM_0+((randn/8)*uFWHMb)) + (S_C+((randn/8)*uFWHMm))*(T(iiF,jt)+ ((randn/8)*uT));
                        D       = (D_0+((randn/8)*uDb)) + (S_D+((randn/8)*uDm))*(T(iiF,jt)+ ((randn/8)*uT));     
                        
                        %create gaussian spectra
                        for kWL = 1:nWL
                            Spec(iiF,jt,kWL) = (A_PI*exp(- ((WL(kWL)+((randn/16)*uWL)-B_PWL).^2) *(4*log(2)) / (C_FWHM)^2))+((randn/8)*uSN) + D + uNoise;
                        end
                        
                        %Create Green and ~780 nm Reference, then add to Spectra
                        Spec_temp(1,1,:) = ((maketrap(WL_index,1275,1285,1290,1300)) * ((Ref(iiF,jt)+1)/2))*3000;% + uNoise;
%                         Spec_temp(1,1,:) = ((maketrap(WL_index,4,10,45,55) + maketrap(WL_index,1300,1310,1345,1355)) * ((Ref(iiF,jt)+1)/2))*1500;% + uNoise;
%                         Spec_temp(1,1,:) = ((maketrap(WL_index,[4 10 45 55]) + maketrap(WL_index,[1300 1310 1345 1355])) * ((Ref(iiF,jt)+1)/2))*1500;% + uNoise;
%                         Spec(iiF,jt,4:55) = 5000; % Makes a clean reference signal
                        Spec(iiF,jt,:) = Spec(iiF,jt,:) + Spec_temp(1,1,:) + D;    
                        Ref_check(jt) = Spec(iiF,jt,30);                            
                    end
                    
                    t_vector(iiF,:) = t;
                    
%                         if iiF ==number
%                             fftTime = fft(t);
%                             fftT_AC = fft(T_mod(iiF,:));
%                             fftT_Mod = fft(T(iiF,:));
%                             
%                             
%                             for i=2:length(t)
%                                 dt(i-1)=t(i)-t(i-1);
%                             end
%                             N=length(t);
%                             N2=ceil(N/2);
%                             dt=mean(dt);
%                             Fs=1/dt;
%                             fftFrequencies = 0:Fs/length(t):Fs/2;
%                         end
%                         fftFreq_vector(iiF,fftFrequencies);
%                 
            
% Uncomment this if you want to save files
            cd([Folder.SpectraSaveParent,'\',str])            
            tempName = [str3,'_01_T',num2str(T0),'K_Motor+0',zstr,'um-QD(',tstr,'ms),SinMod',fstr,'(Temp).ref'];
            TempSave = [t(1:nt)' Ref(iiF,1:nt)'];
            dlmwrite(tempName,TempSave,'precision', 16)
            clear TempSave
                        
            tempName = [str3,'_01_T',num2str(T0),'K_Motor+0',zstr,'um-QD(',tstr,'ms),SinMod',fstr,'(Temp).temp'];
            TempSave = [t(1:nt)' T(iiF,1:nt)'];
            dlmwrite(tempName,TempSave, 'precision', 16)
            clear TempSave
            
            if iiF == 1
                TemperFile=['T',num2str(T0),'K'];
                mkdir(TemperFile)
                disp(freq)
            else
                disp(freq)
            end
            cd(TemperFile)
            
            spectraName = [str3,'_01_T',num2str(T0),'K_Motor+0',zstr,'um-QD(',tstr,'ms),SinMod',fstr,'_IR105_Vac_1Aver(DATA).txt'];
            SpecSave = [t(1:nt)' squeeze(Spec(iiF,:,:));];
            dlmwrite(spectraName,SpecSave, 'precision', 16);
               

            end
    T0
end            
            

%%
% close all
for i=1:nF
    % spacing of t for fft frequencies
    t = t_vector(i,:);
    for ii=2:length(t)
        junk(ii-1)=t(ii)-t(ii-1);
    end
    N=length(t);
    N2=ceil(N/2);
    dt=mean(junk);
    Fs=1/dt;
    fftFreq(i,:) = 0:Fs/length(t):Fs/2;
%     pause
    
    for j=1:2001
        Peak(i,j) = max(Spec(i,j,:));
        II(i,j) = 0;
        for k=250:850
            II(i,j) = Spec(i,j,k) + II(i,j);
        end
    end
    
    % fft
    fftjunk = fft(Peak(i,:))/N;
        fftPeak(i,:) = fftjunk(1:ceil(N/2));
    fftjunk = fft(II(i,:))/N;
        fftII(i,:) = fftjunk(1:ceil(N/2));
%     Mag(i) = abs(fftPeak(i,101));
%     Phase(i) = angle(fftPeak(i,101));
    Mag(i) = abs(fftII(i,101));
    Phase(i) = angle(fftII(i,101));
    fftjunk = fft(Ref(i,:))/N;
        fftRef(i,:) = fftjunk(1:ceil(N/2));
    PhaseRef(i) = angle(fftRef(i,101));

        
    
end

% figure;plot(f_vector,Mag)
%%
fsize = 16;
figure;loglog(f_vector,abs(T_f)/max(abs(T_f)),f_vector,Mag/max(Mag),'o')
title('Magnitude of Frequency Scan','FontWeight','bold','FontSize', fsize)
xlabel('Frequency (Hz)','FontWeight','bold','FontSize', fsize)
ylabel('Ln(Mag) (a.u.)','FontWeight','bold','FontSize', fsize)
set(gca,'FontWeight','bold','FontSize', fsize)
legend('Model','Simulation')

figure;semilogx(f_vector,(angle(T_f)-pi)*(180/pi),f_vector,(unwrap(PhaseRef-Phase)-pi)*(180/pi),'o')
title('Phase of Frequency Scan','FontWeight','bold','FontSize', fsize)
xlabel('Frequency (Hz)','FontWeight','bold','FontSize', fsize)
ylabel('Phase (rad)','FontWeight','bold','FontSize', fsize)
set(gca,'FontWeight','bold','FontSize', fsize)
legend('Model','Simulation')

figure;
loglog(fftFreq(1,:),abs(fftRef(1,:)')/100)
hold all
for i = 1:16
    loglog(fftFreq(i,:),abs(fftRef(i,:)')/100)
end
for i = 1:16
    loglog(fftFreq(i,:),abs(fftII(i,:)'))
end
title('FFT Magnitude','FontWeight','bold','FontSize', fsize)
xlabel('Frequency (Hz)','FontWeight','bold','FontSize', fsize)
ylabel('Mag (a.u.)','FontWeight','bold','FontSize', fsize)
set(gca,'FontSize',16,'FontWeight','bold','XMinorTick','on','XScale',...
    'log','XTick',[1000 10000 100000 1000000],'YMinorTick','on','YScale','log');
annotation(gcf,'textbox',[0.605 0.417 0.227 0.0657],'String',{'Reference Signal'},...
    'LineStyle','none','FontWeight','bold','FontSize', fsize);
annotation(gcf,'textbox',[0.587 0.774 0.247 0.0657],'String',{'Integrated Intensity'},...
    'LineStyle','none','FontWeight','bold','FontSize', fsize);
xlim([900,1e6])
ylim([1e-6,5e5])

figure;plot(Phase)
hold all
plot(PhaseRef)
title('Phase','FontWeight','bold','FontSize', fsize)
ylabel('Phase ()','FontWeight','bold','FontSize', fsize)
set(gca,'FontWeight','bold','FontSize', fsize)
legend('Ref','Simulation','Location','SouthEast')

figure;semilogy(Mag)
title('Mag','FontWeight','bold','FontSize', fsize)
ylabel('Mag ()','FontWeight','bold','FontSize', fsize)
set(gca,'FontWeight','bold','FontSize', fsize)
legend('Ref','Simulation','Location','SouthEast')

%%
figure;
hold all
for iiF=1:nWL
    test(iiF)=Spec(Index_z0,1,iiF);
end
plot(WL,test)
for iiF=1:nWL
    test(iiF)=Spec(Index_z0,2,iiF);
end
plot(WL,test)
for iiF=1:nWL
    test(iiF)=Spec(Index_z0,5,iiF);
end
plot(WL,test)
for iiF=1:nWL
    test(iiF)=Spec(Index_z0,10,iiF);
end
plot(WL,test)
for iiF=1:nWL
    test(iiF)=Spec(Index_z0,50,iiF);
end
plot(WL,test)
for iiF=1:nWL
    test(iiF)=Spec(Index_z0,100,iiF);
end
plot(WL,test)
xlabel('Wavelength (nm)','FontWeight','bold','FontSize', fsize)
ylabel('Counts','FontWeight','bold','FontSize', fsize)
set(gca,'FontWeight','bold','FontSize',11)
saveas(gcf,'Fig - Select Spectras','png')
saveas(gcf,'Fig - Select Spectras','fig')

%%
figure;
temp=squeeze(Spec(1,:,:));
plot(WL,temp/1000)
% clear temp
xlabel('Wavelength (nm)','FontWeight','bold','FontSize', fsize)
ylabel('Counts (thousands)','FontWeight','bold','FontSize', fsize)
set(gca,'FontWeight','bold','FontSize', fsize)
xlim([550,800])
saveas(gcf,'Fig - All Spectras','png')
saveas(gcf,'Fig - All Spectras','fig')


%%
figure;plot(t,Spec(1,:,600))
xlabel('Time (s)','FontWeight','bold','FontSize', fsize)
ylabel('Counts','FontWeight','bold','FontSize', fsize)
title(['Intensity vs t at wavelength, ',num2str(WL(600)),'nm'])
saveas(gcf,'Fig - I vs t','png')
saveas(gcf,'Fig - I vs t','fig')

%%
figure;plot(t,Spec(1,:,30))
xlabel('Time (s)','FontWeight','bold','FontSize', fsize)
ylabel('Counts','FontWeight','bold','FontSize', fsize)
title('Intensity vs t for Green Reference')
saveas(gcf,'Fig - Ref','png')
saveas(gcf,'Fig - Ref','fig')

%%
figure;plot(t,T(1,:))
xlabel('Time (s)','FontWeight','bold','FontSize', fsize)
ylabel('Temp (K)','FontWeight','bold','FontSize', fsize)
title('Temperature')
saveas(gcf,'Temperature','png')
saveas(gcf,'Temperature','fig')

%%
cd(Folder.Current)
errorStr
beep
beep
beep

% %%
% % cd('C:\Users\Troy\Dropbox\3 - Belgium Spider Silk\Programs\Troy''s spectra ones\')
% cd('C:\Users\User\Dropbox\3 - Belgium Spider Silk\Programs\Troy''s spectra ones\')
% LeuvenSpectralAnalysisProgram_V6freq_good
