% Select a folder with calibrated and modulated spectra
%       Version 4 is meant to reduce memory usage by deleting SpectraDataMod info and working toward slope
%       calculation and diffusivity calculation
%       Version 5 includes NN application for temperature
%%
%   1 - The first part finds all the calibration files which will then each be
%       processed and joined together to make a theta to train the NN with.
%                   DONE
%   2 - Find the green, fluoro, and bump indices based on the lowest
%       temperature and closest to the origin file
%   3 - The NN is trained to based on PeakIntensity, IntegratedIntensity,
%       FWHM, PeakWL, Ratio, and Normalized Locations from theta.  
%   4 - Then the modulated spectra at each temperature will be processed with a
%       theta created for each motor position.
%   5 - The NN is applied to each motor position to get a temperature.
%   6 - The fft of the reference signal is then determined, followed by the fft
%       of each spectra characteristic, followed by the NN temperature, 
%       and then the ln(MAG(T)) and PHASE(T) at the modulation frequency at
%       each position is stored for Peak, Integrated, PeakWL, and NN temperature.
%   7 - The slope of an absolute value function is then fit to the data for the
%       MAG and PHASE.
%   8 - The diffusivity is then calculated as the [-Pi*Freq]/[slope(MAG)xslope(PHASE)]
%   
clc
% close all
clear all

%%
%[ Folder ] = SharedSpectraFolders;
%cd(Folder.Input)
%open('LeuvenSpectralSupportProgramInputParameters.m')
%pause
%[ Inputs ] = LeuvenSpectralSupportProgramInputParameters;

%Folder.SpectraExpHome=uigetdir([Folder.Doc,Folder.mainDir,Folder.fspec]);
%cd(Folder.SpectraExpHome)
%dirs = regexp(genpath(pwd),['[^;]*'],'match');
%for i = 2:length(dirs)
%    Folder.TempCell{1,i-1}=dirs{1,i};
%end
%xstruct=dir;



%% 1 - Calibration files processing
FileName.numCalibCell = length(FileName.CalibCell);
FileName.SpecExpHomeFolder = Folder.SpectraExpHome;
tempAll_1  = [];
tempNorm_1 = [];
[y,Fs] = audioread([Folder.Sound,'\Calibrating Files.mp3']);sound(y,Fs);

for CalibI = 1:FileName.numCalibCell    
    FileName.TxtName = FileName.CalibCell{1,CalibI};
    %Save current FileName without .calib at end
    FileName.FullName = FileName.TxtName(1:end-6);
    %Save full path of current FileName 
    FileName.FullPath = horzcat(FileName.SpecExpHomeFolder,'\',FileName.TxtName);
    %Save Time Stamp of current FileName
    scrap = dir(FileName.FullPath);
    FileName.Time=scrap.date;    
    FileName.TimeNum=scrap.datenum;
%     FileName.DateStr=FileName.Time(1:6);
    FileName.DateStr=[FileName.DateExp(4:5),'-',FileName.DateExp(1:3)];
    FileName.YearStr=scrap.date(8:11);    clear scrap
    cd(Folder.Programs)
    [SpectraData] = LeuvenSpectralAnalysisProgramSpecCalib_V2(Folder,FileName,Inputs,IndicesFluoro,IndicesGrn,IndicesBump,FileData);
    
    tempAll_2  = SpectraData.ThetaAll;
    tempAll_3  = [tempAll_1;tempAll_2];
    tempAll_1  =  tempAll_3;
        
    tempNorm_2  = SpectraData.ThetaNorm;
    tempNorm_3  = [tempNorm_1;tempNorm_2];
    tempNorm_1  =  tempNorm_3;
    
    TempTC(:,CalibI) = SpectraData.TC;
    TempK(:,CalibI)  = SpectraData.TempK;
    TempC(:,CalibI)  = SpectraData.TempC;
    
end
SpectraData.ThetaAll  = tempAll_3;
SpectraData.ThetaNorm = tempNorm_3;

cd(Folder.Temperature)
TempFolderName = [FileData.monthStr,FileData.dateStr];
mkdir(TempFolderName);
cd([Folder.Temperature,'\',TempFolderName])
dlmwrite([TempFolderName,', TC.txt'], TempTC)
dlmwrite([TempFolderName,', TempK.txt'], TempK)
dlmwrite([TempFolderName,', TempC.txt'], TempC)

%% Run NN Calibration Process
[y,Fs] = audioread([Folder.Sound,'\Training Neural Network.mp3']);sound(y,Fs);
cd(Folder.Programs)
LeuvenSpectralAnalysisProgramNNCalib_V1;

%% Process Modulated Spectra
%remove certain SpectraDataMod fields to reduce memory usage
fields = {'BaseValues','Raw','BaseSubtracted','SmoothedNorm','Bump','SumBump',...
    'TC','TempC','TempK','Smoothed','ThetaNorm'};
[y,Fs] = audioread([Folder.Sound,'\Smoothing Spectra.mp3']);
[y2,Fs2] = audioread([Folder.Sound,'\Done Processing Spectra.mp3']);

for FolderLoop = 1:length(Folder.TempCell)
    sound(y,Fs);
    %Counts the number of 
    files = dir(Folder.TempCell{1,FolderLoop}); % equal
    files(strncmp({files.name}, '.', 1)) = []; % new, no exceptions
    
    for FileLoop = 1:length(files)
        cd(Folder.Programs)
        FileName.FullPath = [Folder.TempCell{1,FolderLoop},'\',files(FileLoop,1).name];
        FileName.CurrentFileFullPath = FileName.FullPath;
        
        [extra] = LeuvenSpectralAnalysisProgramSpecModulated_V1...
            (Folder,FileName,Inputs,IndicesFluoro,IndicesGrn,IndicesBump,FolderLoop);
%         extra = rmfield(extra,fields);
        SpectraDataMod(FolderLoop,FileLoop) = extra;
%         Spectra temp too
        clear extra
%         SpectraDataMod = rmfield(SpectraDataMod,fields);
        close all
    end
end
sound(y2,Fs2);
fprintf('            Done Smoothing Spectra')

%%  Apply NN to each motorposition
%       remove certain SpectraDataMod fields to reduce memory usage
% fields = {'ThetaNorm'};   %'ThetaAll',
clear extra
for FolderLoop = 1:length(Folder.TempCell)
    files = dir(Folder.TempCell{1,FolderLoop}); % equal
    files(strncmp({files.name}, '.', 1)) = []; % new, no exceptions
    for FileLoop = 1:length(files)
        SpectraDataMod(FolderLoop,FileLoop).NNTemp={[],[],[],[]};
        SpectraDataMod(FolderLoop,FileLoop).ThetaNorm=[];
    end
end
[y,Fs] = audioread([Folder.Sound,'\Applying Neural Network.mp3']);sound(y,Fs);
for FolderLoop = 1:length(Folder.TempCell)
    %Counts the number of 
    files = dir(Folder.TempCell{1,FolderLoop}); % equal
    files(strncmp({files.name}, '.', 1)) = []; % new, no exceptions
    
    
    for FileLoop = 1:length(files)
                
        cd(Folder.Programs)
        FileName.FullPath = [Folder.TempCell{1,FolderLoop},'\',files(FileLoop,1).name];
        FileName.CurrentFileFullPath = FileName.FullPath;
        
        [extra, NNs] = LeuvenSpectralAnalysisProgramNNApplication_V1...
            (Folder,FileName,IndicesFluoro,SpectraDataMod,FolderLoop,FileLoop);
        SpectraDataMod(FolderLoop,FileLoop) = extra;
        
%         extra = rmfield(extra,fields);
        clear extra
    end
end

%% Perform FFT
[y,Fs] = audioread([Folder.Sound,'\Perform FFT.mp3']);sound(y,Fs);
% Folder loop is temperature indices
for FolderLoop = 1:length(Folder.TempCell)
    %Counts the number of 
    files = dir(Folder.TempCell{1,FolderLoop}); % equal
    files(strncmp({files.name}, '.', 1)) = []; % new, no exceptions
    
    for FileLoop = 1:length(files)
        
        %Temperature storage
        SpectraDataMod(FolderLoop,FileLoop).Temperature = Folder.TempCell{FolderLoop}(end-3:end-1);
        SpectraDataMod(FolderLoop,FileLoop).TempStr = Folder.TempCell{FolderLoop}(end-3:end-1);
        SpectraDataMod(FolderLoop,FileLoop).TempNum = str2double(SpectraDataMod(FolderLoop,FileLoop).TempStr);
        
        %Perform FFT
        SpectraDataMod(FolderLoop,FileLoop).fftRefGrn      = ...
            fft(SpectraDataMod(FolderLoop,FileLoop).SumGrn);
        SpectraDataMod(FolderLoop,FileLoop).fftPeak        = ...
            fft(SpectraDataMod(FolderLoop,FileLoop).PeakIntensity);
        SpectraDataMod(FolderLoop,FileLoop).fftIntegrated   = ...
            fft(SpectraDataMod(FolderLoop,FileLoop).IntegIntensity);
        SpectraDataMod(FolderLoop,FileLoop).fftPeakWL      = ...
            fft(SpectraDataMod(FolderLoop,FileLoop).PeakWL);
        for iNNs = 1:NNs
            SpectraDataMod(FolderLoop,FileLoop).fftNNTemp{1,iNNs}  = ...
            fft(SpectraDataMod(FolderLoop,FileLoop).NNTemp{1,iNNs});
        end
        
                %% spacing of t
                t=SpectraDataMod(FolderLoop,FileLoop).Time;
                for i=2:length(t)
                    dt(i-1)=t(i)-t(i-1);
                end
                N=length(SpectraDataMod(FolderLoop,FileLoop).Time);
                N2=ceil(N/2);
        SpectraDataMod(FolderLoop,FileLoop).dt=mean(dt);
                Fs=1/SpectraDataMod(FolderLoop,FileLoop).dt;
        SpectraDataMod(FolderLoop,FileLoop).fftFrequencies = 0:Fs/length(t):Fs/2;
        
        %Find modulation frequency based on SpectraDataMod.fftRefGrn
            [pks,locs] = findpeaks((abs(SpectraDataMod(FolderLoop,FileLoop).fftRefGrn)),...
                'MINPEAKDISTANCE',10);%,'THRESHOLD',1,'NPEAKS',2)
            [t,b] = max(pks);
            ModFreqIndex = locs(b);
            ModFreq = SpectraDataMod(FolderLoop,FileLoop).fftFrequencies(ModFreqIndex);
        SpectraDataMod(FolderLoop,FileLoop).ModFreqIndex = ModFreqIndex;
        SpectraDataMod(FolderLoop,FileLoop).ModFreq = ModFreq;
        
        
%         %Plot Phase of reference signal, Peak, Integ, PeakWL, and NN
%         figure;  hold all
%         plot(SpectraDataMod(FolderLoop,FileLoop).fftFrequencies,...
%             unwrap(angle(SpectraDataMod(FolderLoop,FileLoop).fftRefGrn(1:N2,1))));
%         plot(SpectraDataMod(FolderLoop,FileLoop).fftFrequencies,...
%             unwrap(angle(SpectraDataMod(FolderLoop,FileLoop).fftPeak(1:N2,1))));
%         plot(SpectraDataMod(FolderLoop,FileLoop).fftFrequencies,...
%             unwrap(angle(SpectraDataMod(FolderLoop,FileLoop).fftIntegrated(1:N2,1))));
%         plot(SpectraDataMod(FolderLoop,FileLoop).fftFrequencies,...
%             unwrap(angle(SpectraDataMod(FolderLoop,FileLoop).fftPeakWL(1:N2,1))));
% %         plot(SpectraDataMod(FolderLoop,FileLoop).fftFrequencies,...
% %             unwrap(angle(SpectraDataMod(FolderLoop,FileLoop).fftNNTemp(1:N2,1))));
%         legend('Reference','Peak','Integ','Peak WL','NN Temp')
%         Pos     = num2str(SpectraDataMod(FolderLoop,FileLoop).MotorPosition);
%         Temp    = SpectraDataMod(FolderLoop,FileLoop).TempStr;
%         TempPos = [Temp,'K, Position ',Pos,' mm'];
%         title(TempPos)        
%         xlabel('Frequencies (Hz)','fontsize',20)
%         ylabel('Phase (rads)','fontsize',20)
%         axis('tight') 
%         h10=gcf;
%         set(gca,'fontsize',15)
%         saveas(h10,['(',FileName.DateStr,')Phase, ',TempPos],'jpg')
%         hold off
%         
        
        
        %Find values at lock-in frequency
        SpectraDataMod(FolderLoop,FileLoop).LockInPeak = ...
            SpectraDataMod(FolderLoop,FileLoop).fftPeak(ModFreqIndex);
        SpectraDataMod(FolderLoop,FileLoop).LockInIntegrated = ...
            SpectraDataMod(FolderLoop,FileLoop).fftIntegrated(ModFreqIndex);
        SpectraDataMod(FolderLoop,FileLoop).LockInPeakWL = ...
            SpectraDataMod(FolderLoop,FileLoop).fftPeakWL(ModFreqIndex);
        for iNNs = 1:NNs
            SpectraDataMod(FolderLoop,FileLoop).LockInNNTemp{1,iNNs}  = ...
            SpectraDataMod(FolderLoop,FileLoop).fftNNTemp{1,iNNs}(ModFreqIndex);
        end

        %MAGnitude
        SpectraDataMod(FolderLoop,FileLoop).LnMagPeak = ...
            log(abs(SpectraDataMod(FolderLoop,FileLoop).LockInPeak));
        SpectraDataMod(FolderLoop,FileLoop).LnMagIntegrated = ...
            log(abs(SpectraDataMod(FolderLoop,FileLoop).LockInIntegrated));
        SpectraDataMod(FolderLoop,FileLoop).LnMagPeakWL = ...
            log(abs(SpectraDataMod(FolderLoop,FileLoop).LockInPeakWL));
        for iNNs = 1:NNs
            SpectraDataMod(FolderLoop,FileLoop).LnMagNNTemp{1,iNNs}  = ...
            log(abs(SpectraDataMod(FolderLoop,FileLoop).LockInNNTemp{1,iNNs}));
        end

        %PHASE
        SpectraDataMod(FolderLoop,FileLoop).PhaseRef = ...
            angle(SpectraDataMod(1,1).fftRefGrn(ModFreqIndex));
        SpectraDataMod(FolderLoop,FileLoop).PhasePeak = ...
            angle(SpectraDataMod(FolderLoop,FileLoop).LockInPeak);
        SpectraDataMod(FolderLoop,FileLoop).PhaseIntegrated = ...
            angle(SpectraDataMod(FolderLoop,FileLoop).LockInIntegrated);
        SpectraDataMod(FolderLoop,FileLoop).PhasePeakWL = ...
            angle(SpectraDataMod(FolderLoop,FileLoop).LockInPeakWL);
        for iNNs = 1:NNs
            SpectraDataMod(FolderLoop,FileLoop).PhaseNNTemp{1,iNNs}  = ...
            angle(SpectraDataMod(FolderLoop,FileLoop).LockInNNTemp{1,iNNs});
        end
    end
end

%% Put into file for fitmain that can be read
for FolderLoop = 1:length(Folder.TempCell)
    %Counts the number of 
    files = dir(Folder.TempCell{1,FolderLoop}); % equal
    files(strncmp({files.name}, '.', 1)) = []; % new, no exceptions
    
    for FileLoop = 1:length(files)
        %Position in m
        z(FileLoop,1)           = SpectraDataMod(FolderLoop,FileLoop).MotorPosition/1000;
        
        omega(FileLoop,1)       = SpectraDataMod(FolderLoop,FileLoop).ModFreq;
        
        MagPeak(FileLoop,1)     = SpectraDataMod(FolderLoop,FileLoop).LnMagPeak;
        MagInteg(FileLoop,1)    = SpectraDataMod(FolderLoop,FileLoop).LnMagIntegrated;
        MagPeakWL(FileLoop,1)   = SpectraDataMod(FolderLoop,FileLoop).LnMagPeakWL;
        for iNNs = 1:NNs
            MagNNTemp(FileLoop,iNNs) = SpectraDataMod(FolderLoop,FileLoop).LnMagNNTemp{1,iNNs};
        end
        
        PhaseRef(FileLoop,1)    = SpectraDataMod(FolderLoop,FileLoop).PhaseRef;
        ref = PhaseRef(FileLoop,1);
        PhasePeak(FileLoop,1)   = ref - SpectraDataMod(FolderLoop,FileLoop).PhasePeak;
        PhaseInteg(FileLoop,1)  = ref - SpectraDataMod(FolderLoop,FileLoop).PhaseIntegrated;
        PhasePeakWL(FileLoop,1) = ref - SpectraDataMod(FolderLoop,FileLoop).PhasePeakWL;
        for iNNs = 1:NNs
            PhaseNNTemp(FileLoop,iNNs) = ref - ...
            SpectraDataMod(FolderLoop,FileLoop).PhaseNNTemp{1,iNNs};
        end
    end
    
    cd(Folder.ForFitMain)
    TempDate=[SpectraDataMod(FolderLoop,FileLoop).Temperature,'K - ',FileData.monthStr,FileData.dateStr];
    MonthDate = [FileData.monthStr,FileData.dateStr];
    
    mkdir(MonthDate)
    cd([Folder.ForFitMain,'\',MonthDate])
    mkdir(TempDate)
    cd([Folder.ForFitMain,'\',MonthDate,'\',TempDate])
    
    [b,i] = sort(z);
    oldz = z;
    z=z(i);
        
    % write position and Magnitude into files for each temperature
    MagPeak=MagPeak(i);
    MagIngeg=MagInteg(i);
    MagPeakWL=MagPeakWL(i);
    dlmwrite([TempDate,', LnMagPeak.txt'], [z,MagPeak])
    dlmwrite([TempDate,', LnMagInteg.txt'], [z,MagInteg])
    dlmwrite([TempDate,', LnMagPeakWL.txt'], [z,MagPeakWL])
    for iNNs = 1:NNs
        MagNNTemp(:,iNNs)=MagNNTemp(i,iNNs);
        dlmwrite([TempDate,', LnMagNNTemp',num2str(iNNs),'.txt'], [z,MagNNTemp(i,iNNs)])
    end
    
    % write position and Phase into files
    PhasePeak=PhasePeak(i);
    PhaseIngeg=PhaseInteg(i);
    PhasePeakWL=PhasePeakWL(i);
    dlmwrite([TempDate,', PhasePeak.txt'], [z,PhasePeak])
    dlmwrite([TempDate,', PhaseInteg.txt'], [z,PhaseInteg])
    dlmwrite([TempDate,', PhasePeakWL.txt'], [z,PhasePeakWL])
    for iNNs = 1:NNs
        PhaseNNTemp(:,iNNs)=PhaseNNTemp(i,iNNs);
        dlmwrite([TempDate,', PhaseNNTemp',num2str(iNNs),'.txt'], [z,PhaseNNTemp(i,iNNs)])
    end
    
    % write position and frequencies
    dlmwrite([TempDate,', Frequencies.txt'], [z,omega])
end

cd(Folder.SpectraExpHome)
filename = [FileName.DateStr,'.mat'];
save(filename,'-v7.3')
aveFreq = mean(omega);
close all





%% Call Fitmain to calculate slope of MAG then PHASE
% Make thing to be saved that goes, mag of peak, integ, WL then phase, then
% temperature, then motor position, then frequency
cd(Folder.ForFitMain)
filename
beep;beep;beep;beep;beep;
[y,Fs] = audioread([Folder.Sound,'\Ready for Fitting.mp3']);sound(y,Fs);
pause
clearvars -except omega filename

%%
params = [1600      1.75    0.0004      11];
fitmain

%%
mLnMag = params(1);
aveFreq = mean(omega);
aveFreq=f;

alpha = abs((pi*aveFreq)/mLnMag^2)

