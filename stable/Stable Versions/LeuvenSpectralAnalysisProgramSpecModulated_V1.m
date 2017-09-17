function [SpectraData] = LeuvenSpectralAnalysisProgramSpecModulated_V1(Folder,FileName,Inputs,IndicesFluoro,IndicesGrn,IndicesBump,FolderLoop)

    %Save current FileName with .txt at end
    s=regexp(FileName.CurrentFileFullPath,' Spectras\','split');
    s2=regexp(s{1,end},'\','split');
    
    FileName.TxtName = s2{1,end}; clear s s2
    
    %Save current FileName without .txt at end
    FileName.FullName = FileName.TxtName(1:end-4);
        
    %Save Date and experiment number of current FileName
    FileName.DateExp = FileName.FullName(1:8);
    
    %Save Time Stamp of current FileName
    scrap = dir(FileName.CurrentFileFullPath);
    FileName.Time=scrap.date;    
    FileName.TimeNum=scrap.datenum;
    FileName.YearStr=scrap.date(8:11);
    clear scrap
    
    data0=load(FileName.CurrentFileFullPath);
    theta0=[];
    theta0=[theta0;data0];

    

    %%
    % Find starting indices for each type of data
    %   Function in Folder.Programs
    [ Index, IndexFlag ] = indexFinding( theta0 );
    Index.Time      = 1;
    Index.SpecStart = 2;

    %%
    % Find out how many non-spectra columns by selecting the last number in the string, 
    %   assumes that it's not a multiple of 10
    [SpectraData.NumSpec,columns]=size(theta0);
    lastnum=num2str(columns);
    lastnum=str2num(lastnum(end));

    % Find out how many spectra columns
    %       spectraColumns=columns - lastnum;
    %       spectraColumns=columns - mod(columns,10);
    SpectraData.NumWLs = columns - (Index.SpecStart-1);

    % Call case structure function to find out the wavelength ranges 
    %   of the spectra (in Folder.Programs)
    [IndicesVis,IndicesIR] = SpectraIndicesFinder(SpectraData.NumWLs,Index);
    
    
    % Assign WLs in Fluoro and IR range from the wl file
    cd(Folder.Sources);
    wavelengths=load('wl-3648.txt');   
    WLs.Vis  = wavelengths(1,IndicesVis.USB4000);
    WLs.IR   = wavelengths(1,IndicesIR.USB4000);
    WLs.Raw  = [WLs.Vis,WLs.IR];


    %% 
    % Call Function that gets information about experiment from the name of
    % the file
    % 	Function in Folder.Programs
    cd(Folder.Programs)
    [ FileData, FileDataFlags ] = infoFromFilename( FileName.TxtName );
    
    % Clean up naming conventions
    if FileDataFlags.Standard==1 && FileDataFlags.Tfile==0
        FileName.DateExp=[FileData.monthStr,FileData.dateStr,'_',FileData.ExpStr];
    end

    %%
    % Import SpectraData from theta0
    SpectraData.Raw=theta0(:,Index.SpecStart:end);
    
    SpectraData.Grn=SpectraData.Raw(:,IndicesGrn.Spec);
    SpectraData.SumGrn=sum(SpectraData.Grn,2);
    
    SpectraData.Bump=SpectraData.Raw(:,IndicesBump.Spec);
        SpectraData.SumBump=sum(SpectraData.Bump,2);
    
    %% Put Temperature Data into SpectraData
    if IndexFlag.TC==1
        SpectraData.TC  = theta0(:,Index.TC);
    else
        SpectraData.TC  = zeros(SpectraData.NumSpec);
    end

    if IndexFlag.Res==1 && IndexFlag.PT==1
        Res=theta0(:,Index.Res);
        SpectraData.TempK=Inputs.a.*Res.^2 + Inputs.b.*Res + Inputs.c + 273.15;
        SpectraData.TempC=SpectraData.TempK-273.15;
    elseif IndexFlag.PT==1
        if min(theta0(1,Index.PT))>150
            SpectraData.TempK=theta0(:,Index.PT);
            SpectraData.TempC=SpectraData.TempK - 273.15;
        else
            SpectraData.TempC=theta0(:,Index.PT);
            SpectraData.TempK=SpectraData.TempC + 273.15;
        end
    else
        SpectraData.TempK   = zeros(SpectraData.NumSpec);
        SpectraData.TempC   = zeros(SpectraData.NumSpec);
    end


    %% Put Time Data into SpectraData
    if IndexFlag.Time==1
        SpectraData.Time = (theta0(:,Index.Time) - theta0(1,Index.Time))/1000;
    elseif (FileDataFlags.timeSpace || FileDataFlags.intTime) == 1
        dt=max([FileData.timeSpace,FileData.intTime]);
        SpectraData.Time = dt*linspace(0,SpectraData.NumSpec-1,SpectraData.NumSpec);
    else
        SpectraData.Time = Inputs.dt*linspace(0,SpectraData.NumSpec-1,SpectraData.NumSpec);
    end 
        
        
    %% Subtract baseline, smooth spectra 
    %Subtracting baseline of spectra based on minimal value
    SpectraData.BaseValues=min(SpectraData.Raw,[],2);
    SpectraData.BaseSubtracted=SpectraData.Raw-SpectraData.BaseValues*ones(1,size(SpectraData.Raw,2));

    %Smooth spectra
    fprintf(horzcat('Currently Smoothing Spectra','\n'));
    for ii=1:SpectraData.NumSpec
        SpectraData.Smoothed(ii,:)=filter_mean(SpectraData.BaseSubtracted(ii,1:end),Inputs.Nx); %
    end
    
    %% Calculate peak intensity, integrated intensity,fwhm, emission maximum
    % Call to get findpeak,FWHM programs
    cd(Folder.Sources)

    % Loop through each spectra collected
    for isp=1:SpectraData.NumSpec;
        %At each temperature, the peak wavelength and the intensity at that
        %spectra are found
        [SpectraData.PeakIntensity(isp,:),SpectraData.PeakWL(isp,:)]=...
            findpeak(WLs.Vis(IndicesFluoro.Spec),...
            SpectraData.Smoothed(isp,IndicesFluoro.Spec),30);
        tempPeak=SpectraData.PeakIntensity(isp,:);
        SpectraData.IntegIntensity(isp,:)=...
            sum(SpectraData.Smoothed(isp,IndicesFluoro.Spec));
        SpectraData.Ratio(isp,:)=...
            SpectraData.IntegIntensity(isp,:)/SpectraData.PeakIntensity(isp,:);
        SpectraData.SmoothedNorm(isp,:)=...
            SpectraData.Smoothed(isp,:)./SpectraData.PeakIntensity(isp);
        
        %to do FWHM, check if Peak value is above a Threshold set in the input parameters program
        if (tempPeak>=Inputs.PeakThreshold)
            SpectraData.FWHM(isp,:)=...
                fwhm(WLs.Vis(IndicesFluoro.Spec),...
                SpectraData.Smoothed(isp,IndicesFluoro.Spec));
        else
            SpectraData.FWHM(isp,:)=0;
        end
    end

    SpectraData.TempStr = Folder.TempCell{FolderLoop}(end-3:end-1);

    

    %% create thetas for NN and save smoothed spectral characteristics
    % Temperature is put at the end
    SpectraData.ThetaAll=[SpectraData.PeakIntensity,SpectraData.IntegIntensity,...
        SpectraData.Ratio,SpectraData.FWHM,SpectraData.PeakWL,...
        SpectraData.SmoothedNorm(:,IndicesFluoro.Spec),SpectraData.TempK];
    SpectraData.ThetaNorm=[SpectraData.SmoothedNorm(:,IndicesFluoro.Spec),SpectraData.TempK];
    
    cd(Folder.Programs)
    if Inputs.ModSaveFiles ==1
        fprintf(horzcat('Currently Saving Text Files','\n'));
        [ OutFile ] = SaveSpectralData ( SpectraData, Folder, FileName );   
    end
    
    
    %% Directory for Figures
    cd(Folder.Programs)
    if Inputs.ModGraphSelect==1
        [name] = SpectraPlotting(SpectraData, WLs, Inputs, Folder, FileName, FileData...
            ,IndicesFluoro,IndicesBump,IndicesGrn)
    end

    %% Documentation of Results, date taken, and program version used
    counter = 0;

    PeakWLMean = mean(SpectraData.PeakWL);  counter=counter+1;
    PeakWLStd = std(SpectraData.PeakWL);    counter=counter+1;
    PeakWLMin = min(SpectraData.PeakWL);    counter=counter+1;
    PeakWLMax = max(SpectraData.PeakWL);    counter=counter+1;

    TempMean = mean (SpectraData.TempK);    counter=counter+1;
    TempStd = std(SpectraData.TempK);       counter=counter+1;
    TempMin = min(SpectraData.TempK);       counter=counter+1;
    TempMax = max(SpectraData.TempK);       counter=counter+1;

    IPeakMean = mean(SpectraData.PeakIntensity); counter=counter+1;
    IPeakStd = std(SpectraData.PeakIntensity);   counter=counter+1;
    IPeakMin = min(SpectraData.PeakIntensity);   counter=counter+1;
    IPeakMax = max(SpectraData.PeakIntensity);   counter=counter+1;

    IIntegMean = mean(SpectraData.IntegIntensity); counter=counter+1;
    IIntegStd = std(SpectraData.IntegIntensity);   counter=counter+1;
    IIntegMin = min(SpectraData.IntegIntensity);   counter=counter+1;
    IIntegMax = max(SpectraData.IntegIntensity);   counter=counter+1;

    RatioMean = mean(SpectraData.Ratio); counter=counter+1;
    RatioStd = std(SpectraData.Ratio);   counter=counter+1;
    RatioMin = min(SpectraData.Ratio);   counter=counter+1;
    RatioMax = max(SpectraData.Ratio);   counter=counter+1;

    version1 = mfilename;                            counter=counter+1;

    timestamp = strcat(datestr(clock,'mmmm dd, yyyy HH:MM:SS.FFF AM'),...
    'm',datestr(clock,'ss'),'s');                    counter=counter+1;

    FileName.FullName; counter=counter+1;
    Inputs.Nx; counter=counter+1;

    dataSum = cell(1,counter);
    dataSum{1,1}=timestamp; dataSum{1,2}=FileName.FullName; dataSum{1,3}=FileName.Time; dataSum{1,4}=version1; dataSum{1,5}=Inputs.Nx; 
    dataSum{1,6}=TempMean; dataSum{1,7}=TempStd; dataSum{1,8}=TempMin; dataSum{1,9}=TempMax; 
    dataSum{1,10}=IPeakMean; dataSum{1,11}=IPeakStd; dataSum{1,12}=IPeakMin; dataSum{1,13}=IPeakMax; 
    dataSum{1,14}=IIntegMean; dataSum{1,15}=IIntegStd; dataSum{1,16}=IIntegMin; dataSum{1,17}=IIntegMax; 
    dataSum{1,18}=PeakWLMean; dataSum{1,19}=PeakWLStd; dataSum{1,20}=PeakWLMin; dataSum{1,21}=PeakWLMax; 
    dataSum{1,22}=RatioMean; dataSum{1,23}=RatioStd; dataSum{1,24}=RatioMin; dataSum{1,25}=RatioMax;

    %%
    if FileDataFlags.Motor==1
        s =regexp(FileData.Motor,'[0-9]*','match');
        s2=regexp(FileData.Motor,'(+|-)','match');
        if s2{1,1}=='+'
            %Position in mm
            SpectraData.MotorPosition=str2double(s{1,1})*0.001;
        elseif s2{1,1}=='-'
            %Position in mm
            SpectraData.MotorPosition=str2double(s{1,1})*-0.001;
        else
            SpectraData.MotorPosition=9999;
        end
    else
        SpectraData.MotorPosition=9999;
    end
    
    cd(Folder.Documentation)
    excelFile = 'Spectra_summaries.xls';
    %Title string
    % titlestr = ['Timestamp , Filename , TimeFileDataTaken , Version , SmoothingNx, TempMean , TempStd , TempMin , TempMax , PeakIntMean , PeakIntStd , PeakIntMin , PeakIntMax , IntegIntMean , IntegIntStd , IntegIntMin , IntegIntMin , PeakWLMean , PeakWLStd , PeakWLMin , PeakWLMax , SpectraData.RatioMean , SpectraData.RatioStd , SpectraData.RatioMin , SpectraData.RatioMax']
    [success,message] = xlsappend(excelFile,dataSum,'Sheet1');


    FileName.FullName
    cd(Folder.Programs)
    beep
    
end