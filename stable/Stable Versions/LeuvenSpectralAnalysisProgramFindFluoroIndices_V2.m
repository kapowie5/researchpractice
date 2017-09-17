%% Load Experimental Spectral and Resistance FileData

FileStructTemp = dir(Folder.TempCell{1,1});

for i = 1:length(FileStructTemp)
    filenamecell{1,i}=FileStructTemp(i,1).name;
    %21
    %Motor position in microns
    s=regexpi(filenamecell{1,i},'Motor(+|-)?[0-9]*um','match');
    if isempty(s)~=1
        s2=s{1,1};
        s2=regexp(s2,'[0-9]*','match');
        FileData.Motor=s2{1,1};
        FileDataFlags.Motor=1;
    else
        FileData.Motor=0;
        FileDataFlags.Motor=0;
    end
    motorPosition(i)=str2double(FileData.Motor);
end

scratch = abs(motorPosition-max(motorPosition));
[~,MotorIndex] = max(scratch);
scratch = filenamecell{1,i};
FileName.FullPath = [Folder.TempCell{1,1},'\',scratch];

%%
%Save current FileName with .txt at end 
    FileName.TxtName = scratch;
    
    %Save current FileName without .txt at end
    FileName.FullName = FileName.TxtName(1:end-4);
    
    %Save Date and experiment number of current FileName
    FileName.DateExp = FileName.FullName(1:8);
    
    %Save Time Stamp of current FileName
    scrap = dir(FileName.FullPath);
    FileName.Time=scrap.date;    
    FileName.TimeNum=scrap.datenum;
    FileName.YearStr=scrap.date(8:11);
    clear scrap

%%
    data0=load(FileName.FullPath);
    theta0=[];
    theta0=[theta0;data0];




    %%
    % Find starting indices for each type of data
    %   Function in Folder.Programs
    [ Index, IndexFlag ] = indexFinding( theta0 );

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
    
    %%
    % Find where green index starts and add the values of all the
    % wavelengths in that range in preparation for fft finding of the
    % modulation frequency
    clear tempsum tempGrnPeak tempGrnPeakIndex tempNormGrnBase tempI tempIgrad tempIfind
    cd(Folder.Programs)
    
    for i=IndicesVis.Spec
        tempsum(i) = sum(SpectraData.Raw(:,i));
    end
    figure;plot(tempsum)
        tempsum=tempsum-min(tempsum);
        [tempGrnPeak,tempGrnPeakIndex] = findpeaks(tempsum,'MINPEAKDISTANCE',10,'THRESHOLD',10,'NPEAKS',1);
        if tempGrnPeakIndex <=20
            p0 = tempGrnPeak; %old values
            pI0 = tempGrnPeakIndex;
            [tempGrnPeak,tempGrnPeakIndex] = findpeak2(1:200,tempsum(1:60),20,p0,pI0);
        end
        tempNormGrnBase = (tempsum-min(tempsum))/tempGrnPeak;
        tempI=find(tempNormGrnBase>0.4);
        tempIgrad=gradient(tempI);
        tempIfind=find(tempIgrad>1);
        
    IndicesGrn.Spec =    tempI(1:tempIfind(1));
    IndicesGrn.USB4000 = IndicesGrn.Spec+IndicesVis.USB4000(1);
    IndicesGrn.Theta =   IndicesGrn.Spec+IndicesVis.Theta(1); 
        
    SpectraData.Grn=SpectraData.Raw(:,IndicesGrn.Spec);
    SpectraData.SumGrn=sum(SpectraData.Grn,2);
%     clear tempsum tempGrnPeak tempGrnPeakIndex tempNormGrnBase tempI tempIgrad tempIfind
    
    
    
    %% 
    % Find where Bumps in IR range index starts and add the values of all the
    % wavelengths in that range in preparation for fft finding of the
    % modulation frequency
    if (isempty(IndicesIR.USB4000)==0)
    clear tempsum tempNormGrnBase tempI tempIgrad tempIfind
    
        for i=IndicesIR.Spec-IndicesVis.Spec(end)
            tempsum(i) = sum(SpectraData.Raw(:,i+IndicesVis.Spec(end)));
        end
        figure;plot(tempsum)
            tempsum=tempsum-min(tempsum);
            [tempBumpPeak,tempBumpPeakIndex] = findpeak(1:length(tempsum),tempsum,20);
            tempNormBumpBase = (tempsum-min(tempsum))/tempBumpPeak;
            tempI=find(tempNormBumpBase>0.4);
            tempIgrad=gradient(tempI);
            tempIfind=find(tempIgrad>1);
            if (isempty(tempIgrad))==0
                tempIfind=length(tempI);
            end
            
        IndicesBump.Spec =    tempI(1:tempIfind(1));
        IndicesBump.USB4000 = IndicesBump.Spec+IndicesIR.USB4000(1);
        IndicesBump.Theta =   IndicesBump.Spec+IndicesIR.Theta(1); 

        SpectraData.Bump=SpectraData.Raw(:,IndicesBump.Spec);
        SpectraData.SumBump=sum(SpectraData.Bump,2);
%         clear tempsum tempBumpPeak tempBumpPeakIndex tempNormBumpBase tempI tempIgrad tempIfind
    end
    
    %%  defining wavelength indices of Fluorescence that used to be Former Spectra_range
    IndicesFluoro.Spec      = (IndicesGrn.Spec(end)+1:IndicesVis.Spec(end));
    IndicesFluoro.Theta     = (IndicesGrn.Theta(end)+1:IndicesVis.Theta(end));
    IndicesFluoro.USB4000   = (IndicesGrn.USB4000(end)+1:IndicesVis.USB4000(end));
    
    
   