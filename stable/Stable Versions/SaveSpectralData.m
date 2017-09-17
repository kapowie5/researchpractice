function [ OutFile ] = SaveSpectralData ( SpectraData, Folder, FileName );

    %Theta all has all spectral information and properties
    cd(Folder.ThetasAll);temp=SpectraData.ThetaAll;
    OutFile.All=horzcat(FileName.DateExp,'_theta_all');
    save([OutFile.All,'.nn'],'temp','-ascii','-tabs');clear temp

    %Theta all norm has only normalized spectra
    cd(Folder.ThetasNorm);temp=SpectraData.ThetaNorm;
    OutFile.Norm=horzcat(FileName.DateExp,'_theta_norm');
    save([OutFile.Norm,'.nn'],'temp','-ascii','-tabs');clear temp

    %Green part of the spectra
    cd(Folder.GrnSignal);temp=SpectraData.Grn;
    OutFile.Grn=horzcat(FileName.DateExp,'_green_signal');
    save([OutFile.Grn,'.spec'],'temp','-ascii','-tabs');clear temp

    %Bumps part of the spectra
    cd(Folder.IRSignal);temp=SpectraData.Bump;
    OutFile.IR=horzcat(FileName.DateExp,'_IR_signal');
    save([OutFile.IR,'.spec'],'temp','-ascii','-tabs');clear temp

    %Integrated Intensity
    cd(Folder.IntegratedIntensity);temp=SpectraData.IntegIntensity;
    OutFile.IntIntensity=horzcat(FileName.DateExp,'_IntIntensity');
    save([OutFile.IntIntensity,'.specdata'],'temp','-ascii','-tabs');clear temp

    %Peak Intensity
    cd(Folder.PeakIntensity);temp=SpectraData.PeakIntensity;
    OutFile.PeakIntensity=horzcat(FileName.DateExp,'_PeakIntensity');
    save([OutFile.PeakIntensity,'.specdata'],'temp','-ascii','-tabs');clear temp

    %Emission Wavelength
    cd(Folder.EmissionWL);temp=SpectraData.PeakWL;
    OutFile.EmissionWL=horzcat(FileName.DateExp,'_EmissionWL');
    save([OutFile.EmissionWL,'.specdata'],'temp','-ascii','-tabs');clear temp

    %Temperature
    cd(Folder.PtTemperature);temp=SpectraData.TempK;
    OutFile.PtTemperature=horzcat(FileName.DateExp,'_PtTemperature');
    save([OutFile.PtTemperature,'.specdata'],'temp','-ascii','-tabs');clear temp

    %Time
    cd(Folder.Time);temp=SpectraData.Time;
    OutFile.Time=horzcat(FileName.DateExp,'_Time');
    save([OutFile.Time,'.specdata'],'temp','-ascii','-tabs');clear temp
    
end