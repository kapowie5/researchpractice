function [ Folder ] = SharedSpectraFolders( ~ )
%SharedSpectraFolders 
%   Selects where the dropbox folder is automatically selected (kinda)
%       Used in spectra_from_liwangV4_ui.m, TFscanCalibCGNN_V2.m,
%       applyNN2SingleTest_V2.m
%   OUTPUT -   Folder (structure)
laptop='C:\Users\colto\Desktop\research practice\repo\stable\dropbox\';
laptopDoc = 'C:\Users\colto\Desktop\research practice\repo\stable\documents\';
desktop='C:\Users\User\Dropbox\';
desktopDoc = 'C:\Users\User\Documents\';
    A=exist(laptop,'dir');
    B=exist(desktop,'dir');
    C=exist(laptopDoc,'dir');
    D=exist(desktopDoc,'dir');
    
    if A==7
        Dboxpath=laptop;
    elseif B==7
        Dboxpath=desktop;
    else
        msgbox('Dropbox Folder doesn''t exist')
        Dboxpath=uigetdir('C:\Users\');
    end
    
    if C==7
        Docpath=laptopDoc;
    elseif D==7
        Docpath=desktopDoc;
    else
        msgbox('Document Folder doesn''t exist')
        Docpath=uigetdir('C:\Users\');
    end

Folder.Dpth = Dboxpath;
Folder.Doc = Docpath;
Folder.mainDir = '3 - Belgium Spider Silk\';
Folder.fspec= 'Fluoresence spectra\';

Folder.Sources= [Folder.Dpth,Folder.mainDir,Folder.fspec,'Sources'];
Folder.Input= [Folder.Dpth,Folder.mainDir,Folder.fspec,'Input Files'];
Folder.Documentation = [Folder.Dpth,Folder.mainDir,Folder.fspec,'Results'];

Folder.ThetasNorm = [Folder.Dpth,Folder.mainDir,Folder.fspec,'Theta_norm'];
Folder.ThetasAll = [Folder.Dpth,Folder.mainDir,Folder.fspec,'Theta_all'];
Folder.GrnSignal = [Folder.Dpth,Folder.mainDir,Folder.fspec,'Green_signal'];
Folder.IRSignal = [Folder.Dpth,Folder.mainDir,Folder.fspec,'IR_signal'];

Folder.IntegratedIntensity = [Folder.Dpth,Folder.mainDir,Folder.fspec,'IntegratedIntensity'];
Folder.PeakIntensity = [Folder.Dpth,Folder.mainDir,Folder.fspec,'PeakIntensity'];
Folder.EmissionWL = [Folder.Dpth,Folder.mainDir,Folder.fspec,'EmissionWL'];
Folder.PtTemperature = [Folder.Dpth,Folder.mainDir,Folder.fspec,'PtTemperature'];
Folder.Time = [Folder.Dpth,Folder.mainDir,Folder.fspec,'Time'];

Folder.FiguresParent = [Folder.Dpth,Folder.mainDir,Folder.fspec,'Spectra Figures'];
Folder.NNFiguresParent = [Folder.Dpth,Folder.mainDir,Folder.fspec,'NN Temp Figures'];
Folder.NNparamsLoad = [Folder.Dpth,Folder.mainDir,Folder.fspec,'NN parameters'];
Folder.NNparamsSave = [Folder.Dpth,Folder.mainDir,Folder.fspec,'NN parameters'];

Folder.Programs = [Folder.Dpth,Folder.mainDir,'Programs\Troy''s spectra ones'];
Folder.Sound =[Folder.Dpth,Folder.mainDir,Folder.fspec,'Sound Files'];

Folder.ForFitMain = [Folder.Dpth,Folder.mainDir,Folder.fspec,'ForFitMain'];

Folder.Temperature = [Folder.Dpth,Folder.mainDir,Folder.fspec,'Temperature'];
Folder.Temperature = [Folder.Dpth,Folder.mainDir,Folder.fspec,'Phase Figures'];
        
end