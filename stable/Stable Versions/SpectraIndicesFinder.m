function [ IndicesVis, IndicesIR ] = SpectraIndicesFinder( spectraColumns, Index )
%SpectraIndicesFinder
%   Searches through the data in theta and then assigns index starting
%   numbers for each column for a different type of data.  Also saves a
%   flag that tells if that data exists
%   Inputs:     spectraColumns (matrix) - Number of columns in the spectra
%                   data file

%   Outputs:    IndicesVis (structure) - 
%                        USB4000- indices for getting wavelengths from file
%                        Spec-  indices for getting the data from the
%                                spectra file
%               IndicesIR     (structure)

switch spectraColumns
    case 1200
        %The most common one (850-2050) or 
        IndicesVis.USB4000=(850:850+1199);
        IndicesIR.USB4000=[];        
    case 1150
        %The very uncommon one (850-2000) or 
        IndicesVis.USB4000=(850:850+1149);
        IndicesIR.USB4000=[];        
    case 800
        %The Beginning one (800-1600)
        IndicesVis.USB4000=(800:800+799);
        IndicesIR.USB4000=[];
    case 1400
        %The most common one (850-2050) with a bump in it (2600-2800)
        IndicesVis.USB4000=(850:850+1199);
        IndicesIR.USB4000=(2600:2600+199);
    case 2000
        %The most common one (850-2050) with All Bumps in it (2200-3000)
        IndicesVis.USB4000=(850:850+1199);
        IndicesIR.USB4000=(2200:2200+799);
    case 1500
        %The most common one (850-2050) with IR in it (3348-3648)
        IndicesVis.USB4000=(850:850+1200);
        IndicesIR.USB4000=(3348:3348+299);
    otherwise
        %The most common one (850-2050) or 
        IndicesVis.USB4000=(850:850+1199);
        IndicesIR.USB4000=[];
end

    IndicesVis.Spec = (1:length(IndicesVis.USB4000));
    IndicesIR.Spec = (IndicesVis.Spec(end)+1 : IndicesVis.Spec(end)+length(IndicesIR.USB4000));
    
    IndicesVis.Theta = IndicesVis.Spec+(Index.SpecStart-1);
    IndicesIR.Theta = IndicesVis.Theta(end)+1;

%     SpectraData.FluoroRange=[Inputs.FluoroIndexMin:Inputs.FluoroIndexMax];%[150:750];
end

