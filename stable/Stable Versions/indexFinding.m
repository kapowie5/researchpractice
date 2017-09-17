function [ Index, IndexFlag ] = indexFinding( theta0 )
%indexFinding
%   Searches through the data in theta and then assigns index starting
%   numbers for each column for a different type of data.  Also saves a
%   flag that tells if that data exists
%   Inputs:     theta0 (matrix)

%   Outputs:    Index (structure)
%                        PT-  indices for getting PT from data file
%                        TC-  indices for getting TC from data file
%                        t-   indices for getting time from data file
%                        Res- indices for getting resistance from data file
%                        SpecStart- indices for where spectra starts in data file
%               IndexFlag (structure)
%                        PT-  flag that PT is in data file
%                        TC-  flag that TC is in data file
%                        t-   flag that time is in data file
%                        Res- flag that Resistance is in data file
%                        SpecStart- flag that Spectra is in data file

%%
%ranges for each potential data thing
Ltime=1e8;  Utime=1e9;
LTC=20;     UTC=40;
LPT=273;    UPT=340;
LRes=1000;  URes=1260;
LSpec=1200; USpec=66000;
IndexFlag=struct('TC',0,'PT',0,'Time',0,'SpecStart',0,'Res',0); 
Index=struct('TC',0,'PT',0,'Time',0,'SpecStart',0,'Res',0); 

for i=1:10
    minTest{i,1}=min(theta0(:,i));
    maxTest{i,1}=max(theta0(:,i));

    if minTest{i,1}==0 && maxTest{i,1}==0
        %case for PT and TC not turned on
        %i-1 is because there usually are two 0s so it'll set it the second
        %time the loop executes
        Index.PT=i-1;  
        Index.TC=i;
        fprintf(horzcat('Zero Case ',num2str(i),'\n'));
        IndexFlag.TC=1; IndexFlag.PT=1;
    elseif minTest{i,1} >=Ltime && maxTest{i,1} <=Utime
        % time is column #i
        Index.Time = i;
        fprintf(horzcat('Time Case ',num2str(i),'\n'));
        IndexFlag.Time=1;
    elseif minTest{i,1} >=LTC && maxTest{i,1} <=UTC
        % TC is column #i
        Index.TC=i;
        fprintf(horzcat('TC Case ',num2str(i),'\n'));
        IndexFlag.TC=1;
    elseif minTest{i,1} >=LPT && maxTest{i,1} <=UPT
        % PT is column #i
        Index.PT=i;
        fprintf(horzcat('PT Case ',num2str(i),'\n'));
        IndexFlag.PT=1;
    elseif minTest{i,1} >=LRes && maxTest{i,1} <=URes
        % Resistance is column #i
        Index.Res=i;
        fprintf(horzcat('Res Case ',num2str(i),'\n'));
        IndexFlag.Res=1;
    elseif minTest{i,1} >=LSpec && maxTest{i,1} <=USpec && Index.SpecStart==0
        % Spectra begins at column #i
        Index.SpecStart=i;
        fprintf(horzcat('Spec Case ',num2str(i),'\n'));
        IndexFlag.SpecStart=1;
    elseif Index.SpecStart ~=0
        %you've found the column for the spec starting
        fprintf(horzcat('Spectra start at ',num2str(Index.SpecStart),'\n'))
    else
        fprintf('Error in importing and sorting data\n'); 
        IndexFlag=struct('TC',0,'PT',0,'Time',0,'SpecStart',0,'Res',0); 
        beep
    end

end

%%

end

