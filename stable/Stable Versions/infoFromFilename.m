function [ FileData, FileDataFlags ] = infoFromFilename( filename )
%infoFromFilename
%   Detailed explanation goes here
%   Inputs:     filename (str) - the filename of the data

%   Outputs:    FileDataFlags (structure) - 1 exists, 0 doesn't
%                        Theta- indices for getting wavelengths from file
%                        Spec-  indices for getting the data from the
%                                spectra file
%               Data     (structure)

    %remove extension from filename string
    s=regexp(filename,'\.txt','split');
    filename=s{1,1};

    %initialize Data Stuctures
    FileData      = struct('Standard','empty','ModFreq',0,'IRint',0,'monthStr','empty', ...
                'monthNum',0','monthNumStr','empty','dateStr','empty','dateNum',0,...
                'flourophore','empty','tempBehavior','empty','modType','empty','voltage',0,...
                'intTime',0,'notes','empty','timeSpace',0,'MotorOld',0,'ExpStr','empty',...
                'Aver',0,'Tfile','empty','Vac','empty','Motor','empty','DATA','empty',...
                'TempK',0,'Calib','empty');
    FileDataFlags = struct('Standard',0,'ModFreq',0,'IRint',0,'monthStr',0, ...
                'monthNum',0','monthNumStr',0,'dateStr',0,'dateNum',0,...
                'flourophore',0,'tempBehavior',0,'modType',0,'voltage',0,...
                'intTime',0,'notes',0,'timeSpace',0,'MotorOld',0,'ExpStr',0,...
                'Aver',0,'Tfile',0,'Vac',0,'Motor',0,'DATA',0,...
                'TempK',0,'Calib',0);
    
    
    %1
      %create flag for nonstandard names and flag
      s=regexp(filename(4),'[0-9]','match');
      if isempty(s)
          FileFileData.Standard='Non-Standard';
          FileDataFlags.Standard=0;
      else
          FileFileData.Standard='Standard';
          FileDataFlags.Standard=1;
      end
      
    
    %2      
    %Modulation frequency in hertz
    s=regexp(filename,'Mod[0-9]*','match');
    if isempty(s)~=1
        FileData.ModFreq=str2double(s{1,1}(4:end))/1000;
        FileDataFlags.ModFreq=1;
    else
        FileData.ModFreq=0;
        FileDataFlags.ModFreq=0;
    end
    
    
    %3
    %IR intensity
    s=regexp(filename,'IR[0-9]*','match');
    if isempty(s)~=1
        FileData.IRint=str2double(s{1,1}(3:end));
        FileDataFlags.IRint=1;
    else
        FileData.IRint=0;
        FileDataFlags.IRint=0;
    end
        

    %4,5,6
    %Month
    FileData.monthStr=filename(1:3);
        if      FileData.monthStr=='Jan'
                FileData.monthNum=1; FileData.monthNumStr='01';
                FileDataFlags.monthStr=1;FileDataFlags.monthNum=1;
                FileDataFlags.monthNumStr=1;
            elseif FileData.monthStr=='Feb'
                FileData.monthNum=2; FileData.monthNumStr='02';
                FileDataFlags.monthStr=1;FileDataFlags.monthNum=1;
                FileDataFlags.monthNumStr=1;
            elseif FileData.monthStr=='Mar'
                FileData.monthNum=3; FileData.monthNumStr='03';
                FileDataFlags.monthStr=1;FileDataFlags.monthNum=1;
                FileDataFlags.monthNumStr=1;
            elseif FileData.monthStr=='Apr'
                FileData.monthNum=4; FileData.monthNumStr='04';
                FileDataFlags.monthStr=1;FileDataFlags.monthNum=1;
                FileDataFlags.monthNumStr=1;
            elseif FileData.monthStr=='May'
                FileData.monthNum=5; FileData.monthNumStr='05';
                FileDataFlags.monthStr=1;FileDataFlags.monthNum=1;
                FileDataFlags.monthNumStr=1;
            elseif FileData.monthStr=='Jun'
                FileData.monthNum=6; FileData.monthNumStr='06';
                FileDataFlags.monthStr=1;FileDataFlags.monthNum=1;
                FileDataFlags.monthNumStr=1;
            elseif FileData.monthStr=='Jul'
                FileData.monthNum=7; FileData.monthNumStr='07';
                FileDataFlags.monthStr=1;FileDataFlags.monthNum=1;
                FileDataFlags.monthNumStr=1;
            elseif FileData.monthStr=='Aug'
                FileData.monthNum=8; FileData.monthNumStr='08';
                FileDataFlags.monthStr=1;FileDataFlags.monthNum=1;
                FileDataFlags.monthNumStr=1;
            elseif FileData.monthStr=='Sep'
                FileData.monthNum=9; FileData.monthNumStr='09';
                FileDataFlags.monthStr=1;FileDataFlags.monthNum=1;
                FileDataFlags.monthNumStr=1;
            elseif FileData.monthStr=='Oct'
                FileData.monthNum=10; FileData.monthNumStr='10';
                FileDataFlags.monthStr=1;FileDataFlags.monthNum=1;
                FileDataFlags.monthNumStr=1;
            elseif FileData.monthStr=='Nov'
                FileData.monthNum=11; FileData.monthNumStr='11';
                FileDataFlags.monthStr=1;FileDataFlags.monthNum=1;
                FileDataFlags.monthNumStr=1;
            elseif FileData.monthStr=='Dec'
                FileData.monthNum=12; FileData.monthNumStr='12';
                FileDataFlags.monthStr=1;FileDataFlags.monthNum=1;
                FileDataFlags.monthNumStr=1;
            else
                FileData.monthNum=0; FileData.monthNumStr='empty';
                FileDataFlags.monthStr=0;FileDataFlags.monthNum=0;
                FileDataFlags.monthNumStr=0;
        end
     
        
    %7,8
    %Date
    s=regexp(filename,[FileData.monthStr,'\w?[0-9]{2}'],'match');
    if isempty(s)~=1
        FileData.dateStr=(s{1,1}(end-1:end));
        FileData.dateNum=str2double(FileData.dateStr);
        FileDataFlags.dateStr=1;
        FileDataFlags.dateNum=1;
    else
        FileData.dateStr='empty';
        FileData.dateNum=0;
        FileDataFlags.dateStr=0;
        FileDataFlags.dateNum=0;
    end
        
    
    %9    
    %Fluorophore
    s=regexp(filename,'(QD|RhB)','match');
    if isempty(s)~=1
        FileData.flourophore=s{1,1};
        FileDataFlags.flourophore=1;
    else
        FileData.flourophore='empty';
        FileDataFlags.flourophore=0;
    end
    
   
    %10
    %Temp Behavior
    s=regexp(filename,'Temp\w*','match');
    if isempty(s)~=1
        FileData.tempBehavior=s{1,1};
        FileDataFlags.tempBehavior=1;
    else
        FileData.tempBehavior='empty';
        FileDataFlags.tempBehavior=0;
    end
    
    
    %11
    %Modulation Type
    s=regexpi(filename,'(Sin|Sqr|CW)','match');
    if isempty(s)~=1
        FileData.modType=s{1,1};
        FileDataFlags.modType=1;
    else
        FileData.modType='Sqr';
        FileDataFlags.modType=0;
    end
     
    
    %12
    %Voltage
    s=regexp(filename,'[0-9]{1,2}(?=V)','match');
    if isempty(s)~=1
        FileData.voltage=str2double(s{1,1});
        FileDataFlags.voltage=1;
    else
        FileData.voltage=0;
        FileDataFlags.voltage=0;
    end
    

    %13
    %Integration Time, in seconds
    % Regular expression looks for a number followed by ms in parenthesis
    % that can include a .
    s3=regexp(filename,'\([0-9]*\.?[0-9]*ms\)','match');
    [rows, columns]=size(s3);
    if isempty(s3)~=1
        s2=regexp(s3{1,1},'[0-9]*\.?[0-9]*','match');
        s3=s3{1,1}(2:end-1);
        FileData.intTime=str2double(s2{1,1})/1000;
        FileDataFlags.intTime=1;
    else
        s3='';
        FileData.intTime=0;
        FileDataFlags.intTime=0;
    end        
        
    
    %14
    %Notes
    % Regular expression looks for anything between parenthesis that isn't
    % already found in integration time (s3).  s2 basically only selects
    % the second set of parenthesis contents
    s=regexp(filename,['\((?!',s3,').*\)'],'match');
    [rows,columns]=size(s);
    
    if isempty(s)~=1
        FileData.notes=s{1,1}(2:end-1);
        FileDataFlags.notes=1;
    else
        FileData.notes='empty';
        FileDataFlags.notes=0;
    end
        
    
    %15
    %TimeSpacing
    % series of regular expressions that look in the notes variable,
    % matches anything that would be similar to seconds, in the case of
    % where the notes also pick up 
    s=regexp(FileData.notes,'[0-9]{1,3}(\s?sec|s\s|ms)','match');
    if isempty(s)~=1
        numS=regexp(s{1,1},'[0-9]+','match');
        strS=regexp(s{1,1},'[0-9]+','split');
            strS=strS{1,end};
        switch strS
            case 'sec'
                FileData.timeSpace=str2double(numS{1,1});
                FileDataFlags.timeSpace=1;
            case 's'
                FileData.timeSpace=str2double(numS{1,1});
                FileDataFlags.timeSpace=1;
            case 'ms'
                FileData.timeSpace=str2double(numS{1,1})/1000;
                FileDataFlags.timeSpace=1;
        end
    else
        FileData.timeSpace=0;
        FileDataFlags.timeSpace=0;
    end
    
    
    %16
    %Motor position
    % Searchs for spot that says motor in string
    s=regexp(filename,'Motor(+|-)[0-9]+um','match');
    if isempty(s)~=1
        s2=s{1,1}(6);
        switch s2
            case '+'
                FileData.MotorOld=str2double(s{1,1}(7:end-2));
                DataFlag.MotorOld=1;
            case '-'
                FileData.MotorOld=str2double(s{1,1}(7:end-2))*(-1);
                DataFlag.MotorOld=1;
        end
    else
        FileData.MotorOld=0;
        FileDataFlags.MotorOld=0;
    end        
        
    
    %17
    %Need experimental number
    %   Date between date#_  -QD
    s=regexp(filename,'_[0-9]{2}(-|_)','match');
    if isempty(s)~=1
        FileData.ExpStr=s{1,1}(2:3);
        FileDataFlags.ExpStr=1;
    else
        FileData.ExpStr='empty';
        FileDataFlags.ExpStr=0;
    end

    
    %18
    %Number of averages
    s=regexpi(filename,'[0-9]*aver','match');
    if isempty(s)~=1
        FileData.Aver=s{1,1}(1:end-4);
        FileDataFlags.Aver=1;
    else
        FileData.Aver=1;
        FileDataFlags.Aver=0;
    end
    
    
    %19
    %Is it a T file from the temperature controller program?
    s=regexpi(filename,'T_file','match');
    if isempty(s)~=1
        FileData.Tfile=s{1,1};
        FileDataFlags.Tfile=1;
    else
        FileData.Tfile='empty';
        FileDataFlags.Tfile=0;
    end
        
    
    %20
    %Is it under vacuum?
    s=regexpi(filename,'(Vac|Air)','match');
    if isempty(s)~=1
        FileData.Vac=s{1,1};
        FileDataFlags.Vac=1;
    else
        FileData.Vac='empty';
        FileDataFlags.Vac=0;
    end
   

    %21
    %Motor position in microns
    s=regexpi(filename,'Motor(+|-)?[0-9]*um','match');
    if isempty(s)~=1
        s2=s{1,1};
        s2=regexp(s2,'[0-9]','match');
        FileData.Motor=s{1,1};
        FileDataFlags.Motor=1;
    else
        FileData.Motor=0;
        FileDataFlags.Motor=0;
    end
    
    
    %22
    %Is it one of the actual good tests with the label of DATA
    s=regexp(filename,'\(DATA\)','match');
    if isempty(s)~=1
        FileData.DATA=s{1,1};
        FileDataFlags.DATA=1;
    else
        FileData.DATA='empty';
        FileDataFlags.DATA=0;
    end
   

    %23
    %Temperature that test occured at
    s=regexp(filename,'T[0-9]{3}K','match');
    if isempty(s)~=1
        s2=regexp(s{1,1},'[0-9]+','match');
        if isempty(s2)~=1
            FileData.TempK=str2double(s2{1,1});
            FileDataFlags.TempK=1;
        else
            FileData.TempK=0;
            FileDataFlags.TempK=0;
        end
    else
        FileData.TempK=0;
        FileDataFlags.TempK=0;
    end
    
    
    %24
    %Flag for calibrated Spectras
    s=regexpi(filename,'calib','match');
    if isempty(s)~=1
        FileData.Calib=FileData.TempK;
        FileDataFlags.Calib=1;
    else
        FileData.Calib=0;
        FileDataFlags.Calib=0;
    end
    
end