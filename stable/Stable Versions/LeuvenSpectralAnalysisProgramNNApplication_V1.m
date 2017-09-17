function [SpectraDataOut, NNs] = LeuvenSpectralAnalysisProgramNNApplication_V1...
    (Folder,FileName,IndicesFluoro,SpectraDataMod,FolderLoop,FileLoop);

close all

%% Number of neural networks that will be applied
NNs = length(Folder.NNsavedParams);
SpectraDataOut = SpectraDataMod(FolderLoop,FileLoop);
MotorStr = num2str(SpectraDataOut.MotorPosition*1000);
TempStr = SpectraDataOut.TempStr;

name0=[FileName.DateStr,'-',TempStr,'K, ',MotorStr,'um'];

%% Loop through different 
for iNNs=1:NNs;

    clear Time FWHM Ratio theta output_exp1 thetam_exp1
    clear SpectPeak ExpSpect Exp_data Peak_WL ExpSpectNorm
    clear ExpDataIP ExpDataNorm ExpSpect PumpInteg PumpPeak
    clear T_recon T_ACDC T_acdc
    clear exp_data Time PumpLaser ExpSpect PumpInteg PumpPeak SpectInteg
    clear SpectPeak1 SpectPI SpectIP ExpSpect_sort IS FWHM_s ExpSpectNorm
    clear Peak_WL SpectPeak3 T_out ExpData ExpDataNorm ExpDataNorm2 TimeDepen
    clear T_recon fft_T T_ACDC T_fit
    
    theta=SpectraDataMod(FolderLoop,FileLoop).ThetaAll;
    Time=SpectraDataMod(FolderLoop,FileLoop).Time;

    %Spectra Plotting
    h100=figure(100+iNNs);
    plot(theta(10, 6+IndicesFluoro.Theta-IndicesFluoro.Theta(1)))
    cd(Folder.NNFiguresParent)
    title(['Spectra, ',name0])
    
    saveas(h100,[name0,', Spectra'],'jpg')
    saveas(h100,[name0,', Spectra'],'fig') 

    %% choose proper NN
    FolderNN=Folder.NNsavedParams{1,iNNs};
    

    cd(FolderNN) 
    load input_Q.mat
    input_train1=theta(:,input_Q);
    load param.nn
    shiftt=load('shift.nn');
    factorr=load('factor.nn');
    thetam_exp1=apply_NN(param,shiftt,factorr,input_train1);
    output_exp1=thetam_exp1(:,end);
    output_exp1=output_exp1';
    SpectraDataOut.NNTemp{1,iNNs} = output_exp1;

end

cd(Folder.NNFiguresParent)
h7=figure(20+iNNs);
hold all
for iNNs = 1:NNs
    plot(Time,SpectraDataOut.NNTemp{1,iNNs})
    s1 = regexp(Folder.NNsavedParams{1,iNNs},'Inputs.*','match');
    legendCell{1,iNNs} = s1{1,1};
end
hold off
ylabel('Predicted Temperature (K)','fontsize',15)
xlabel('Time (Sec)','Fontsize',15)
title(['NN Temp, ',name0])
legend(legendCell)
saveas(h7,[name0,',NN Temp'],'jpg')
saveas(h7,[name0,',NN Temp'],'fig')    

end