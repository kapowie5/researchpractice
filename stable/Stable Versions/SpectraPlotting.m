function [ name ] = SpectraPlotting (SpectraData,WLs,Inputs,Folder,FileName,FileData,IndicesFluoro,IndicesBump,IndicesGrn)

    cd(Folder.FiguresParent)
    name = horzcat(FileName.YearStr,'-',FileData.monthNumStr,'-',FileData.dateStr,', Exp',FileData.ExpStr);
    mkdir(Folder.FiguresParent,name)
    FolderFiguresSave = [Folder.FiguresParent,'\',name];
    cd(FolderFiguresSave)

%% Plot spectra
        h10=figure(10);
        box on
        grid on
        hold on
        plot(WLs.Raw,SpectraData.Raw,'linewidth',2);
        xlabel('wavelength (nm)','fontsize',20)
        ylabel('intensity (counts)','fontsize',20)
        axis('tight') 
        set(gca,'fontsize',15)
        saveas(h10,[name,'_RawSpectra'],'jpg')

        %Smoothed spectra
        h11=figure(11);
        box on
        grid on
        hold on
        plot(WLs.Vis(IndicesFluoro.Spec),SpectraData.Smoothed(:,IndicesFluoro.Spec),'linewidth',2);
        xlabel('wavelength (nm)','fontsize',20)
        ylabel('intensity (counts)','fontsize',20)
        title(['Smoothed Spectra, n=',num2str(Inputs.Nx)])
        axis('tight') 
        set(gca,'fontsize',15)
        saveas(h11,[name,'_SmoothedSpectra'],'jpg')
        saveas(h11,[name,'_SmoothedSpectra'],'fig')

%         %Peak Normalized smoothed spectra
%         h12=figure(12);
%         box on
%         grid on
%         hold on
%         plot(WLs.Vis(IndicesFluoro.Spec),SpectraData.SmoothedNorm(:,IndicesFluoro.Spec),'linewidth',2);
%         xlabel('wavelength (nm)','fontsize',20)
%         ylabel('normalized intensity (a.u.)','fontsize',20)
%         title(['Smoothed, Peak Normalized Spectra, n=',num2str(Inputs.Nx)])
%         axis('tight') 
%         set(gca,'fontsize',15)
%         saveas(h12,[name,'_NormSmoothedSpectra'],'jpg')
%         saveas(h12,[name,'_NormSmoothedSpectra'],'fig')

        % Plot Green spectra
        h13=figure(13);
        box on
        grid on
        hold on
        plot(WLs.Vis(IndicesGrn.Spec),SpectraData.Grn,'linewidth',2);
        xlabel('wavelength (nm)','fontsize',20)
        ylabel('intensity (counts)','fontsize',20)
        axis('tight') 
        set(gca,'fontsize',15)
        saveas(h13,[name,'_GreenSpectra'],'jpg')
        saveas(h13,[name,'_GreenSpectra'],'fig')

        % Plot "IR" spectra
        h14=figure(14);
        box on
        grid on
        hold on
        plot(WLs.IR(IndicesBump.Spec),SpectraData.Bump,'linewidth',2);
        xlabel('wavelength (nm)','fontsize',20)
        ylabel('intensity (counts)','fontsize',20)
        axis('tight') 
        set(gca,'fontsize',15)
        saveas(h14,[name,'_OtherSpectra'],'jpg')
        saveas(h14,[name,'_OtherSpectra'],'fig')

        %% Plot Spectral Characteristics, Red are with respect to Temperature
        % Integrated Intensity vs Temperature
        h20=figure(20);
        grid on
        box on
        hold on
        plot(SpectraData.TempK,SpectraData.IntegIntensity,'-*r','linewidth',3);
        xlabel('Temperature (K)','fontsize',20)
        ylabel('Integrated Intensity (counts)','fontsize',20)
        set(gca,'fontsize',15)
        saveas(h20,[name,'_TempDepenIntI'],'jpg')
        saveas(h20,[name,'_TempDepenIntI'],'fig')

        % Integrated Intensity vs Time
        h21=figure(21);
        grid on
        box on
        hold on
        plot(SpectraData.Time,SpectraData.IntegIntensity,'-*b','linewidth',3);
        xlabel('Time (sec)','fontsize',20) %Really recorded spectra so if 1 spectra/sec, then this is time in seconds
        ylabel('Integrated Intensity (counts)','fontsize',20)
        set(gca,'fontsize',15)
        saveas(h21,[name,'_IntegratedIntensity'],'jpg')
        saveas(h21,[name,'_IntegratedIntensity'],'fig')

        % Peak Intensity vs Temperature--------------------------------------
        h22=figure(22);
        grid on
        box on
        hold on
        plot(SpectraData.TempK,SpectraData.PeakIntensity,'-*r','linewidth',3);
        xlabel('Temperature (K)','fontsize',20)
        ylabel('Peak Intensity (counts)','fontsize',20)
        set(gca,'fontsize',15)
        saveas(h22,[name,'_TempDepenPeakI'],'jpg')
        saveas(h22,[name,'_TempDepenPeakI'],'fig')

        % Peak Intensity vs Time
        h23=figure(23);
        grid on
        box on
        hold on
        plot(SpectraData.Time,SpectraData.PeakIntensity,'-*b','linewidth',3);
        ylabel('Peak Intensity (counts)','fontsize',20)
        xlabel('Time (sec)','fontsize',20) %Really recorded spectra so if 1 spectra/sec, then this is time in seconds
        set(gca,'fontsize',15)
        saveas(h23,[name,'_PeakIntensity'],'jpg')
        saveas(h23,[name,'_PeakIntensity'],'fig')

        % Ratio vs Temp---------------------------------------------------
        h24=figure(24);
        box on
        grid on
        hold on
        plot(SpectraData.TempK,SpectraData.IntegIntensity./SpectraData.PeakIntensity,'-*r','linewidth',3);
        xlabel('Temperature (K)','fontsize',20) 
        ylabel('Ratio (a.u.)','fontsize',20)
        set(gca,'fontsize',15)
        saveas(h24,[name,'_TempDepenRatio'],'jpg')
        saveas(h24,[name,'_TempDepenRatio'],'fig')

        % Ratio vs Time
        h25=figure(25);
        box on
        grid on
        hold on
        plot(SpectraData.Time,SpectraData.IntegIntensity./SpectraData.PeakIntensity,'-*b','linewidth',3);
        xlabel('Time (sec)','fontsize',20) %Really recorded spectra so if 1 spectra/sec, then this is time in seconds
        ylabel('Ratio (a.u.)','fontsize',20)
        set(gca,'fontsize',15)
        saveas(h25,[name,'_Ratio'],'jpg')
        saveas(h25,[name,'_Ratio'],'fig')

        % Fixed spectral position vs Temp-----------------------------------
        h26=figure(26);
        box on
        grid on
        hold on
        plot(SpectraData.TempK,SpectraData.SmoothedNorm(:,600),'-*r','linewidth',3);
        xlabel('Temperature (k)','fontsize',20) 
        ylabel('Normalized Spec Location (a.u.)','fontsize',20)
        set(gca,'fontsize',15)
        saveas(h26,[name,'_TempDepenLoc'],'jpg')
        saveas(h26,[name,'_TempDepenLoc'],'fig')

        % Fixed spectral position vs Time
        h27=figure(27);
        box on
        grid on
        hold on
        plot(SpectraData.Time,SpectraData.SmoothedNorm(:,600),'-*b','linewidth',3);
        xlabel('Time (sec)','fontsize',20) %Really recorded spectra so if 1 spectra/sec, then this is time in seconds
        ylabel('Normalized Spec Location (a.u.)','fontsize',20)
        set(gca,'fontsize',15)
        saveas(h27,[name,'_FixedLocation'],'jpg')
        saveas(h27,[name,'_FixedLocation'],'fig')

        % Emission wavelength vs Temperature ----------------------------------
        h32=figure(32);
        box on
        grid on
        hold on
        plot(SpectraData.TempK,SpectraData.PeakWL,'-*r','linewidth',3);
        xlabel('Temperature (K)','fontsize',20)
        ylabel('Emission Maximum (nm)','fontsize',20)
        set(gca,'fontsize',15)
        saveas(h32,[name,'_TempPeakWavelength'],'jpg')
        saveas(h32,[name,'_TempPeakWavelength'],'fig')

        % Emission wavelength vs Time
        h33=figure(33);
        box on
        grid on
        hold on
        plot(SpectraData.Time,SpectraData.PeakWL,'-*b','linewidth',3);
        xlabel('Time (sec)','fontsize',20) %Really recorded spectra so if 1 spectra/sec, then this is time in seconds
        ylabel('emission maximum (nm)','fontsize',20)
        set(gca,'fontsize',15)
        saveas(h33,[name,'_PeakWavelength'],'jpg')
        saveas(h33,[name,'_PeakWavelength'],'fig')

        % Green Wavelength vs Time
        h34=figure(34);
        box on
        grid on
        hold on
        plot(SpectraData.Time,SpectraData.SumGrn,'-*b','linewidth',3);
        xlabel('Time (sec)','fontsize',20) %Really recorded spectra so if 1 spectra/sec, then this is time in seconds
        ylabel('Green Intensity (counts)','fontsize',20)
        set(gca,'fontsize',15)
        saveas(h34,[name,'_GreenIntensity'],'jpg')
        saveas(h34,[name,'_GreenIntensity'],'fig')

        % IR wavelength vs Time
        h35=figure(35);
        box on
        grid on
        hold on
        plot(SpectraData.Time,SpectraData.SumBump,'-*b','linewidth',3);
        xlabel('Time (sec)','fontsize',20) %Really recorded spectra so if 1 spectra/sec, then this is time in seconds
        ylabel('IR Intensity (counts)','fontsize',20)
        set(gca,'fontsize',15)
        saveas(h35,[name,'_IRIntensity'],'jpg')
        saveas(h35,[name,'_IRIntensity'],'fig')
        
        
        %%
        % Peak and Integ vs signal
        h36=figure(36);
        box on
        grid on
        hold on
        blsigBump=(SpectraData.SumBump-min(SpectraData.SumBump));
        blsigGrn=(SpectraData.SumGrn-min(SpectraData.SumGrn));
        blp=(SpectraData.PeakIntensity-min(SpectraData.PeakIntensity));
        bli=(SpectraData.IntegIntensity-min(SpectraData.IntegIntensity));
        normsigBump=blsigBump/max(blsigBump);
        normsigGrn=blsigGrn/max(blsigGrn);
        normp=blp/max(blp);
        normi=bli/max(bli);
        plot(SpectraData.Time,normsigBump,'-*g','linewidth',3);
        plot(SpectraData.Time,normsigGrn,'-*b','linewidth',3);
        plot(SpectraData.Time,normp,'-*k','linewidth',3);    
        plot(SpectraData.Time,normi,'-*r','linewidth',3);
        xlim([0.2*max(SpectraData.Time),0.4*max(SpectraData.Time)]);
        xlabel('Time (sec)','fontsize',20) %Really recorded spectra so if 1 spectra/sec, then this is time in seconds
        ylabel('Normalized (counts)','fontsize',20)
        legend('IR Signal','Grn Signal','Peak','Integrated')
        set(gca,'fontsize',15)
        saveas(h36,[name,'_ref_and_signal'],'jpg')
        saveas(h36,[name,'_ref_and_signal'],'fig')


        % Temperature vs Time
        h40=figure(40);
        box on
        grid on
        hold on
        plot(SpectraData.Time,SpectraData.TempK,'-*g','linewidth',3);
        xlabel('Time (sec)','fontsize',20) %Really recorded spectra so if 1 spectra/sec, then this is time in seconds
        ylabel('temperature (K)','fontsize',20)
        set(gca,'fontsize',15)
        saveas(h40,[name,'_Temp'],'jpg')
        saveas(h40,[name,'_Temp'],'fig')
        
end