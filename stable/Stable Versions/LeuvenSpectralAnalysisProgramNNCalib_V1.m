% function [works] = LeuvenSpectralAnalysisProgramNNCalib_V1(SpectraData,Folder,FileName,Inputs)
for NNloop = 1:Inputs.NNVariations
 close all
    %% Choose which thetas folder to select
    if Inputs.ThetaChoice==1
        theta = SpectraData.ThetaAll;
    else
        theta = SpectraData.ThetaNorm;
    end


    TempLow=FileName.CalibCell{1,1}(2:4);
    TempHigh=FileName.CalibCell{1,end}(2:4);

    

    %% Reorganize thetas if that choice was made
    if Inputs.ThetaArrange==1
       theta=[theta(:,2:end),theta(:,1)];%% [I,P,Ratio,FWHM,PWL,Norm]
       theta=theta(1:end,:);%% [I,P,Ratio,FWHM,PWL,Norm]
    end

    %% Chose theta columns for Neural network training
    acronyms = {'PI','II','R','F','PW','Norm'};  %Indexed to each thing  
    if Inputs.ThetaChoice==1
        p0 = 'You selected all spectral data.';
        p1 = '1 is Peak Intensity (PI)';
        p2 = '2 is Integrated Intensity (II)';
        p3 = '3 is Ratio (R)';
        p4 = '4 is FWHM (F)';
        p5 = '5 is Peak Wavelength (PW)';
        p6 = '6 to 606 is peak normalized spectra';
        p7 = '607 is Temperature';
        fprintf('%s\r\n %s\r\n %s\r\n %s\r\n %s\r\n %s\r\n %s\r\n ',p0,p1,p2,p3,p4,p5,p6,p7);
        prompt2= 'Which aspects to use [1,2,3,5 are default] ';
        input_Q = Inputs.NNParams{1,NNloop};      
        num_inputs = length(input_Q);
        if num_inputs > 6
            tempindex = [input_Q(input_Q<=5),6];
            extraname = acronyms(tempindex);
        else
            tempindex = input_Q;
            extraname = acronyms(tempindex);
        end
    else
        fprintf('You selected norm sectral data');
        prompt2= 'Select beginning index, ending index, and number of spectra [150,600,40 is default] ';
        input_Q = input(prompt2);
            if isempty(input_Q)
                input_Q=[round(linspace(150,600,40))];
            else
                s=input_Q(1); e=input_Q(2); num=input_Q(3);
                input_Q=[round(linspace(s,e,num))];
            end 

    end
    
    name = [FileName.DateStr,', Temps',TempLow,'K-',TempHigh,'K, Inputs-',strjoin(extraname,',')]

    %% Make folders fo placement of neural networks files
    mkdir(Folder.NNparamsSave,[name,', ',num2str(Inputs.HiddenNodes),' Hidden Nodes'])

    savingFolder=[Folder.NNparamsSave,'\',name,', ',num2str(Inputs.HiddenNodes),' Hidden Nodes'];
    cd(savingFolder)
    copyfile([Folder.Sources,'\train_NN2.m'],savingFolder);
    copyfile([Folder.Sources,'\apply_NN.m'],savingFolder);
    copyfile([Folder.Sources,'\param.nn'],savingFolder);

    % input_Q=[1,2,3];  %%5,80-------------------------------------------
    save input_Q.mat input_Q
    % load input_Q
    factor4train=0.80;%%5,80-------------------------------------------
    SamplingF=50;
    train_Q=[];test_Q=[];
    % for itQ=1:size(theta,1)/SamplingF;
    %     train_Q=[train_Q,(itQ-1)*SamplingF+1:itQ*SamplingF-round((1-factor4train)*SamplingF)];
    %     test_Q=[test_Q,itQ*SamplingF-round((1-factor4train)*SamplingF)+1:SamplingF*itQ];
    % end
    % load train_Q.mat
    % load test_Q.mat
    RandIndex=randperm(size(theta,1));
    train_Q=RandIndex(1:size(theta,1)*factor4train);
    test_Q=RandIndex(size(theta,1)*factor4train+1:end);
    % save train_Q.mat train_Q
    % save test_Q.mat test_Q
    input_train=theta(train_Q,input_Q);
    target_train=theta(train_Q,end);
    input_test=theta(test_Q,input_Q);
    target_test=theta(test_Q,end);

    %% train CG NN 
    cd(Folder.NNparamsLoad)

    input_cg=[input_train,target_train;input_test,target_test];
    target_T=input_cg(:,end);
    cg1=size(input_train,1);
    cg2=size(input_cg,2);
    load param.nn

    cd(savingFolder)

    param([1,2])=[2,2];% configuration for NN
    param([3,4])=[1,cg2-1];     % input colums
    param([5,6])=[cg2,cg2];     % target colums
    param([7,8])=[cg2,cg2];     % output colums
    param(12)=Inputs.NNIterations;               % iteration times
    param([13,14])=[cg1,size(input_cg,1)]; % test rows,from param 13 to 14
    param([15,16])=[Inputs.HiddenNodes,1];       % hidden layers,output units  %%5,80-------------------------------------------
    param(17)=1;       % norm
    save param.nn param -ascii -tabs

    thetam=train_NN2(param,input_cg);
    h1=gcf;
    saveas(h1,[name,'_NN_hidden_nodes.png'])
    saveas(h1,[name,'_NN_hidden_nodes.fig'])
    close
    h2=gcf;
    saveas(h2,[name,'_NN_Iterations.png'])
    saveas(h2,[name,'_NN_Iterations.fig'])

    output_T=thetam(:,end);
    output_train=output_T(1:cg1);
    output_test=output_T(cg1+1:end);
    error_T=output_T-target_T;
    error_train=output_train-target_train;
    error_test=output_test-target_test;
    rms_train=(sum(error_train.*error_train)/length(error_train)).^(0.5)
    rms_test=(sum(error_test.*error_test)/length(error_test)).^(0.5)

    cd(Folder.NNparamsLoad)
    load shift.nn
    load factor.nn

    cd(savingFolder)
    thetam2=apply_NN(param,shift,factor,input_cg(:,1:end-1));
    T_cg(:,1)=target_T;T_cg(:,2)=output_T;T_cg(:,3)=error_T;


    %%
    N_scan=size(theta,1)/SamplingF;
    [sortTrain,idTrain]=sort(target_train); 
    [sortTest,idTest]=sort(target_test);
    T_av_train=ones(N_scan,2);
    T_av_test=ones(N_scan,2);
    for iN=0:N_scan-1
        range_1=(1:SamplingF)+iN*SamplingF;
        T_AV(iN+1,:)=[mean(theta(range_1,end))]; %
    end
    for iscan=1:length(T_AV)

    idTrain0=find(target_train<(T_AV(iscan)+3)&(target_train>=(T_AV(iscan)-3)));
    T_av_train(iscan,1)=mean(target_train(idTrain0));
    T_av_train(iscan,2)=mean(output_train(idTrain0));

    idTest0=find(target_test<(T_AV(iscan)+3)&(target_test>=(T_AV(iscan)-3)));
    T_av_test(iscan,1)=mean(target_test(idTest0));
    T_av_test(iscan,2)=mean(output_test(idTest0));

    end
    error_T_train_av=T_av_train(:,2)-T_av_train(:,1);
    error_T_test_av=T_av_test(:,2)-T_av_train(:,1);

    rms_train_av=(sum(error_T_train_av.*error_T_train_av)/length(error_T_train_av)).^(0.5);
    rms_test_av=(sum(error_T_test_av.*error_T_test_av)/length(error_T_test_av)).^(0.5);

    cd(savingFolder)
    % close all
    figure(33)
    subplot(211)
    plot(target_train,error_train,'ob'); hold on
    plot(target_test,error_test,'.r'); hold on
    % xlim([260,325]);
    ylabel('Absolute error (K)','Fontsize',10);xlabel('Target temperature (K)','Fontsize',10)
    subplot(212)
    saveas(gcf,[name,'_NN_train_error.png'])
    saveas(gcf,[name,'_NN_train_error.fig'])

    plot(target_train,output_train,'ob'); hold on
    plot(target_test,output_test,'.r'); hold on
    name_rms_train_cg=['RMSE of train data: ',num2str(rms_train)];
    name_rms_test_cg=['RMSE of test data: ',num2str(rms_test)];
    text(290,320,name_rms_train_cg);text(290,310,name_rms_test_cg);
    % ylim([215,325])
    % xlim([215,325]);
    ylabel('Output temperature (K)','Fontsize',10);xlabel('Target temperature (K)','Fontsize',10)
    name33=['T_NN'];
    saveas(gcf,[name,'_NN_error_unaveraged.png'])
    saveas(gcf,[name,'_NN_error_unaveraged.fig'])

    figure(44)
    subplot(211)
    plot(T_av_train(:,2),error_T_train_av,'ob'); hold on
    plot(T_av_train(:,1),error_T_test_av,'.r'); hold on
    % xlim([215,325]);
    ylabel('Absolute error (K)','Fontsize',10);xlabel('Target temperature (K)','Fontsize',10);
    subplot(212)
    plot(T_av_train(:,1),T_av_train(:,2),'ob'); hold on
    plot(T_av_test(:,1),T_av_test(:,2),'.r'); hold on
    name_rms_train_cg_av=['RMSE of train data: ',num2str(rms_train_av)];
    name_rms_test_cg_av=['RMSE of test data: ',num2str(rms_test_av)];
    text(290,320,name_rms_train_cg_av);text(290,310,name_rms_test_cg_av);
    % ylim([215,325]);
    % xlim([215,325]);
    ylabel('Output temperature (K)','Fontsize',10);xlabel('Target temperature (K)','Fontsize',10)
    name44=['T_NN_AV'];
    saveas(gcf,[name,'_NN_error_averaged.png'])
    saveas(gcf,[name,'_NN_error_averaged.fig'])
    error_T_train=[target_train,output_train,error_train];
    error_T_test=[target_test,output_test,error_test];
    error_T_trainAV=[T_av_train,error_T_train_av];
    error_T_testAV=[T_av_test,error_T_test_av];
    save error_T_train.txt error_T_train -ascii -tabs
    save error_T_test.txt error_T_test -ascii -tabs
    save error_T_trainAV.txt error_T_trainAV -ascii -tabs
    save error_T_testAV.txt error_T_testAV -ascii -tabs
    
    Folder.NNsavedParams{1,NNloop} = savingFolder;
    Folder.NNsavedParams

    % end
end