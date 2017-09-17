%%%randomize w1, w2 weigth matrices before each iteration??? check rr
function P=train_NN2(param,inputs)
% paramFile='D:\liwang liu\PhD study\experiments\Fluore_exp_NN\withPump\param.nn';
cc=fix(clock);
cc1=fix(clock);
cc2=fix(clock);

units=2;%=input('Give number of hidden units ');
outunit=1;%=input('Press 1 for linear output unit, 2 for tangent hyperbolic ');
norm=1;%=input('Press 1 to normalize signals ');
red=0;%=input('Press 1 for extra signal size reduction ');
stop_crit=1e-20;
%------ inladen van i2/2 amplitudes en fazes per profiel en opp k-waarde,itr data
%i=input('Press 1 to load scaled file theta.nn, 2 to load file using scaling files, 3 to load theta0.nn ','s');
theta=inputs;
[N,kolommen]=size(theta);
i1=1;
i2=1;
i3=2;
i4=2;
i5=2;
i6=2;
itr=round(N/2);
maxiter=102;
inprof='knik';
outprof='knik';
PHI=theta(:,i1:i2)';
scaling=0;
if scaling==0,
    m=mean(PHI');% gemiddelde per signaal over profielen
    sd=sqrt(dot(PHI',PHI')/N-m.*m);% standaard deviatie per signaal over profielen
    a=min(theta(:,i3:i4));
    b=max(theta(:,i3:i4));
    if b-a~=0,
        factor=2*ones(1,i4-i3+1)./(b-a);
        shift=(a+b)./(a-b);
    else
        factor=0.404*ones(1,i4-i3+1);
        shift=-1.002*ones(1,i4-i3+1);
    end;
end;
Y=theta(:,i3:i4)';
for i=1:(i4-i3+1),
    Y(i,:)=factor(i)*Y(i,:)+shift(i);
end;
iter=0;
%clc;
npar=4;
y=ones(4);
xmin=0;
xmax=150e-6;
nx=50;
knn=1:nx;
kreal=1:nx;
diepte=1:nx;
inprof='knik ';
outprof='knik ';

%i=inputnnn('Press 1 to load parameterfile param.nn ','num');
i=1;
if i==1,
%     fprintf('Parameters are loaded from file param.nn ...\n');
    %fid1=fopen('param.nn');
A=param;
%     fclose(fid1);
    i=A(1);
    if i==1,
        inprof='tanh ';
    elseif i==2,
        inprof='knik ';
    end;
    i=A(2);
    if i==1,
        outprof='tanh ';
    elseif i==2,
        outprof='knik ';
    end;
    i1=round(A(3));
    i2=round(A(4));
    i3=round(A(5));
    i4=round(A(6));
    i5=round(A(7));
    i6=round(A(8));
    xmin=A(9);
    xmax=A(10);
    nx=A(11);
    maxiter=A(12);
    if scaling==0,
        itr=A(13);
        N=A(14);
    else
        itr=round(N/2);
    end;
    units=A(15);
    outunit=A(16);
    norm=A(17);
    red=A(18);
    ruis=A(19);
    clear A;
    [o,p]=size(theta);
    if N>o,
        N=o;
        itr=round(N/2);
        fprintf('matrix shorter than expected: reduced to N=%i and itr=%i\n',N,itr);
    end;
    if N<o,
        N=o;
        fprintf('matrix longer than expected: increased to N=%i \n',N);
    end;
    q=[i1,i2,i3,i4,i5,i6];
    if (max(q)>p),
        fprintf('Matrix too large, please give new sizes \n');
        fprintf('  (old sizes: i1: %i i2: %i i3: %i i4: %i i5: %i i6: %i) \n',i1,i2,i3,i4,i5,i6);
        i1=inputnnn('new i1 ? ','num');
        i2=inputnnn('new i2 ? ','num');
        i3=inputnnn('new i3 ? ','num');
        i4=inputnnn('new i4 ? ','num');
        i5=inputnnn('new i5 ? ','num');
        i6=inputnnn('new i6 ? ','num');
    end;
    clear o p q;
    PHI=theta(:,i1:i2)';
    clear m sd;
    if scaling==0,
        m=mean(PHI');% gemiddelde per signaal over profielen
        sd=sqrt(dot(PHI',PHI')/N-m.*m);% standaard deviatie per signaal over profielen
        a=min(theta(:,i3:i4));
        b=max(theta(:,i3:i4));
        if b-a~=0,
            factor=2*ones(1,i4-i3+1)./(b-a);
            shift=(a+b)./(a-b);
        else
            factor=0.404*ones(1,i4-i3+1);
            shift=-1.002*ones(1,i4-i3+1);
        end;
    else,
        load m.nn
        load sd.nn
        load factor.nn
        load shift.nn
    end;
    Y=theta(:,i3:i4)';
    for i=1:(i4-i3+1)
        Y(i,:)=factor(i)*Y(i,:)+shift(i);
    end;
end;
%clc;
fprintf('In the current configuration: \n');
fprintf('Number of hidden units : %i\n',units);
if outunit==1,
    fprintf('outunit=linear\n');
else
    fprintf('outunit=hyperbolic tangent\n');
end;
fprintf('Data matrix theta.nn contains %i rows (spectra) of \n',N);
fprintf(' input  columns %i ... %i\n',i1,i2);
fprintf(' output columns %i ... %i.\n',i3,i4);
fprintf('weight files: w1.nn w2.nn\n');
fprintf('error convergence file: pimat.nn\n');
fprintf('used configuration numbers: param.nn\n');
%fprintf('factor:  %4.3e \n',factor);
%fprintf('shift:   %4.3e \n',shift);
fprintf('input range: %i ... %i  output range: %i ... %i \n',i1,i2,i3,i4);
fprintf('# spectra: %i\n',N);
fprintf('# training spectra: %i\n',itr);
fprintf('max. # iterations: %i\n',maxiter);
fprintf('noise level (percent or absolute): %4.3e',ruis);
fprintf('\nPress a key to continue\n');

%---------- Define network structure and initialize weights ----------
rand('seed',0);
w1         = rand(units,(i2-i1+2));  % weights to hidden layer
w2         = rand(1,units+1);  % weights to output

if outunit==1,
    NetDef=['H'
        'L'];
elseif outunit==2,
    NetDef=['H'
        'H'];
end;
a=['H'
    '-'];
for i=1:units-1,
    NetDef=[NetDef,a];
end;
clear a;
rc=0;
%rc=input('Druk 1 voor rekencentrum automatische versie');

% ----------- Menu ----------
%   k = menu('Main menu nnn.m (last modification 9-7-2000):',...
% 	    'load new theta.nn or similar file',...
% 	    '***trained NN analysis***',...
% 	    'Change configuration',...
% 	    'Keyboard Control',...
%             '***Load/save data/weights/add noise to theta0.nn***',...
%             'Save scaling files',...
%             '***NN training***(based on theta.nn)',...
%             '***Apply NN***    (calculates thetam.nn)',...
%             '***A posteriori errors*** (requires thetaq.nn)',...
% 	    'quit program');
% k=5
%        l=menu('Load/save menu nnn.m',...
%             'Save theta.nn',...
% 	    'Load w1,w2 for specific parameter column',...
% 	    '***Add noise***',...
%             'Save w1,w2  to  files w1.nn and w2.nn',...
% 	    'Load w1,w2 from files w1.nn and w2.nn',...
% 	    'Create random w1,w2',...
%             'Compare noisy theta.nn with theta0.nn',...
% %             'Return to main menu');
% 
% %%%%%%%%%%%%%%%%Add noise to data
% l=3;
% if l==3,
%     theta00=theta;
%     nn=floor((i2-i1+1)/2);
%     theta0=theta;
%     [N,breedte]=size(theta0);
%     % --- Add gaussian noise to data
%     %     ll=1;
%     %     while ll~=5,
%     %     ll = menu(['Noise menu: spectra colums: i1: ',int2str(i1),' i2: ',int2str(i2),' nn: ',int2str(nn)],...
%     %             'Gaussian complex noise (% of amp)',...
%     %             'Gaussian amplitude noise (% of amp)',...
%     %             'Gaussian absolute additive noise',...
%     %             ['noise level : ',num2str(ruis)],...
%     %             'return');
%     
%     %ruis=inputnnn('noise level (% or absolute) ? ','num');
%     %%%%%%%%%Normalize theta0 before adding noise theta0(i,:)=theta0(i,:)/theta0(i,end);
%     %     load theta0.nn
%     %     plot(theta0(1:5,1:51)')
%     [r,c]=size(theta0);
% %     for i=1:r
% %         theta0(i,1:51)=theta0(i,1:51)/theta0(i,51);
% %         theta0(i,1:51)=theta0(i,1:51)/theta0(i,51);
% %     end
%     %     plot(theta0(1:5,1:51)');
%     clear r c i
%     [a,b]=size(theta);
%     fprintf('Size calculated theta: %i arrays and %i columns\n',a,b);
%     clear amp faze param;
%     %end;%end while ll~=
%     
%     clear o p q;
%     PHI=theta(:,i1:i2)';
%     clear m sd;
%     m=mean(PHI');% gemiddelde per signaal over profielen
%     sd=sqrt(dot(PHI',PHI')/N-m.*m);% standaard deviatie per signaal over profielen
%     a=min(theta(:,i3:i4));
%     b=max(theta(:,i3:i4));
%     if b-a~=0,
%         factor=2*ones(1,i4-i3+1)./(b-a);
%         shift=(a+b)./(a-b);
%     else
%         factor=0.404*ones(1,i4-i3+1);
%         shift=-1.002*ones(1,i4-i3+1);
%     end;
%     Y=theta(:,i3:i4)';
%     for i=1:(i4-i3+1)
%         Y(i,:)=factor(i)*Y(i,:)+shift(i);
%     end;
% elseif l==4,
%     % --- Save gewichten w1 en w2 op file
%     fprintf('w1 and w2 are saved to w1.nn and w2.nn ...\n');
%     save w1.nn w1 -ascii -tabs;
%     save w2.nn w2 -ascii -tabs;
%     fprintf('\nPress a key to continue\n');
%     pause;
%     
% elseif l==5,
%     % --- Load gewichten van file w1.mat en w2.mat
%     fprintf('w1 and w2 are loaded from w1.nn and w2.nn ...\n');
%     load w1.nn
%     load w2.nn
%     fprintf('\nPress a key to continue\n');
%     pause;
%     
% elseif l==6,
%     % --- Definieer random gewichten w1 en w2
%     fprintf('w1 and w2 are randomized ...\n');
%     w1         = rand(units,(i2-i1+2));  % weights to hidden layer
%     w2         = rand(1,units+1);  % weights to output
%     fprintf('\nPress a key to continue\n');
%     pause;
% end

%%%%%%%%%%%%%%%%NN training
k=7
if k==7
    fprintf('start NN training')
    % ----- Marquardt algorithm -----
    %clc;
    fprintf('In the current configuration: \n');
    fprintf('Number of hidden units : %i\n',units);
    if outunit==1,
        fprintf('outunit=linear\n');
    else
        fprintf('outunit=hyperbolic tangent\n');
    end;
    fprintf('Data matrix theta.nn contains %i rows (spectra) of \n',N);
    fprintf(' input  columns %i ... %i\n',i1,i2);
    fprintf(' output columns %i ... %i.\n',i3,i4);
    fprintf('# training spectra: %i\n',itr);
    fprintf('maximum # iterations: %i\n',maxiter);
    fprintf('Parameters are saved in file param.nn ...\n');
    fid1=fopen('param.nn','w');
    if inprof=='tanh ',
        i=1;
    elseif inprof=='knik ',
        i=2;
    end;
    if outprof=='tanh ',
        j=1;
    elseif outprof=='knik ',
        j=2;
    end;
    a=0;
    A=[i,j,i1,i2,i3,i4,i5,i6,xmin,xmax,nx,maxiter,itr,N,units,outunit,norm,red,ruis,a,a,a,a,a,a];
    fprintf(fid1,'%4.3e\n',A);
    status=fclose(fid1);
    cc1=fix(clock);
    D=0;
    fprintf('Enter start and stopparameter\n');
    %istart=inputnnn('Start parameter # ? ','num');
    %istop=inputnnn('Stop parameter # ? ','num');
    istart=i3;
    istop=i4;
    fprintf(['Start parameter ',num2str(istart),'\n']);
    fprintf(['Stop parameter ',num2str(istop),'\n']);
    if istop<istart,
        istap=-1;
    else
        istap=1;
    end;
    fprintf('nn trainings are performed from parameter %i to %i\n',i4,i3);
    rr=1;
    fprintf('randomize w1, w2 weigth matrices before each iteration');
    rand('seed',0);
    w1         = rand(units,(i2-i1+2));  % weights to hidden layer
    w2         = rand(1,units+1);  % weights to output
    for par=istart:istap:istop,
        fprintf('Training starts for parameter %i of %i:%i.\n',par,istart,istop);
        clear a b factor shift Y;
        a=min(theta(:,par)');
        b=max(theta(:,par)');
        if scaling==0,
            factor=2/(b-a);
            shift=(a+b)/(a-b);
        else
            load factor.nn
            load shift.nn
            factor=factor(par-i3+1);
            shift=shift(par-i3+1);
        end;
        Y=factor*theta(:,par)'+shift;
        if norm==0,
            PHI1=PHI(:,1:itr);
        else
            a=0*(1:itr)+1;
            PHI1=(PHI(:,1:itr)-m'*a)./(sd'*a);
            if red==1,
                PHI1=PHI1/(i2-i1+2);
            end;
            clear a;
        end;
        if norm==0,
            PHItest=PHI(:,itr+1:N);
        else
            a=0*(1:N-itr)+1;
            PHItest=(PHI(:,itr+1:N)-m'*a)./(sd'*a);
            if red==1,
                PHItest=PHItest/(i2-i1+2);
            end;
            clear a;
        end;
        Y1=Y(:,1:itr);
        Ytest=Y(:,itr+1:N);
        if rr==1,
            rand('seed',0);
            w1         = rand(units,(i2-i1+2));  % weights to hidden layer
            w2         = rand(1,units+1);  % weights to output
        end;
        if par==istart,
            %ss=inputnnn('Press 1 to load w1, w2 weigth matrices for first iteration','num');
            ss=0;
            if ss==1,
                load w1.nn;
                load w2.nn;
            end;
        end;
        lambda=1;
        trparms=[maxiter stop_crit lambda D];
        [w1,w2,PI_vector,iter,lambda]=nnmarq(NetDef,w1,w2,PHI1,Y1,PHItest,Ytest,trparms,par);
        w1naam=sprintf('w1%i.nn',par);
        fid=fopen(w1naam,'w');
        for i=1:units,
            fprintf(fid,'%1.4f   ',w1(i,:));
            fprintf(fid,'\n');
        end;
        fclose(fid);
        w2naam=sprintf('w2%i.nn',par)
        fid=fopen(w2naam,'w');
        fprintf(fid,'%1.4f   ',w2(1,:));
        fclose(fid);
    end;
    %
    cc2=fix(clock);
    fid=fopen('conf.nn','w')
    fprintf(fid,'During the performed nn training: \n');
    fprintf(fid,'Data matrix theta.nn contained %i rows (spectra) of \n',N);
    fprintf(fid,' input  columns %i ... %i\n',i1,i2);
    fprintf(fid,' output columns %i ... %i.\n',i3,i4);
    fprintf(fid,'weight files: w1.nn w2.nn\n');
    fprintf(fid,'error convergence file: pimat.nn\n');
    fprintf(fid,'used configuration numbers: param.nn\n');
    %fprintf(fid,'factor:  %4.3e  shift: %4.3e \n',factor,shift);
    fprintf(fid,'# training spectra: %i\n',itr);
    fprintf(fid,'max. # iterations: %i\n',maxiter);
    fprintf(fid,'eff. # iterations: %i\n',iter);
    fprintf(fid,'Network training started at %2i.%2i.%2i\n',cc1(4),cc1(5),cc1(6));
    fprintf(fid,'Network training ended at %2i.%2i.%2i\n',cc2(4),cc2(5),cc2(6));
    fprintf(fid,'Date: %2i.%2i.%4i\n',cc2(3),cc2(2),cc2(1));
    fprintf('During the performed nn training: \n');
    fprintf('Data matrix theta.nn contained %i rows (spectra) of \n',N);
    fprintf(' input  columns %i ... %i\n',i1,i2);
    fprintf(' output columns %i ... %i.\n',i3,i4);
    fprintf('weight files: w1.nn w2.nn\n');
    fprintf('error convergence file: pimat.nn\n');
    fprintf('used configuration numbers: param.nn\n');
    %fprintf('factor:  %4.3e  shift: %4.3e \n',factor,shift);
    fprintf('# training spectra: %i\n',itr);
    fprintf('max. # iterations: %i\n',maxiter);
    fprintf('eff. # iterations: %i\n',iter);
    fprintf('noise level (percent or absolute): %4.3e\n',ruis);
    fprintf('Network training started at %2i.%2i.%2i\n',cc1(4),cc1(5),cc1(6));
    fprintf('Network training ended at %2i.%2i.%2i\n',cc2(4),cc2(5),cc2(6));
    fprintf('Date: %2i.%2i.%4i\n',cc2(3),cc2(2),cc2(1));
    status=fclose(fid);
    if rc==1,
        pause;
    end;
end

%%%%%%%%%%%%%%%%Saving scaling files
k=6;
if k==6
    fprintf('Saving scaling files m.nn sd.nn factor.nn shift.nn \n');
    PHI=theta(:,i1:i2)';
    clear m sd;
    m=mean(PHI');% gemiddelde per signaal over profielen
    sd=sqrt(dot(PHI',PHI')/N-m.*m);% standaard deviatie per signaal over profielen
    a=min(theta(:,i3:i4));
    b=max(theta(:,i3:i4));
    if b-a~=0,
        factor=2*ones(1,i4-i3+1)./(b-a);
        shift=(a+b)./(a-b);
    else
        factor=0.404*ones(1,i4-i3+1);
        shift=-1.002*ones(1,i4-i3+1);
    end;
    save m.nn m -ascii -tabs
    save sd.nn sd -ascii -tabs
    save factor.nn factor -ascii -tabs
    save shift.nn shift -ascii -tabs
end;

%%%%%%%%%applying trained NN
k=8
if k==8,
    % ----- Calculate all nn outputs -----
    %clc;
    fprintf('In the current configuration: \n');
    fprintf('Number of hidden units : %i\n',units);
    if outunit==1,
        fprintf('outunit=linear\n');
    else
        fprintf('outunit=hyperbolic tangent\n');
    end;
    fprintf('Data matrix theta.nn contains %i rows (spectra) of \n',N);
    fprintf(' input  columns %i ... %i\n',i1,i2);
    fprintf(' output columns %i ... %i.\n',i3,i4);
    for i=1:(i4-i3+1),
        %fprintf('factor(%i) : %4.3e\n',i,factor(i));
        %fprintf('shift(%i) : %4.3e\n',i,shift(i));
    end;
    fprintf('# parameters: %i\n',i4-i3+1);
    fprintf('profile range: %3.4e...%3.4e   # steps: %i \n',xmin,xmax,nx);
    fprintf('Actual type of true profile:  \n');
    if inprof=='tanh ',
        fprintf('Hyperbolic tangent (1)\n');
    elseif inprof=='knik ',
        fprintf('Sloped staircase (2)\n');
    end;
    fprintf('Actual type of nn output profile: \n');
    if outprof=='tanh ',
        fprintf('Hyperbolic tangent (1)\n');
    elseif outprof=='knik ',
        fprintf('Sloped staircase (2)\n');
    end;
    %fprintf('Press 1 to calculate and save new thetam.nn\n');
    %letter=inputnnn(' (old spectra with new nn profile) ','num');
    letter=1;
    if letter==1,
        fprintf('calculate and save new thetam.nn\n');
        % bereken  nn parameters voor spectra van PHI
        fprintf('thetam.nn is calculated...\n');
        tic
        thetam=theta;
        spectrum=PHI(:,:);
        [rij,kolom]=size(PHI);
        s=ones(kolom,1);
        if norm==0,
            spectrum1=spectrum;
        else
            spectrum1=(spectrum-(s*m)')./(s*sd)';
        end;
        if red==1,
            spectrum1=spectrum1/(i2-i1+2);
        end;
        for i=i5:i6,
            thetam(:,[i])=(nncalc(spectrum1,outprof,i,units,outunit,i2-i1+1))';
            thetam(:,[i])=(thetam(:,[i])-shift(:,i-i5+1)*s)./(factor(:,i-i5+1)*s);	%gecorrigeerd 6 juni 2004
        end;
        toc
        P=thetam;
        fprintf('thetam.nn is saved...\n\n');
        save thetam.nn thetam -ascii -tabs;

    end;
end

%%%%%%%%%%%%%%%%analyse results
k=2;
if k==2
    % --- Berekenen van nn outputs en figuur maken
    % --- kiezen tussen figuren
    h=figure;
    set(h,'PaperType','a4letter');
    % get(h,'Position');
    % A=get(h,'Position');
    % set(h,'Position',[round(2.5*A(1)),round(2.5*A(2)),round(0.7*A(3)),round(0.7*A(4))]);
    s1=['N=',int2str(N)];
    s2=['i1..i2:',int2str(i1),'..',int2str(i2)];
    s3=['i3..i4:',int2str(i3),'..',int2str(i4)];
    s4=['fac(1): ',num2str(factor(1)),' shift(1):',num2str(shift(1))];
    s5=['# train.spectra: ',int2str(itr)];
    s6=['max # iter.:',int2str(maxiter),' norm(1/0): ',int2str(norm),' red(1/0): ',int2str(red)];
    s7=['#units: ',int2str(units)];
    s8=['type output unit(1=L,2=tan): ',int2str(outunit)];
    s9=['time: ',int2str(cc2(4)),'.',int2str(cc2(5)),'.',int2str(cc2(6))];
    s10=['date: ',int2str(cc2(3)),'.',int2str(cc2(2)),'.',int2str(cc2(1))];
    
    %     while j~=6,
    %     j = menu('View menu nnn.m:',...
    % 	    'View network architecture',...
    % 	    'View iteration convergence pimat.nn',...
    % 	    'View histogram test and training errors',...
    % 	    'Plot correlation between parameters and nn error',...
    %             'Calculate characteristic performances',...
    % 	    'Return to main menu');
    
    %%%%%%%%%%%%%%%%%%%%Calculate characteristic performances
    
    if j==1,
        % --- figuur van network architecture
        subplot(1,1,1);
        drawnet(w1,w2,1e-9);
        title('Network Architecture');
        drawnow;
        %pause;
    end
    
    if j==2
        load pimat.nn
        subplot(1,1,1);
        plot(pimat(:,1),pimat(:,3),'+');
        title('Mean Training Error vs Training Time');
        drawnow;
        %pause;
        clear pimat
    end
    
    if j==3,
        % --- plot van histogram van Y1-y21 en Y2-y22
        ianalyse=inputnnn('Give number of parameter column to analyze ','num');
        % --- Load gewichten van specifieke parameter
        fprintf('w1 and w2 are loaded from files w1%i.nn w2%i.nn\n',ianalyse,ianalyse);
        a=min(theta(:,ianalyse)');
        b=max(theta(:,ianalyse)');
        if scaling==0,
            factor=2/(b-a);
            shift=(a+b)/(a-b);
        else
            load factor.nn
            load shift.nn
            factor=factor(ianalyse-i3+1);
            shift=shift(ianalyse-i3+1);
        end;
        Y=factor*theta(:,ianalyse)'+shift;
        w1naam=sprintf('w1%i.nn',ianalyse);
        w2naam=sprintf('w2%i.nn',ianalyse);
        %w1T en w2T worden opgevuld per kolom tijdens fscanf
        fid=fopen(w1naam,'r');
        [w1T,count]=fscanf(fid,'%e',[i2-i1+2,units]);
        fclose(fid);
        w1=w1T';
        fid=fopen(w2naam,'r');
        [w2T,count]=fscanf(fid,'%e',[units+1,1]);
        w2=w2T';
        fclose(fid);
        clear w1T w2T;
        % --- vergelijking trainingsdata Y1 van theta.nn met voorspelling
        if norm==0,
            PHI0=PHI(:,1:itr);
        else
            a=0*(1:itr)+1;
            PHI0=(PHI(:,1:itr)-m'*a)./(sd'*a);
            clear a;
        end;
        if red==1,
            PHI0=PHI0/(i2-i1+2);
        end;
        Y1=Y(:,1:itr);
        [outputs,itr] = size(Y1);                 % # of outputs and # of data
        y21=nnetout(NetDef,w1,w2,PHI0);
        E        = (Y1- y21);                      % Training error
        PI1      = (sum(sum(E.*E))/itr)^(0.5);       % Sum of squared errors
        clear PHI0;
        % --- vergelijking testdata Y2 van theta.nn met voorspelling
        [outputs,N] = size(Y);
        if norm==0,
            PHI0=PHI(:,itr+1:N);
        else
            a=0*(1:N-itr)+1;
            PHI0=(PHI(:,itr+1:N)-m'*a)./(sd'*a);
            clear a;
        end;
        if red==1,
            PHI0=PHI0/(i2-i1+2);
        end;
        Y2=Y(:,itr+1:N);
        [outputs,itest] = size(Y2);                 % # of outputs and # of data
        y22=nnetout(NetDef,w1,w2,PHI0);
        E        = (Y2- y22);                      % Test error
        PI2      = (sum(sum(E.*E))/itest)^(0.5);       % Sum of squared errors
        clear PHI0;
        subplot(2,2,1);
        [a,b]=size(Y1);
        b=round(b/20);
        if b<10,
            b=10;
        end;
        hist(Y1-y21,b);
        s=['Train.Err.(kol ',int2str(ianalyse),'):',num2str(PI1)];
        title(s);
        subplot(2,2,2);
        [a,b]=size(Y2);
        b=round(b/20);
        if b<10,
            b=10;
        end;
        hist(Y2-y22,b);
        s=['Test.Err.(kol ',int2str(ianalyse),'):',num2str(PI2)];
        title(s);
        subplot(2,2,3);
        t1=text(0,1.0,s1);
        t2=text(0,0.9,s2);
        t3=text(0,0.8,s3);
        t4=text(0,0.7,s4);
        t5=text(0,0.6,s5);
        t6=text(0,0.5,s6);
        t7=text(0,0.4,s7);
        t8=text(0,0.3,s8);
        t9=text(0,0.2,s9);
        t10=text(0,0.1,s10);
        subplot(2,2,4);
        drawnet(w1,w2,1e-9);
        title('Network Architecture');
        drawnow;
        %pause;
    end
    
    if j==4,
        % --- plot van correlaties tussen nn fout en parameters
        fprintf('\nCorrelation plot nn error vs parameter \n');
        ianalyse=inputnnn('Give number of error column ','num');
        % --- Load gewichten van specifieke parameter
        fprintf('w1 and w2 are loaded from files w1%i.nn w2%i.nn\n',ianalyse,ianalyse);
        a=min(theta(:,ianalyse)');
        b=max(theta(:,ianalyse)');
        if scaling==0,
            factor=2/(b-a);
            shift=(a+b)/(a-b);
        else
            load factor.nn
            load shift.nn
            factor=factor(ianalyse-i3+1);
            shift=shift(ianalyse-i3+1);
        end;
        Y=factor*theta(:,ianalyse)'+shift;
        w1naam=sprintf('w1%i.nn',ianalyse);
        w2naam=sprintf('w2%i.nn',ianalyse);
        %w1T en w2T worden opgevuld per kolom tijdens fscanf
        fid=fopen(w1naam,'r');
        [w1T,count]=fscanf(fid,'%e',[i2-i1+2,units]);
        fclose(fid);
        w1=w1T';
        fid=fopen(w2naam,'r');
        [w2T,count]=fscanf(fid,'%e',[units+1,1]);
        w2=w2T';
        fclose(fid);
        clear w1T w2T;
        % --- vergelijking trainingsdata Y1 van theta.nn met voorspelling
        if norm==0,
            PHI0=PHI(:,1:itr);
        else
            a=0*(1:itr)+1;
            PHI0=(PHI(:,1:itr)-m'*a)./(sd'*a);
            clear a;
        end;
        if red==1,
            PHI0=PHI0/(i2-i1+2);
        end;
        Y1=Y(:,1:itr);
        [outputs,itr] = size(Y1);                 % # of outputs and # of data
        y21=nnetout(NetDef,w1,w2,PHI0);
        E        = (Y1- y21);                      % Training error
        PI1      = (sum(sum(E.*E))/itr)^(0.5);       % Sum of squared errors
        clear PHI0;
        % --- vergelijking testdata Y2 van theta.nn met voorspelling
        [outputs,N] = size(Y);
        if norm==0,
            PHI0=PHI(:,itr+1:N);
        else
            a=0*(1:N-itr)+1;
            PHI0=(PHI(:,itr+1:N)-m'*a)./(sd'*a);
            clear a;
        end;
        if red==1,
            PHI0=PHI0/(i2-i1+2);
        end;
        Y2=Y(:,itr+1:N);
        [outputs,itest] = size(Y2);                 % # of outputs and # of data
        y22=nnetout(NetDef,w1,w2,PHI0);
        E        = (Y2- y22);                      % Test error
        PI2      = (sum(sum(E.*E))/itest)^(0.5);       % Sum of squared errors
        clear PHI0;
        
        ipar=inputnnn('Give number of parameter ','num');
        [N,kolommen]=size(theta);
        parc=theta(:,ipar);
        y01=theta(1:itr,ipar)';
        y11=(y21-Y1)/factor;
        y02=theta(itr+1:N,ipar)';
        y12=(y22-Y2)/factor;
        yte=[y11';y12'];
        parcc=[parc,yte];
        subplot(2,2,1);
        plot(y01,y11,'+');
        s=['abs nn training error ',int2str(ianalyse),' vs abs parameter ',int2str(ipar)];
        title(s);
        subplot(2,2,2);
        plot(y02,y12,'+');
        s=['abs nn test error ',int2str(ianalyse),' vs abs parameter ',int2str(ipar)];
        title(s);
        subplot(2,2,3);
        t1=text(0,1.0,s1);
        t2=text(0,0.9,s2);
        t3=text(0,0.8,s3);
        t4=text(0,0.7,s4);
        t5=text(0,0.6,s5);
        t6=text(0,0.5,s6);
        t7=text(0,0.4,s7);
        t8=text(0,0.3,s8);
        t9=text(0,0.2,s9);
        t10=text(0,0.1,s10);
        subplot(2,2,4);
        drawnet(w1,w2,1e-9);
        title('Network Architecture');
        drawnow;
        fprintf('Press a key to continue \n');
        pause;
        i=inputnnn('Press 1 to save result to file corr.nn ','num');
        if i==1,
            save corr.nn parcc -ascii -tabs;
        end;
        %fprintf('Press a key to continue \n');
        %pause;
        clear y0 y1;
    end
    
    j=5;
    if j==5,
        % plot van nn error voor verschillende parameters
        fprintf('Calculation of characteristic performances \n');
        [x,y]=size(theta);
        fprintf('Size theta: %i x %i\n',x,y);
        fprintf('Data range: %i ... %i\n',i1,i2);
        clear x y;
        fprintf('Give range i5..i6 of parameter columns to analyze\n');
        %     i5=inputnnn('value for i5? ','num');
        %     i6=inputnnn('value for i6? ','num');
        %     uu=inputnnn('Press 1 to save all true-NN parameter pairs on files train*.nn test*.nn ','s');
        %     vv=inputnnn('Press 1 to view all performances ','s');
        i5=i3;
        i6=i4;
        uu='1';
        vv='1';
        %fprintf(['value for i5  ',num2str(i5),'\n'])
        %fprintf(['value for i6  ',num2str(i6),'\n'])
        fprintf('save all true-NN parameter pairs on files train*.nn test*.nn \n')
        fprintf('View all performances');
        %vv=inputnnn('View all performances ');
        if isempty(uu),
            uu='0';
        end;
        if isempty(vv),
            vv='0';
        end;
        trainf=1:(i6-i5+1);
        testf=1:(i6-i5+1);
        for nr=1:(i6-i5+1),
            clear w1 w2;
            i3a=i3;
            i4a=i4;
            i3=i5+nr-1;
            i4=i3;
            fprintf('Errors for parameter column %i : \n',i3);
            %fprintf('    (Loading w1, w2 ...)\n');
            a=min(theta(:,i3)');
            b=max(theta(:,i3)');
            if scaling==0,
                factor=2/(b-a);
                shift=(a+b)/(a-b);
            else,
                load factor.nn
                load shift.nn
                factor=factor(i3-i3a+1);
                shift=shift(i3-i3a+1);
            end;
            Y=factor*theta(:,i3)'+shift;
            % inladen w1 w2 voor aktuele parameter
            w1naam=sprintf('w1%i.nn',i3);
            w2naam=sprintf('w2%i.nn',i3);
            %w1T en w2T worden opgevuld per kolom tijdens fscanf
            fid=fopen(w1naam,'r');
            [w1T,count]=fscanf(fid,'%e',[(i2-i1+2),units]);
            fclose(fid);
            w1=w1T';
            fid=fopen(w2naam,'r');
            [w2T,count]=fscanf(fid,'%e',[units+1,1]);
            w2=w2T';
            fclose(fid);
            clear w1T w2T;
            % berekenen van trainingsfouten voor aktuele parameter
            if norm==0,
                PHI0=PHI(:,1:itr);
            else
                a=0*(1:itr)+1;
                PHI0=(PHI(:,1:itr)-m'*a)./(sd'*a);
                clear a;
            end;
            if red==1,
                PHI0=PHI0/(i2-i1+2);
            end;
            Y1=Y(:,1:itr);
            [outputs,N] = size(Y1);                 % # of outputs and # of data
            y21=nnetout(NetDef,w1,w2,PHI0);
            E        = (Y1- y21);                      % Training error
            trainf(nr)  = (sum(sum(E.*E))/N)^(0.5);       % Sum of squared errors
            fprintf('Training error: %4.3e \n',trainf(nr));
            clear PHI0;
            % --- vergelijking testdata Y2 van theta.nn met voorspelling
            [outputs,N] = size(Y);
            if norm==0,
                PHI0=PHI(:,itr+1:N);
            else
                a=0*(1:N-itr)+1;
                PHI0=(PHI(:,itr+1:N)-m'*a)./(sd'*a);
                clear a;
            end;
            if red==1,
                PHI0=PHI0/(i2-i1+2);
            end;
            Y2=Y(:,itr+1:N);
            [outputs,N] = size(Y2);                 % # of outputs and # of data
            y22=nnetout(NetDef,w1,w2,PHI0);
            E        = (Y2- y22);                      % Test error
            testf(nr)  = (sum(sum(E.*E))/N)^(0.5);       % Sum of squared errors
            fprintf('Test error    : %4.3e \n',testf(nr));
            clear PHI0;
            Y1=(Y1-shift)./factor;
            y21=(y21-shift)./factor;
            Y2=(Y2-shift)./factor;
            y22=(y22-shift)./factor;
            if vv=='1',
                subplot(2,2,1);
                plot(Y1,y21,'+',Y1,Y1,'-');
                s=['Abs. training error (',int2str(nr),'):',num2str(trainf(nr))];
                title(s);
                subplot(2,2,2);
                plot(Y2,y22,'+',Y2,Y2,'-');
                s=['Abs. test error (',int2str(nr),'):',num2str(testf(nr))];
                title(s);
                subplot(2,2,3);
                t1=text(0,1.0,s1);
                t2=text(0,0.9,s2);
                t3=text(0,0.8,s3);
                t4=text(0,0.7,s4);
                t5=text(0,0.6,s5);
                t6=text(0,0.5,s6);
                t7=text(0,0.4,s7);
                t8=text(0,0.3,s8);
                t9=text(0,0.2,s9);
                t10=text(0,0.1,s10);
                subplot(2,2,4);
                drawnet(w1,w2,1e-9);
                title('Network Architecture');
                drawnow;
                %fprintf('Press a key to continue \n');
                %pause;
            end;
            if uu=='1',
                fprintf('saving train*.nn \n')
                trnaam=sprintf('train%i.nn',i3);
                TR=[Y1',y21']';
                fid = fopen(trnaam,'w');
                fprintf(fid,'%4.3e  %4.3e\n',TR);
                status=fclose(fid);
                fprintf('saving test*.nn \n')
                tenaam=sprintf('test%i.nn',i3);
                fid = fopen(tenaam,'w');
                TE=[Y2',y22']';
                fprintf(fid,'%4.3e  %4.3e\n',TE);
                status=fclose(fid);
                clear Y1 Y1 y21 y22 TR TE PHI0;
            end;
            i3=i3a;
            i4=i4a;
        end;
        %ii=inputnnn('Press 1 to view Training/Test error vs parameter ','num');
        ii=1;
        if ii==1,
            fprintf('View Training/Test error vs parameter\n')
            % --- plot van fout ifv parameter
            close(h);
            h=figure;
            subplot(2,2,1);
            plot(i5:i6,trainf,'+');
            s=['Rel.Train.Err. vs parameter:'];
            title(s);
            subplot(2,2,2);
            plot(i5:i6,testf,'+');
            s=['Rel.Test.Err. vs parameter:'];
            title(s);
            subplot(2,2,3);
            drawnet(w1,w2,1e-9);
            title('Network Architecture');
            subplot(2,2,4);
            t1=text(0,1.0,s1);
            t2=text(0,0.9,s2);
            t3=text(0,0.8,s3);
            t4=text(0,0.7,s4);
            t5=text(0,0.6,s5);
            t6=text(0,0.5,s6);
            t7=text(0,0.4,s7);
            t8=text(0,0.3,s8);
            t9=text(0,0.2,s9);
            t10=text(0,0.1,s10);
            drawnow;
        end;
        %i=inputnnn('Press 1 to save errors to files trainf.nn testf.nn totalf.nn ','num');
        i=1;
        if i==1,
            fprintf('save errors to files trainf.nn testf.nn totalf.nn \n')
            fid1=fopen('trainf.nn','w');
            fprintf(fid1,'%4.3e\n',trainf);
            status=fclose(fid1);
            fid1=fopen('testf.nn','w');
            fprintf(fid1,'%4.3e\n',testf);
            status=fclose(fid1);
        end;
        %fprintf('Press a key to continue\n');
        %pause;
    end
    
end
%   %close(h);
%   clear PHI0 y1 y2 h1 h2 y21 y22 y23 Y1 Y2 Y3;
%   clear s1 s2 s3 s4 s5 s6 s7 s8 s9 s10
%   clear t1 t2 t3 t4 t5 t6 t7 t8 t9 t10


%%%%%%%%%%%%%%%generate report and save figs
