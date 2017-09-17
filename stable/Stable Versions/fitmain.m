fprintf('program fitmain.m \n');
fprintf('correction of error bar calculation: July 24, 2002\n');
fprintf('last change: October 30, 2009\n');
global Texp Sexp counter Sstart h chi2vec n0 Ifitpar Ifixpar i1 i2 chi2vecext datatype xscalestring yscalestring functionstring ifitstring
global Sact parstr1 parstr2 parstr3 parstr4 parstr5 parstr6 parstr7 parstr8 parstr9
global Ifitparnew Ifixparnew
global iplotfit nplot
global chihistory chilimhistory

if ~exist('iexecuted'),
nplot=5;
iplotfit=1;
n0=9;			%total number of parameters
I0=1:n0;		%parameter index array
Ifitpar=1:4;		%fit parameter index array
a=ones(1,n0);
a(Ifitpar)=0;
Ifixpar=find(a);	%fix parameter index array
npar=length(Ifitpar);
nfix=n0-npar;
startparam=zeros(1,n0);	%values of parameters
startparam(Ifitpar)=1;	%values of fit parameters
fitresult=startparam(Ifitpar); %fitresult
name=[];
chi2vec=[];
chi2vecext=[];
datatype='real    ';%datatype='realimag';datatype='complex ';

ifitstring=0;

ipar=1;			%fit parameter index for error test
isqrt=0;
parmin=0;		%fit parameter range for error test
parmax=1;
ipar1=1;		%fit parameter index for contour test
ipar2=2;		%fit parameter index for contour test
par1min=0;		%fit parameter range for contour test
par1max=1;
par2min=0;		%fit parameter range for contour test
par2max=1;
ntest=20;		%number of error analysis values
ndata=100;		%number of data
Texp=(1:ndata)';	%T-data
Sexp=(1:ndata)';	%S-data
i1=1;			%start fit range
i2=ndata;		%end fit range
Iact=[i1:i2];
iexecuted=1;
xscale=0;xscalestring='lin';
yscale=0;yscalestring='lin';
functionstring='funmain';
parstr1='';
parstr2='';
parstr3='';
parstr4='';
parstr5='';
parstr6='';
parstr7='';
parstr8='';
parstr9='';

end;


m=0;
while m~=21,
 m=menu(['fitmain.m 3-5-2002- datatype: ',datatype],...
        ['starting values - par1:',parstr1,': ',num2str(startparam(1)),' par2: ',parstr2,': ',num2str(startparam(2)),' par3: ',parstr3,': ',num2str(startparam(3))],...
        ['starting values - par4:',parstr4,': ',num2str(startparam(4)),' par5: ',parstr5,': ',num2str(startparam(5)),' par6: ',parstr6,': ',num2str(startparam(6))],...
        ['starting values - par7:',parstr7,': ',num2str(startparam(7)),' par8: ',parstr8,': ',num2str(startparam(8)),' par9: ',parstr9,': ',num2str(startparam(9))],...
        ['fit parameter index array: ',int2str(Ifitpar),'- fix parameter index array: ',int2str(Ifixpar)],...
        ['starting fit parameters: ',num2str(startparam(Ifitpar)),' - press here to set fit parameters to start parameters'],...
        ['best fit parameters    : ',num2str(fitresult),' - press here to accept for next fit'],...
        'load/save',...
        ['LOAD DATA (2 columns (real) or 3 columns (complex amp phase)) ',name,' -  number of data: ',int2str(ndata)],...
        ['choose data range for fitting: from ',int2str(i1),' to ',int2str(i2)],...
        'START FIT',...
        'CALCULATE FOR CURRENT PARAMETERS',...
        'CHI2 ERROR ANALYSIS',...
        'CHI2 CONTOUR PLOT',...
        ['xscale: (0=lin,1=log): ',int2str(xscale)],...
        ['yscale: (0=lin,1=log): ',int2str(yscale)],...
        ['fitfunction name: ',functionstring],...
        ['plot during fit every ',int2str(nplot),' iterations'],...
        'CHI2 ERROR ANALYSIS WITH MINIMIZATION OVER OTHER PARAMETERS FOR EVERY PARAMETER VALUE',...
        ['sqrt x-scale: (0/1): ',int2str(isqrt)],...
        ['fittype (0:complex,1:amp,2:faze,3:real,4:imag,5:log(amp)): ',int2str(ifitstring)],...
        'end');

if m==1,
%initial parameters
startparam(1)=input([parstr1,': starting parameter 1? ']);
startparam(2)=input([parstr2,': starting parameter 2? ']);
startparam(3)=input([parstr3,': starting parameter 3? ']);
elseif m==2,
startparam(4)=input([parstr4,': starting parameter 4? ']);
startparam(5)=input([parstr5,': starting parameter 5? ']);
startparam(6)=input([parstr6,': starting parameter 6? ']);
elseif m==3,
startparam(7)=input([parstr7,': starting parameter 7? ']);
startparam(8)=input([parstr8,': starting parameter 8? ']);
startparam(9)=input([parstr9,': starting parameter 9? ']);
elseif m==4,
Ifitpar=input('fit parameter index array (e.g. [1 5 8])? ');
npar=length(Ifitpar);
a=ones(1,n0);
a(Ifitpar)=0;
Ifixpar=find(a);	%fix parameter index array
nfix=length(Ifixpar);
fitresult=startparam(Ifitpar);
elseif m==5,
fprintf('starting fit parameters: \n');
startparam(Ifitpar)
fitresult=startparam(Ifitpar);
elseif m==6,
fprintf('best fit parameters: \n');
fitresult
if size(Ifitpar)==size(fitresult),
startparam(Ifitpar)=fitresult; 
else,
fprintf('not accepted: sizes incompatible \n');
end;

elseif m==7,
m0=0;
while m0~=6,
 m0=menu('Load/save menu',...
        'save signal for best fit',...
        'save signal for starting parameters',...
        'save settings to fitmain.par',...
        'load settings from fitmain.par',...
        'keyboard',...
        'return');
 if m0==1,
   param=startparam;
   param(Ifitpar)=fitresult;
   Sfit=eval([functionstring,'(Texp,param)']);
   [name,pnaam]=uiputfile('*.*',['choose file to save signal for best fit']);
   if datatype=='real    ',
    result=[Texp,Sfit];
    elseif datatype=='complex ',
    result=[Texp,abs(Sfit),angle(Sfit)*180/pi];
    elseif datatype=='realimag',
    result=[Texp,real(Sfit),imag(Sfit)];
    
    end;%end if datatype==
   savemat([pnaam,name],result)
  elseif m0==2,
   Sstart=eval([functionstring,'(Texp,startparam)']);
   [name,pnaam]=uiputfile('*.*',['choose file to save signal for starting parameters']);
   if datatype=='real    ',
    result=[Texp,Sstart];
    elseif datatype=='complex ',
    result=[Texp,abs(Sstart),angle(Sstart)*180/pi];
    elseif datatype=='realimag',
    result=[Texp,real(Sstart),imag(Sstart)*180/pi];
    end;%end if datatype==
   savemat([pnaam,name],result)
  elseif m0==3,
  A=[n0;npar;nfix;Ifitpar';Ifixpar';startparam';startparam(Ifitpar)';fitresult'];
  save fitmain.par A -ascii -tabs
  elseif m0==4,
  A=load('fitmain.par')';
  n0=A(1);
  npar=A(2);
  nfix=A(3);
  Ifitpar=A(4:npar+3);
  Ifixpar=A(npar+4:npar+nfix+3);
  startparam=A(npar+nfix+4:npar+nfix+3+n0);
  startparam(Ifitpar)=A(npar+nfix+n0+4:npar+nfix+n0+3+npar);
  fitresult=A(npar+nfix+n0+npar+4:end);
  elseif m0==5,
  keyboard;
  end;%end if m0==
 end;%end while m0~=


elseif m==8,
%loading exp data
[name,pnaam]=uigetfile('*.*',['choose 2-column x,y-file or 3-column x,amp,phase-file with data']);
signal=load([pnaam,name]);
[r,k]=size(signal);
if k==2,
Texp=signal(:,1);
Sexp=signal(:,2);
datatype='real    ';
elseif k==3,
keuze=input('press 0 for amplitude-phase(o) or 1 for real-imag format ');
Texp=signal(:,1);
if keuze==0,
Sexp=signal(:,2).*exp(sqrt(-1)*signal(:,3)*pi/180);
datatype='complex ';
else,
Sexp=signal(:,2)+sqrt(-1)*signal(:,3);
datatype='realimag';
end;
end;
ndata=length(Texp);
fprintf('filename: %s \n',name);
fprintf('number of data: %i \n',ndata);
i1=1;
i2=ndata;
Iact=[i1:i2];;
elseif m==9,
%choose data range
Iprev=i1:i2;
Iact=i1:i2;
if datatype=='real    ',
h1=figure;
loglog(Texp,Sexp,Texp(Iprev),Sexp(Iprev),'.',Texp(Iact),Sexp(Iact),'.');
xlabel('all data, previous fit selection, actual selection');
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);

zoom on;
as=axis;
mm=0;
while mm~=7,
figure(h1);
loglog(Texp,Sexp,Texp(Iprev),Sexp(Iprev),'.',Texp(Iact),Sexp(Iact),'.');
xlabel('all data, previous fit selection, actual selection');
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
axis(as);
zoom on;
mm=menu('data range selection',...
        'click left data limit',...
        'click right data limit',...
        'back to full range',...
        'back to previous range',...
        ['manual range input, actual range: ',int2str([Iact(1),Iact(end)])],...
        'zoom out',...
        'accept and continue');
if mm==1,
  zoom off;
  [x,y]=ginput(1);
  [a,i1act]=min(abs(Texp-x));
  Iact=i1act:max(Iact);
  zoom on;
  elseif mm==2,
  zoom off;
  [x,y]=ginput(1);
  [a,i2act]=min(abs(Texp-x));
  Iact=min(Iact):i2act;
  zoom on;
  elseif mm==3,
  Iact=1:ndata;
  elseif mm==4,
  Iact=Iprev;
  elseif mm==5,
  Iact=input('range (e.g. 3:45) ?  ');
  elseif mm==6,
  figure(h1);
  zoom out;
  end;%end if mm==
as=axis;
end;%end while mm~=

elseif datatype=='complex ',
h1=figure;
subplot(2,1,1);
loglog(Texp,abs(Sexp),'.',Texp(Iprev),abs(Sexp(Iprev)),'.',Texp(Iact),abs(Sexp(Iact)),'.');
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
xlabel('all data, previous fit selection, actual selection');
ylabel('amp');
zoom on;
as=axis;
subplot(2,1,2);
semilogx(Texp,angle(Sexp)*180/pi,'.',Texp(Iprev),angle(Sexp(Iprev))*180/pi,'.',Texp(Iact),angle(Sexp(Iact))*180/pi,'.');
set(gca,'Xscale',xscalestring);
xlabel('all data, previous fit selection, actual selection');
ylabel('phase (°)');
zoom on;

mm=0;
while mm~=7,
figure(h1);
subplot(2,1,1);
loglog(Texp,abs(Sexp),'.',Texp(Iprev),abs(Sexp(Iprev)),'.',Texp(Iact),abs(Sexp(Iact)),'.');
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
xlabel('all data, previous fit selection, actual selection');
axis(as);
zoom on;
mm=menu('data range selection',...
        'click left data limit',...
        'click right data limit',...
        'back to full range',...
        'back to previous range',...
        ['manual range input, actual range: ',int2str([Iact(1),Iact(end)])],...
        'zoom out',...
        'accept and continue');
if mm==1,
  zoom off;
  [x,y]=ginput(1);
  [a,i1act]=min(abs(Texp-x));
  Iact=i1act:max(Iact);
  zoom on;
  elseif mm==2,
  zoom off;
  [x,y]=ginput(1);
  [a,i2act]=min(abs(Texp-x));
  Iact=min(Iact):i2act;
  zoom on;
  elseif mm==3,
  Iact=1:ndata;
  elseif mm==4,
  Iact=Iprev;
  elseif mm==5,
  Iact=input('range (e.g. 3:45) ?  ');
  elseif mm==6,
  figure(h1);
  zoom out;
  end;%end if mm==
as=axis;
end;%end while mm~=

elseif datatype=='realimag',
h1=figure;
subplot(2,1,1);
loglog(Texp,real(Sexp),'.',Texp(Iprev),real(Sexp(Iprev)),'.',Texp(Iact),real(Sexp(Iact)),'.');
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
xlabel('all data, previous fit selection, actual selection');
ylabel('real');
zoom on;
as=axis;
subplot(2,1,2);
semilogx(Texp,imag(Sexp),'.',Texp(Iprev),imag(Sexp(Iprev)),'.',Texp(Iact),imag(Sexp(Iact)),'.');
set(gca,'Xscale',xscalestring);
xlabel('all data, previous fit selection, actual selection');
ylabel('imag');
zoom on;

mm=0;
while mm~=7,
figure(h1);
subplot(2,1,1);
loglog(Texp,real(Sexp),'.',Texp(Iprev),real(Sexp(Iprev)),'.',Texp(Iact),real(Sexp(Iact)),'.');
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
xlabel('all data, previous fit selection, actual selection');
axis(as);
zoom on;
mm=menu('data range selection',...
        'click left data limit',...
        'click right data limit',...
        'back to full range',...
        'back to previous range',...
        ['manual range input, actual range: ',int2str([Iact(1),Iact(end)])],...
        'zoom out',...
        'accept and continue');
if mm==1,
  zoom off;
  [x,y]=ginput(1);
  [a,i1act]=min(real(Texp-x));
  Iact=i1act:max(Iact);
  zoom on;
  elseif mm==2,
  zoom off;
  [x,y]=ginput(1);
  [a,i2act]=min(real(Texp-x));
  Iact=min(Iact):i2act;
  zoom on;
  elseif mm==3,
  Iact=1:ndata;
  elseif mm==4,
  Iact=Iprev;
  elseif mm==5,
  Iact=input('range (e.g. 3:45) ?  ');
  elseif mm==6,
  figure(h1);
  zoom out;
  end;%end if mm==
as=axis;
end;%end while mm~=


end;%end if datatype==
i1=Iact(1);
i2=Iact(end);

elseif m==10,
%fitting
Sstart=eval([functionstring,'(Texp,startparam)']);
close all;
h=figure;
counter=0;
chi2vec=[];
chi2vecext=[];

foptions=optimset;
chihistory=[];
fitresult=fminsearch('chimain',startparam(Ifitpar),foptions,startparam(Ifixpar));
save chifithistory.txt chihistory -ascii -tabs
param=startparam;
param(Ifitpar)=fitresult;
Sfit=eval([functionstring,'(Texp,param)']);
h2=figure;
if datatype=='real    ',
loglog(Texp,Sexp,'o',Texp(Iact),Sexp(Iact),'o',Texp,Sstart,Texp,Sfit);
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit');
title(['fitmain.m ',date]);
elseif datatype=='complex ',
subplot(2,1,1);
loglog(Texp,abs(Sexp),'o',Texp(Iact),abs(Sexp(Iact)),'o',Texp,abs(Sstart),Texp,abs(Sfit));
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit amp');
title(['fitmain.m ',date]);
subplot(2,1,2);
semilogx(Texp,angle(Sexp)*180/pi,'o',Texp(Iact),angle(Sexp(Iact))*180/pi,'o',Texp,angle(Sstart)*180/pi,Texp,angle(Sfit)*180/pi);
set(gca,'Xscale',xscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit phase (°)');

elseif datatype=='realimag',
subplot(2,1,1);
loglog(Texp,real(Sexp),'o',Texp(Iact),real(Sexp(Iact)),'o',Texp,real(Sstart),Texp,real(Sfit));
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit real');
title(['fitmain.m ',date]);
subplot(2,1,2);
semilogx(Texp,imag(Sexp),'o',Texp(Iact),imag(Sexp(Iact)),'o',Texp,imag(Sstart),Texp,imag(Sfit));
set(gca,'Xscale',xscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit imag');


end;%end if datatype==

fprintf('fitting parameters: \n');
fprintf('fit parameter index array: %s\n',int2str(Ifitpar));
fprintf('fit parameters: %s\n',num2str(fitresult));
fprintf('fixed parameter index array: %s\n',int2str(Ifixpar));
fprintf('fixed parameters: %s\n',num2str(startparam(Ifixpar)));


elseif m==11,
%show actual result
Sstart=eval([functionstring,'(Texp,startparam)']);
close all;
param=startparam;
param(Ifitpar)=fitresult;
Sfit=eval([functionstring,'(Texp,param)']);

figure;
if datatype=='real    ',
loglog(Texp,Sexp-Sfit,'o',Texp(Iact),Sexp(Iact)-Sfit,'o',Texp,Sstart-Sfit,'o');
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit - residues');
title(['fitmain.m ',date]);
elseif datatype=='complex ',
subplot(2,1,1);
loglog(Texp,abs(Sexp)-abs(Sfit),'o',Texp(Iact),abs(Sexp(Iact))-abs(Sfit(Iact)),'o',Texp,abs(Sstart)-abs(Sfit),'o');
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit amp - residues');
title(['fitmain.m ',date]);
subplot(2,1,2);
semilogx(Texp,angle(Sexp)*180/pi-angle(Sfit)*180/pi,'o',Texp(Iact),angle(Sexp(Iact))*180/pi-angle(Sfit(Iact))*180/pi,'o',Texp,angle(Sstart)*180/pi-angle(Sfit)*180/pi,'o');
set(gca,'Xscale',xscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit phase (°) - residues');

elseif datatype=='realimag',
subplot(2,1,1);
loglog(Texp,real(Sexp)-real(Sfit),'o',Texp(Iact),real(Sexp(Iact))-real(Sfit(Iact)),'o',Texp,real(Sstart)-real(Sfit),'o');
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit real - residues');
title(['fitmain.m ',date]);
subplot(2,1,2);
semilogx(Texp,imag(Sexp)-imag(Sfit),'o',Texp(Iact),imag(Sexp(Iact))-imag(Sfit(Iact)),'o',Texp,imag(Sstart)-imag(Sfit),'o');
set(gca,'Xscale',xscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit imag - residues');


end;%end if datatype==

figure;
if datatype=='real    ',
loglog(Texp,Sexp,'o',Texp(Iact),Sexp(Iact),'o',Texp,Sstart,Texp,Sfit);
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit');
title(['fitmain.m ',date]);
elseif datatype=='complex ',
subplot(2,1,1);
loglog(Texp,abs(Sexp),'o',Texp(Iact),abs(Sexp(Iact)),'o',Texp,abs(Sstart),Texp,abs(Sfit));
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit amp');
title(['fitmain.m ',date]);
subplot(2,1,2);
semilogx(Texp,angle(Sexp)*180/pi,'o',Texp(Iact),angle(Sexp(Iact))*180/pi,'o',Texp,angle(Sstart)*180/pi,Texp,angle(Sfit)*180/pi);
set(gca,'Xscale',xscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit phase (°)');

elseif datatype=='realimag',
subplot(2,1,1);
loglog(Texp,real(Sexp),'o',Texp(Iact),real(Sexp(Iact)),'o',Texp,real(Sstart),Texp,real(Sfit));
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit real');
title(['fitmain.m ',date]);
subplot(2,1,2);
semilogx(Texp,imag(Sexp),'o',Texp(Iact),imag(Sexp(Iact)),'o',Texp,imag(Sstart),Texp,imag(Sfit));
set(gca,'Xscale',xscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit imag');


end;%end if datatype==



if isqrt==1,

h3=figure;
if datatype=='real    ',
semilogy(sqrt(Texp),Sexp,'o',sqrt(Texp(Iact)),Sexp(Iact),'o',sqrt(Texp),Sstart,sqrt(Texp),Sfit);
%set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit');
title(['fitmain.m ',date]);
elseif datatype=='complex ',
subplot(2,1,1);
semilogy(sqrt(Texp),abs(Sexp),'o',sqrt(Texp(Iact)),abs(Sexp(Iact)),'o',sqrt(Texp),abs(Sstart),sqrt(Texp),abs(Sfit));
%set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit amp');
title(['fitmain.m ',date]);
subplot(2,1,2);
plot(sqrt(Texp),angle(Sexp)*180/pi,'o',sqrt(Texp(Iact)),angle(Sexp(Iact))*180/pi,'o',sqrt(Texp),angle(Sstart)*180/pi,sqrt(Texp),angle(Sfit)*180/pi);
%set(gca,'Xscale',xscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit phase (°)');

elseif datatype=='realimag',
subplot(2,1,1);
plot(sqrt(Texp),real(Sexp),'o',sqrt(Texp(Iact)),real(Sexp(Iact)),'o',sqrt(Texp),real(Sstart),sqrt(Texp),real(Sfit));
%set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit real');
title(['fitmain.m ',date]);
subplot(2,1,2);
plot(sqrt(Texp),imag(Sexp),'o',sqrt(Texp(Iact)),imag(Sexp(Iact)),'o',sqrt(Texp),imag(Sstart),sqrt(Texp),imag(Sfit));
%set(gca,'Xscale',xscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit imag');


end;%end if datatype==

end;



fprintf('fit parameter index array: %s\n',int2str(Ifitpar));
fprintf('fit parameters: %s\n',num2str(fitresult));
fprintf('start parameters: %s\n',num2str(startparam(Ifitpar)));
fprintf('fixed parameter index array: %s\n',int2str(Ifixpar));
fprintf('fixed parameters: %s\n',num2str(startparam(Ifixpar)));


elseif m==12,
%error analysis
ipar=min(length(Ifitpar),ipar);
mm=0;
while mm~=10,
 mm=menu('chi2 error analysis (sdv=(interval between chi2-minimum and chi2-doubling)/sqrt(N))',...
         ['parameter (1..',int2str(length(Ifitpar)),') : ',int2str(ipar)],...
         ['range for parameter: ',num2str(parmin),' to ',num2str(parmax)],...
         'calculate and plot error analysis',...
         ['best fit: ',num2str(fitresult(ipar))],...
         ['number of test points: ',int2str(ntest)],...
         'save settings to chi2err.par',...
         'load settings from chi2err.par',...
         'set parameter range to 50%-150%',...
         ['total number of data: ',int2str(length(Texp)),' number of data in fit: ',int2str(i2-i1+1)],...
         'return');
 if mm==1,
   fprintf('fit parameter: %i  \n',ipar);
   ipar=input('new parameter number  ? ');
   elseif mm==2,
   parmin=input('minimum value for fit parameter? ');
   parmax=input('maximum value for fit parameter? ');
   elseif mm==8,
   parmin=0.5*fitresult(ipar);
   parmax=1.5*fitresult(ipar);
   elseif mm==3,
   param=startparam;
   param(Ifitpar)=fitresult;
   Sstart=eval([functionstring,'(Texp,param)']);
   chi2vec=[];
   chi2vecext=[];
   h=figure; 
   counter=0;
   parvec=linspace(parmin,parmax,ntest)';
   chitestvec=[];
   signalvec=[];
   familyvec=[];
   chihistory=[];
   for ii=1:ntest,
    partest=parvec(ii);
    fitresulttest=fitresult;
    fitresulttest(ipar)=partest;
    chitestvec=[chitestvec;chimain(fitresulttest,startparam(Ifixpar))];
    signalvec=[signalvec,Sact];
    end;
  
   if xscalestring=='lin',
      Texp1=linspace(min(Texp),max(Texp),200)';
      else,
      Texp1=logspace(log10(min(Texp)),log10(max(Texp)),200)';
      end;
   Tvec=Texp1*ones(1,ntest);
  save chiparabhistory.txt chihistory -ascii -tabs


   for ii=1:ntest,
    partest=parvec(ii);
    param=startparam;
    param(Ifitpar(ipar))=partest; %variable fitting parameters
    S=eval([functionstring,'(Texp1,param)']);
    familyvec=[familyvec,Sact];
    end;
   
   if length(Sstart)~=length(Sexp),,
    Sstart=Sexp;
    end;
   if length(Sstart)~=length(Sfit),
    Sfit=Sstart;
    end;
   [P,S]=polyfit(parvec,chitestvec,2);
   parabool=polyval(P,parvec);
   sigpar=sqrt(abs(-P(2)^2+4*P(1)*P(3)))/(2*P(1));
   parfitmin=-P(2)/(2*P(1));
   fprintf('par(%i): %4.3e +/- %4.3e/sqrt(%i)\n',ipar,fitresult(ipar),sigpar,i2-i1+1);

   figure;
   comstr=['par. ',int2str(ipar),' varies between ',num2str(parmin),' and ',num2str(parmax)]; 
   if datatype=='real    ',
    loglog(Tvec,familyvec);
    set(gca,'Xscale',xscalestring);
    set(gca,'Yscale',yscalestring);
    zoom on;
    xlabel('x-variable');
    ylabel('S');
    title(['fitmain.m ',date,' ',comstr]);
    elseif datatype=='complex ',
    subplot(2,1,1);
    loglog(Tvec,abs(familyvec));
    set(gca,'Xscale',xscalestring);
    set(gca,'Yscale',yscalestring);
    zoom on;
    xlabel('x-variable');
    ylabel('amp');
    title(['fitmain.m ',date,' ',comstr]);
    subplot(2,1,2);
    semilogx(Tvec,angle(familyvec)*180/pi);
    set(gca,'Xscale',xscalestring);
    zoom on;
    xlabel('x-variable');
    ylabel('phase (°)');
    
    elseif datatype=='realimag',
    subplot(2,1,1);
    loglog(Tvec,real(familyvec));
    set(gca,'Xscale',xscalestring);
    set(gca,'Yscale',yscalestring);
    zoom on;
    xlabel('x-variable');
    ylabel('real');
    title(['fitmain.m ',date,' ',comstr]);
    subplot(2,1,2);
    semilogx(Tvec,imag(familyvec));
    set(gca,'Xscale',xscalestring);
    zoom on;
    xlabel('x-variable');
    ylabel('imag');
    
    
    end;%end if datatype==



   figure;
   plot(parvec,chitestvec,'o',parvec,parabool);
   title(['fitmain.m ',date,' par(',int2str(ipar),'): ',num2str(fitresult(ipar)),'+/-',num2str(sigpar),'/sqrt(',int2str(i2-i1+1),') file: ',name]);
   zoom on;
   figure;
   if datatype=='real    ',
    loglog(Texp,Sexp,'o',Texp(Iact),Sexp(Iact),'o',Texp,Sstart,Texp,Sfit,Tvec,familyvec);
    set(gca,'Xscale',xscalestring);
    set(gca,'Yscale',yscalestring);
    zoom on;
    xlabel('x-variable');
    ylabel('exp,startfit,resultfit,trial curves');
    title(['fitmain.m ',date,' ',comstr]);
    elseif datatype=='complex ',
    subplot(2,1,1);
    loglog(Texp,abs(Sexp),'o',Texp(Iact),abs(Sexp(Iact)),'o',Texp,abs(Sstart),Texp,abs(Sfit),Tvec,abs(familyvec));
    set(gca,'Xscale',xscalestring);
    set(gca,'Yscale',yscalestring);
    zoom on;
    xlabel('x-variable');
    ylabel('exp,startfit,resultfit amp, trial curves amp');
    title(['fitmain.m ',date,' ',comstr]);
    subplot(2,1,2);
    semilogx(Texp,angle(Sexp)*180/pi,'o',Texp(Iact),angle(Sexp(Iact))*180/pi,'o',Texp,angle(Sstart)*180/pi,Texp,angle(Sfit)*180/pi,Tvec,angle(familyvec)*180/pi);
    set(gca,'Xscale',xscalestring);
    zoom on;
    xlabel('x-variable');
    ylabel('exp,startfit,resultfit phase (°), trial curves phase');
    
    elseif datatype=='realimag',
    subplot(2,1,1);
    loglog(Texp,real(Sexp),'o',Texp(Iact),real(Sexp(Iact)),'o',Texp,real(Sstart),Texp,real(Sfit),Tvec,real(familyvec));
    set(gca,'Xscale',xscalestring);
    set(gca,'Yscale',yscalestring);
    zoom on;
    xlabel('x-variable');
    ylabel('exp,startfit,resultfit real, trial curves real');
    title(['fitmain.m ',date,' ',comstr]);
    subplot(2,1,2);
    semilogx(Texp,imag(Sexp),'o',Texp(Iact),imag(Sexp(Iact)),'o',Texp,imag(Sstart),Texp,imag(Sfit),Tvec,imag(familyvec));
    set(gca,'Xscale',xscalestring);
    zoom on;
    xlabel('x-variable');
    ylabel('exp,startfit,resultfit imag, trial curves imag');
    
    end;%end if datatype==
   





   elseif mm==5,
   ntest=input('number of test points? ');
   elseif mm==6,
   %save settings
   A=[ipar;parmin;parmax;ntest];
   save chi2err.par A -ascii -tabs
   elseif mm==7,
      load chi2err.par
      A=chi2err;
      ipar=A(1);parmin=A(2);parmax=A(3);ntest=A(4);
   end;%end if mm
 end;%end while mm


elseif m==13,
ipar1=min(length(Ifitpar),ipar1);
ipar2=min(length(Ifitpar),ipar2);

%contour error analysis
mm=0;
while mm~=10,
 

 mm=menu('contour error analysis',...
         ['parameters (1..',int2str(length(Ifitpar)),') : ',' first: ',int2str(ipar1),' second: ',int2str(ipar2)],...
         ['range for first parameter: ',num2str(par1min),' to ',num2str(par1max)],...
         ['range for second parameter: ',num2str(par2min),' to ',num2str(par2max)],...
         'calculate and plot contour',...
         ['best fit for first parameter: ',num2str(fitresult(ipar1)),' best fit for second parameter: ',num2str(fitresult(ipar2))],...
         ['number of test points: ',int2str(ntest)],...
         'save settings to conterr.par',...
         'load settings from conterr.par',...
         'set parameter range to 50%-150%',...
         'return');
 if mm==1,
   fprintf('fit parameters: %i and %i \n',ipar1,ipar2);
   ipar1=input('new first parameter number  ? ');
   ipar2=input('new second parameter number  ? ');

   elseif mm==2,
   par1min=input('minimum value for first fit parameter? ');
   par1max=input('maximum value for first fit parameter? ');

   elseif mm==3,
   par2min=input('minimum value for second fit parameter? ');
   par2max=input('maximum value for second fit parameter? ');

   elseif mm==9,
   par1min=0.5*fitresult(ipar1);
   par1max=1.5*fitresult(ipar1);
   par2min=0.5*fitresult(ipar2);
   par2max=1.5*fitresult(ipar2);

   elseif mm==4,
   param=startparam;
   param(Ifitpar)=fitresult;
   Sstart=eval([functionstring,'(Texp,param)']);
   chi2vec=[];
   chi2vecext=[];

   counter=0;
   fitvec1=linspace(par1min,par1max,ntest)';
   fitvec2=linspace(par2min,par2max,ntest)';
   h=figure;
   chitestvec=[];
   chihistory=[];
   for i=1:ntest,
    for j=1:ntest,
    fitresulttest=fitresult;
    fitresulttest(ipar1)=fitvec1(i);
    fitresulttest(ipar2)=fitvec2(j);
    chitestvec(i,j)=chimain(fitresulttest,startparam(Ifixpar));
    end;%end for j
    end;%end for i
   figure;
   contour(fitvec1,fitvec2,log10(chitestvec)');
   title(['fitmain.m ',date,' ipar1: ',int2str(ipar1),' ipar2: ',int2str(ipar2),' file: ',name]);
   xlabel(['parameter 1 ',int2str(ipar1)]);
   ylabel(['parameter 2 ',int2str(ipar2)]);
   zoom on;
   figure;
   mesh(fitvec1,fitvec2,log10(chitestvec)');
   title(['fitmain.m ',date,' ipar1: ',int2str(ipar1),' ipar2: ',int2str(ipar2),' file: ',name]);
   xlabel(['parameter 1 ',int2str(ipar1)]);
   ylabel(['parameter 2 ',int2str(ipar2)]);
   zoom on;
   save chicontour chitestvec -ascii -tabs
   save chicontourhistory.txt chihistory -ascii -tabs

   elseif mm==6,
   ntest=input('number of test points? ');
   elseif mm==7,
   %save settings
   A=[ipar1;ipar2;par1min;par1max;par2min;par2max;ntest];
   save conterr.par A -ascii -tabs
   elseif mm==8,
      load conterr.par
      A=conterr;
   ipar1=A(1);ipar2=A(2);par1min=A(3);par1max=A(4);par2min=A(5);par2max=A(6);ntest=A(7);
   end;%end if mm
 end;%end while mm


elseif m==14,
xscale=1-xscale;
if xscale==0,
 xscalestring='lin';
 else,
 xscalestring='log';
 end;
        
elseif m==15,   
yscale=1-yscale;
if yscale==0,
 yscalestring='lin';
 else,
 yscalestring='log';
 end;

elseif m==16,
fprintf('Default fitfunction: funmain\n');
fprintf('Current fitfunction: %s\n',functionstring);
functionstring=input('Give new name for fitfunction ','s');

elseif m==17,
nplot=input('plot during fit every how many iterations? '); 
iplotfit=1-iplotfit;


elseif m==18,
%error analysis with complete iteration for every point of chi2 parabole
ipar=min(length(Ifitpar),ipar);
mm=0;
while mm~=8,
 mm=menu('chi2 error analysis with complete iteration over other variables for each parameter value',...
         ['parameter (1..',int2str(length(Ifitpar)),') : ',int2str(ipar)],...
         ['range for parameter: ',num2str(parmin),' to ',num2str(parmax)],...
         'calculate and plot error analysis',...
         ['best fit: ',num2str(fitresult(ipar))],...
         ['number of test points: ',int2str(ntest)],...
         'save settings to chi2err.par',...
         'load settings from chi2err.par',...
         'return');
 if mm==1,
   fprintf('fit parameter: %i  \n',ipar);
   ipar=input('new parameter number  ? ');
   elseif mm==2,
   parmin=input('minimum value for fit parameter? ');
   parmax=input('maximum value for fit parameter? ');
   elseif mm==3,
   param=startparam;
%   param(Ifitpar)=fitresult;
   Sstart=eval([functionstring,'(Texp,param)']);
   chi2vec=[];
   chi2vecext=[];
   h=figure; 
   counter=0;
   parvec=linspace(parmin,parmax,ntest)';
   chitestvec=[];
   signalvec=[];
   familyvec=[];
   nn=length(Ifitpar);
   ii=find((1:nn)-ipar);
   Ifitparnew=Ifitpar(ii);
   Ifixparnew=[Ifixpar,Ifitpar(ipar)];
   startparamnew=startparam;
   startparamnew(Ifitparnew)=param(Ifitparnew);
   result=[];
   chilimhistory=[];
   for i=1:ntest,
    partest=parvec(i);
    startparamnew(Ifitpar(ipar))=partest;
    fprintf('testing chi2 sensitivity to parameter %i, iteration %i of %i\n',ipar,i,ntest);
    fprintf('minimizing chi2 for this value of parameter, varying all other fit parameters\n');
    foptions=optimset;
    fitresulttest=fminsearch('chimainlim',startparamnew(Ifitparnew),foptions,startparamnew(Ifixparnew));
    param=startparam;
    param(Ifitparnew)=fitresulttest;
    param(ipar)=partest;
    chinow=chimainlim(fitresulttest,param(Ifixparnew));
    result=[result;[param,chinow]];
    chitestvec=[chitestvec;chinow];
    signalvec=[signalvec,Sact];
    end;%end for i
   fprintf('history of fit ([parameters,chi-value] is saved in file iterparab.par\n');
   save chispecialhistory.txt chilimhistory -ascii -tabs
   fprintf('dependence of best fitting parameters and chi2 on fixed parameter (also in file iterparab.par): \n');
   result
   figure;
   plot(parvec*ones(1,length(Ifitparnew)),result(:,1:length(Ifitparnew)),'o-');
   ylabel('best fit parameters ');
   xlabel('fixed parameter value');
   zoom on;
   
   save iterparab.par result -ascii -tabs
   if xscale=='0',
      Texp1=linspace(min(Texp),max(Texp),200)';
      else,
      Texp1=logspace(log10(min(Texp)),log10(max(Texp)),200)';
      end;
   Tvec=Texp1*ones(1,ntest);

   
   if length(Sstart)~=length(Sexp),,
    Sstart=Sexp;
    end;
   if ~exist('Sfit'),
    Sfit=Sstart;
    end;

   if length(Sstart)~=length(Sfit),
    Sfit=Sstart;
    end;
   [P,S]=polyfit(parvec,chitestvec,2);
   parabool=polyval(P,parvec);
   sigpar=sqrt(abs(-P(2)^2+4*P(1)*P(3)))/(2*P(1));
   parfitmin=-P(2)/(2*P(1));
   fprintf('par(%i): %4.3e +/- %4.3e/sqrt(%i) \n',ipar,fitresult(ipar),sigpar,length(Iact));

   comstr=['par. ',int2str(ipar),' varies between ',num2str(parmin),' and ',num2str(parmax),' - each chi2 optimized']; 

   figure;
   plot(parvec,chitestvec,'o',parvec,parabool);
   title(['fitmain.m ',date,' par(',int2str(ipar),'): ',num2str(fitresult(ipar)),'+/-',num2str(sigpar),'/sqrt(',int2str(length(Iact)),') file: ',name]);
   ylabel(comstr);
   zoom on;
   figure;
   if datatype=='real    ',
    loglog(Texp,Sexp,'o',Texp(Iact),Sexp(Iact),'o',Texp,Sstart,Texp,Sfit);
    set(gca,'Xscale',xscalestring);
    set(gca,'Yscale',yscalestring);
    zoom on;
    xlabel('x-variable');
    ylabel('exp,startfit,resultfit');
    title(['fitmain.m ',date,' ',comstr]);
    elseif datatype=='complex ',
    subplot(2,1,1);
    loglog(Texp,abs(Sexp),'o',Texp(Iact),abs(Sexp(Iact)),'o',Texp,abs(Sstart),Texp,abs(Sfit));
    set(gca,'Xscale',xscalestring);
    set(gca,'Yscale',yscalestring);
    zoom on;
    xlabel('x-variable');
    ylabel('exp,startfit,resultfit amp');
    title(['fitmain.m ',date,' ',comstr]);
    subplot(2,1,2);
    semilogx(Texp,angle(Sexp)*180/pi,'o',Texp(Iact),angle(Sexp(Iact))*180/pi,'o',Texp,angle(Sstart)*180/pi,Texp,angle(Sfit)*180/pi);
    set(gca,'Xscale',xscalestring);
    zoom on;
    xlabel('x-variable');
    ylabel('exp,startfit,resultfit phase (°)');
    elseif datatype=='realimag',
    subplot(2,1,1);
    loglog(Texp,real(Sexp),'o',Texp(Iact),real(Sexp(Iact)),'o',Texp,real(Sstart),Texp,real(Sfit));
    set(gca,'Xscale',xscalestring);
    set(gca,'Yscale',yscalestring);
    zoom on;
    xlabel('x-variable');
    ylabel('exp,startfit,resultfit real');
    title(['fitmain.m ',date,' ',comstr]);
    subplot(2,1,2);
    semilogx(Texp,imag(Sexp),'o',Texp(Iact),imag(Sexp(Iact)),'o',Texp,imag(Sstart),Texp,imag(Sfit));
    set(gca,'Xscale',xscalestring);
    zoom on;
    xlabel('x-variable');
    ylabel('exp,startfit,resultfit imag');
    end;%end if datatype==
   

   elseif mm==5,
   ntest=input('number of test points? ');
   elseif mm==6,
   %save settings
   A=[ipar;parmin;parmax;ntest];
   save chi2err.par A -ascii -tabs
   elseif mm==7,
      load chi2err.par
      A=chi2err;
      ipar=A(1);parmin=A(2);parmax=A(3);ntest=A(4);
   end;%end if mm
 end;%end while mm

elseif m==19,
isqrt=1-isqrt;

elseif m==20,
ifitstring=mod(ifitstring+1,6);



end;%end if m==
end;%end while m~=

