function chi2=chimainlim(fitparam,fixparam);
global Texp Sexp counter Sstart h chi2vec n0 Ifitpar Ifixpar i1 i2 chi2vecext datatype xscalestring yscalestring functionstring ifitstring
global Ifitparnew Ifixparnew
global iplotfit nplot
%last change 25-1-06

global chihistory chilimhistory

counter=counter+1; %counts number of iterations
param=zeros(1,n0);
param(Ifitparnew)=fitparam; %variable fitting parameters
param(Ifixparnew)=fixparam; %fixed parameters

param

Sfit=eval([functionstring,'(Texp,param)']);

chi2=sum(real(Sfit(i1:i2)-Sexp(i1:i2)).^2)/length(Sexp(i1:i2)); %chi2 to be minimized
chi2ext=sum(real(Sfit-Sexp).^2)/length(Sexp); %chi2 on complete data range
%fprintf('fitting parameter values: %e\n',fitparam);
chi2vec=[chi2vec;chi2];
chi2vecext=[chi2vecext;chi2ext];
if ((iplotfit==0)|(mod(counter-1,nplot)~=0)),
fprintf('toggle iplot to see fit evolution - current iteration: no : %i  chi2: %e\n',counter,chi2);
else,
fprintf('toggle iplot to see fit evolution - current iteration: no : %i  chi2: %e\n',counter,chi2);
figure(h);
if datatype=='complex ',
subplot(3,1,1);
loglog(Texp,abs(Sexp),'o',Texp,abs(Sstart),Texp,abs(Sfit),Texp(i1:i2),abs(Sfit(i1:i2)));
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit amplitude');
title(['fitmain.m ',date,' fit iteration: ',int2str(counter),'  chi2: ',num2str(chi2)]);
subplot(3,1,2);
semilogx(Texp,angle(Sexp)*180/pi,'o',Texp,angle(Sstart)*180/pi,Texp,angle(Sfit)*180/pi,Texp(i1:i2),angle(Sfit(i1:i2))*180/pi);
set(gca,'Xscale',xscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit phase');
subplot(3,1,3);
semilogy(1:length(chi2vec),chi2vec,1:length(chi2vecext),chi2vecext);
xlabel('iteration');
ylabel('chi2 fit data and chi2 all data');
zoom on;

elseif datatype=='realimag',
subplot(3,1,1);
loglog(Texp,real(Sexp),'o',Texp,real(Sstart),Texp,real(Sfit),Texp(i1:i2),real(Sfit(i1:i2)));
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit amplitude');
title(['fitmain.m ',date,' fit iteration: ',int2str(counter),'  chi2: ',num2str(chi2)]);
subplot(3,1,2);
semilogx(Texp,imag(Sexp),'o',Texp,imag(Sstart),Texp,imag(Sfit),Texp(i1:i2),imag(Sfit(i1:i2)));
set(gca,'Xscale',xscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit phase');
subplot(3,1,3);
semilogy(1:length(chi2vec),chi2vec,1:length(chi2vecext),chi2vecext);
xlabel('iteration');
ylabel('chi2 fit data and chi2 all data');
zoom on;


elseif datatype=='real    ',

subplot(2,1,1);
loglog(Texp,(Sexp),'o',Texp,(Sstart),Texp,(Sfit),Texp(i1:i2),(Sfit(i1:i2)));
set(gca,'Xscale',xscalestring);
set(gca,'Yscale',yscalestring);
zoom on;
xlabel('x-variable');
ylabel('exp,startfit,resultfit ');
title(['fitmain.m ',date,' fit iteration: ',int2str(counter)]);
subplot(2,1,2);
semilogy(1:length(chi2vec),chi2vec,1:length(chi2vecext),chi2vecext);
xlabel('iteration');
ylabel('chi2 fit data and chi2 all data');
zoom on;
end;%end if datatype==
end;%end if iplotfit==
pause(0.1);


chilimhistory=[chilimhistory;[param,chi2,chi2ext]];
