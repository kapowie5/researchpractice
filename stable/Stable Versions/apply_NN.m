function [P]=apply_NN(param,shift,factor,signal_exp)
% paramFile='D:\liwang liu\PhD study\experiments\Fluore_exp_NN\withPump\param.nn';
% fid1=fopen(paramFile);
% [A,count]=fscanf(fid1,'%e',inf);
% fclose(fid1);
A=param;
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
scaling=0;
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
dimE=size(signal_exp,1);

theta=[signal_exp,ones(dimE,i4-i3+1)];
clear r c i

PHI=theta(:,i1:i2)';
clear m sd;
%%%%%%%%%%%%%%%%%%%%important
%%added on 22th Feb. 2011
scaling=1;
%%added on 22th Feb. 2011
%%%%%%%%%%%%%%%%%%%%important
% if scaling==0,
%%deleted
% else,
load m.nn
load sd.nn
% load factor.nn
% shift=load('shift.nn');
load trainf.nn
Y=theta(:,i3:i4)';
for i=1:(i4-i3+1)
    Y(i,:)=factor(i)*Y(i,:)+shift(i);
end;
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
P=thetam;
end

