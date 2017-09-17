function Y=filter_mean(X,nx);
%function Y=filter_mean(X,nx);
%first edition: 20-1-2007
%X=input signal
%Y=median filtered signal with mean area 2*nx+1 
Y=X*0;
n=length(X);
for ii=1:n,
  [i1,I]=max([ii-nx,1]);
  [i2,I]=min([ii+nx,n]);
  A=X(i1:i2);
  Y(ii)=mean(A);
  end;%end for jj


