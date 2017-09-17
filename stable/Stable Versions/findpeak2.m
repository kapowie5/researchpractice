function [PeakSpec,PeakwL]=findpeak2(wL,Spec,c,p0,pI0)
%function w=cutdc(A,dc0);
%dc0:threshold 
[a,b]=size(wL);
if a>1
    wL=wL';
end
[a2,b2]=size(Spec);
if a2>1
    Spec=Spec';
end
[Spect_sort,IS]=sort(Spec,2,'descend');
index=round(mean(IS(1:5)));
indexL=index-c;
indexR=index+c;
if indexR > b2;
    indexR = b2;
    flag = 1;
else
    flag = 0;
end
if indexL < a2;
    indexL = a2;
    flag = 1;
else
    flag = 0;
end
SpecPart=Spec(indexL:indexR);
wLPart=wL(indexL:indexR);
[p_weight,S,MU]=polyfit(wLPart,SpecPart,3);
newwL=linspace(wL(indexL),wL(indexR),100);
SpecFit=polyval(p_weight,newwL,[],MU);
[PeakSpec,PeakwLIndex]=max(SpecFit);
PeakwL=newwL(PeakwLIndex);
if flag==1
    PeakSpec = p0;
    PeakwL = pI0;
end
% figure(100)
% plot(wLPart,SpecPart,'ob',newwL,SpecFit,'-b')
end
