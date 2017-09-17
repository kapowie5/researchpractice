function []=drawnet(w1,w2,CancelVal,instr,outstr)
%  DRAWNET
%  -------
%          Draw a two layer neuarl network.
%
%          drawnet(W1,W2,CancelVal,instring,outstring) draws the network
%          specified by the weights in W1 and W2. Positive weights are
%          represented by a solid line while a dashed line represents
%          a negative weight. Only weights and biases larger than 'CancelVal'
%          are drawn. A bias is represented by a vertical line through the
%          neuron.
%
%  INPUT:
%  W1       : Weights between input and hidden layer. The matrix structure is
%             [(# of hidden units)  *  (inputs + 1)] (the 1 is due to the bias)
%  W2       : Weights between hidden layer and output
%             [(outputs)  *  (# of hidden units + 1)]
%  CancelVal: Draw only weight/bias if it exceeds this value
%  instring : A "string matrix" with as many rows as there are inputs
%             (OPTIONAL). If present it labels the network inputs. Otherwise
%             they are simply numbered.
%  outstring: A "string matrix" with as many rows as there are outputs
%             (OPTIONAL and only if instr exists). If present it labels the
%             network outputs
%
%  See also OBDPRUNE, OBSPRUNE, NNPRUNE

%  Original function programmed by Claus Svarer, EI/CONNECT. Current version is
%  modified by Magnus Norgaard, IAU/EI/IMM, Technical University of Denmark.
%  LastEditDate: July 17, 1995

[N1,N0]=size(w1);
N0=N0-1;
[N2,dummy]=size(w2);
MaxNeu=max([N0 N1 N2]);
cla
LengthTres=0.025*MaxNeu;
axis([-0.1 2.1 0.5 MaxNeu+0.5]);
axis('off')

hold on
for i = 1:N0,
   plot(0,(MaxNeu/(N0+1))*i,'o');
   if nargin<=3,
     text(-0.1,(MaxNeu/(N0+1))*i-0.0,sprintf('%g',i));
   else
     text(-0.4,(MaxNeu/(N0+1))*i-0.0,sprintf(instr(i,:)));
   end
end;
for i = 1:N1,
   plot(1,(MaxNeu/(N1+1))*i,'o');
   if (w1(i,N0+1) ~= 0)   
      plot([1 1],[((MaxNeu/(N1+1))*i-LengthTres) ((MaxNeu/(N1+1))*i+LengthTres)]);
   end;
end;
for i = 1:N2,
   plot(2,(MaxNeu/(N2+1))*i,'o');
   if (w2(i,N1+1) ~= 0)   
      plot([2 2],[((MaxNeu/(N2+1))*i-LengthTres) ((MaxNeu/(N2+1))*i+LengthTres)]);
   end;
   if nargin==5,
      text(2.05,(MaxNeu/(N2+1))*i-0.0,sprintf(outstr(i,:)));
   end  
end;
%
MaxColorNo=7;
for i=1:N0,
   for j=1:N1,
      colour_int = ceil(abs(w1(j,i))*2);
      if (colour_int > MaxColorNo),
         colour_int = MaxColorNo;
      end;
      if (w1(j,i) > 0),
         colour=sprintf('-c%g',colour_int);
      else
         colour=sprintf('-.c%g',colour_int);
      end;
      if (abs(w1(j,i)) >= CancelVal),
         plot([0 1],[(MaxNeu/(N0+1))*i (MaxNeu/(N1+1))*j],'b');
      end;
   end;
end
for i=1:N1,
   for j=1:N2,
      colour_int = ceil(abs(w2(j,i))*2);
      if (colour_int > MaxColorNo),
         colour_int = MaxColorNo;
      end;
      if (w2(j,i) > 0),
         colour=sprintf('-c%g',colour_int);
      else
         colour=sprintf('-.c%g',colour_int);
      end;
      if (abs(w2(j,i)) >= CancelVal),
         plot([1 2],[(MaxNeu/(N1+1))*i (MaxNeu/(N2+1))*j],'b');
      end;
   end;
end
hold off


