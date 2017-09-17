% nncalc.m berekent getrainde parameter uit spectrum en gewichten
% par is berekende parameter
% spectrum is inputdata voor NN
% sstring bepaalt voor welke parameter gewichten w1 en w1 moeten ingeladen worden
% hoort bij programma nnexec.m
% laatste wijziging: 12-10-1997
function [par]=nncalc(spectrum,sstring,nr,units,outunit,ndata)

  %---------- Define network structure ----------
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
 
      w1naam=sprintf('w1%i.nn',nr);
      w2naam=sprintf('w2%i.nn',nr);
      %w1T en w2T worden opgevuld per kolom tijdens fscanf
      fid=fopen(w1naam,'r');
      [w1T,count]=fscanf(fid,'%e',[ndata+1,units]);
      fclose(fid);
      w1=w1T';
      fid=fopen(w2naam,'r');
      [w2T,count]=fscanf(fid,'%e',[units+1,1]);
      w2=w2T';
      fclose(fid);
    %------ initialize network

    L_hidden = find(NetDef(1,:)=='L')';     % Location of linear hidden neurons
    H_hidden = find(NetDef(1,:)=='H')';     % Location of tanh hidden neurons
    L_output = find(NetDef(2,:)=='L')';     % Location of linear output neurons
    H_output = find(NetDef(2,:)=='H')';     % Location of tanh output neurons
    [hidden,inputs] = size(w1);             % # of hidden units

    % --- bereken NN output uit spectrum, w1 en w2

    PHI0     = spectrum;                    % NN inputs
    outputs  = 1;                           % # of outputs 
    [x,N]    = size(PHI0);         	    % # of data
    PHI0     = [PHI0;ones(1,N)];            % Augment PHI with a row containing ones
    y1       = [zeros(hidden,N);ones(1,N)]; % Hidden layer outputs
    y2       = zeros(outputs,N);            % Network output
    h1 = w1*PHI0;
    y1(H_hidden,:) = pmntanh(h1(H_hidden,:));
    y1(L_hidden,:) = h1(L_hidden,:);    
    h2 = w2*y1;
    y2(H_output,:) = pmntanh(h2(H_output,:));
    y2(L_output,:) = h2(L_output,:);
    par=y2;
