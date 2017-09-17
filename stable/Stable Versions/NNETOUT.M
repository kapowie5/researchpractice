function [y2]=nnetout(NetDef,W1,W2,PHI)
%  berekening van getraind neuraal netwerk met gewichten W1,W2 voor 
%     spectra van PHI, resultaat wordt in y2 gezet
%  ----
%          laatste wijziging: 24-3-96

%  INPUT:
%  NetDef  : Network definition
%  W1      : Weights between input and hidden layer. The matrix structure is
%            [(# of hidden units)  *  (inputs + 1)]  (the 1 is due to the bias)
%  W2      : Weights between hidden layer and output
%            [(outputs)  *  (# of hidden units + 1)]
%  PHI     : Input matrix. Structure: [(inputs)  *  (# of data)]
%  Y       : Output data. [(outputs) *  (# of data)]: real data
%  trparms : Vector containing parameters associated with the training
%  itr     : Number of training data


%  OUTPUT:
%  y2       : neural network responsvector


% begin berekening van NN output voor gegevens in PHI

    L_hidden = find(NetDef(1,:)=='L')';     % Location of linear hidden neurons
    H_hidden = find(NetDef(1,:)=='H')';     % Location of tanh hidden neurons
    L_output = find(NetDef(2,:)=='L')';     % Location of linear output neurons
    H_output = find(NetDef(2,:)=='H')';     % Location of tanh output neurons
    [hidden,inputs] = size(W1);             % # of hidden units 
    [x,N]    = size(PHI);                      % # of outputs and # of data
    [outputs,x] = size(W2);
    PHI      = [PHI;ones(1,N)];             % Augment PHI with a row containing ones
    y1       = [zeros(hidden,N);ones(1,N)]; % Hidden layer outputs
    y2       = zeros(outputs,N);            % Network output
    h1 = W1*PHI;  
    y1(H_hidden,:) = pmntanh(h1(H_hidden,:));
    y1(L_hidden,:) = h1(L_hidden,:);    
    h2 = W2*y1;
    y2(H_output,:) = pmntanh(h2(H_output,:));
    y2(L_output,:) = h2(L_output,:);
   



