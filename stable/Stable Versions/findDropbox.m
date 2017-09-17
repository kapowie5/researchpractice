function [ Dboxpath ] = findDropbox( ~ )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
laptop='C:\Users\Troy\Dropbox\';
desktop='C:\Users\User\Dropbox\';
    A=exist(laptop,'dir');
    B=exist(desktop,'dir');
    
    if A==7
        Dboxpath=laptop;
    else if B==7
        Dboxpath=desktop;
    else
        msgbox('Dropbox Folder doesn''t exist')
        Dboxpath=uigetdir('C:\Users\');
    end

end