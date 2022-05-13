
[file,path] = uiputfile('saves\*.mat');
save(strcat(path,file),'a','b')


% [filename, path] = uigetfile('*.mat','Load save file','saves\');
% infile = importdata(strcat(path,filename));
% a = infile.a;
% b = infile.b;

% infile = importdata('C:\Users\Mikael\Dev\Refractive_index_approximation\saves\a.mat');