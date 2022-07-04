% a = table();
var1 = uint16([1 2 3]');
var2 = [true false false]';
var3 = [1.2 1.5 3.2]';
tab = table(var1,var2,var3);
v1 = uint16([1 2 3]');
v2 = [true false false]';
v3 = [1.2 1.5 3.2]';
tab2 = table(v1,v2,v3);
tab = [tab;tab2]
% 
% a.Variables = {var1,var2,var3};
% 
% 
% var1 = uint16([1 2 3 4]');
% var2 = [true false false, true]';
% var3 = [1.2 1.5 3.2 5.2]';
% 
% 
% a.Variables = {var1,var2,var3};


%     t = readtable('tsunamis.xlsx');
    vars = {'Latitude','Longitude','MaxHeight'};
    t2 = t(1:20,vars);
    t3 = t(1:2,[1 2 4]);
    t4 = t(1:2,:);
    
    
%     t5 = table2cell(t3)
%     t5{1,1} = 0;
%     t3
    t6 = cell2table(t5);
    
%     t3(1,1)