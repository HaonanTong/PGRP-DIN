function [] = f_discritize2( myDir_Interp, n_levels)
myFiles2 = dir(fullfile(myDir_Interp,'*.csv'));
% n_levels =10;

for i = 1 : length(myFiles2)
    csv = sprintf('%s/%s',myDir_Interp,myFiles2(i).name);
    f_discritize( csv, n_levels );
end

end

