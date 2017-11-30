function [] = f_Interp_Gen()
outputDir = 'Intrp';
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

%% Locate Files
myDirTFs = './Data/TTFs/';
myFilesTFs = dir(fullfile(myDirTFs,'TFs-DEGs-time-Activation*.csv'));

myDirDEGs = './Data/TDEGs/';
myFilesDEGs = dir(fullfile(myDirDEGs,'DEGs-time-Activation*.csv'));

%% Interpolation
% interpolate expression profiles into outputDir/interp folder
for i = 1 : length(myFilesTFs)
    f_interp(sprintf('%s%s',myDirTFs,myFilesTFs(i).name), outputDir);
end

for i = 1 : length(myFilesDEGs)
    f_interp(sprintf('%s%s',myDirDEGs,myFilesDEGs(i).name), outputDir);
end

end

