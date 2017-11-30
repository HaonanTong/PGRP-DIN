function [target] = f_main( tp_array, n_levels, varargin )
% Example
% tp_array = [ 1, 2, 3, 4]
% for each nTF activated at time point 4 as target
% TF activated at time point 1, 2, 3, 4 as potential regulators

% Default Setting
isInterp = 1;
isInterpGen = 0;
isTPA = 1;


args = varargin;
nargs = length(args);
for i=1:2:nargs
    switch args{i}
        case 'isInterp', isInterp = args{i+1};
        case 'isInterpGen', isInterpGen = args{i+1};
        case 'isTPA', isTPA = args{i+1};
    end
end

if isInterpGen
    f_Interp_Gen
end

% Generate Folder
ntp = length(tp_array);
tartp = tp_array(end);
if ntp == 1 
    pRegtp_array = 1:tartp; % potential regulators are activated by tartp
    % if tp_array = 6, potential TF: 1 2 3 4 5
else
    pRegtp_array = tp_array(1:ntp);
    % if tp_array = [5 6], potential TF: [5, 6]
end

T_array = [.25, .5, 1 , 4, 12, 24 ];

%% Discretize
if isInterp
    myDscDir = 'Intrp';
else
    myDscDir = 'nIntrp';
end

f_discritize2(myDscDir, n_levels);

%% Read Table
myDir_Dsc = strcat(myDscDir,'/Dscrtz/');
T_table = cell(1,6);
TF_at = cell(1,6);
DEGs_at = cell(1,6);

for i = 1 : 6 %TFs
    % T4 TFs-DEGs-time-Activation4-up-Interp-Dscrtz.csv
%     Ethylene-nTFs-DEGs-time-Activation1-Interp-Dscrtz.csv
    csv = sprintf('%s%s',myDir_Dsc,sprintf('TFs-DEGs-time-Activation%d-Interp-Dscrtz.csv',i));
    
    if ~isInterp
        csv = sprintf('%s%s',myDir_Dsc,sprintf('TFs-DEGs-time-Activation%d-Dscrtz.csv',i));
    end
    
    T_table{i} = readtable(csv,'ReadRowNames',true,'ReadVariableNames',true);

    TF_at{i}.tp = i; %T4
    TF_at{i}.it = T_array(i)/.25 + 1; %IT17;
    TF_at{i}.dmtrx = table2array(T_table{i}); 
    TF_at{i}.table = T_table{i}; 
    TF_at{i}.ngenes = size( TF_at{i}.dmtrx , 1 );
    TF_at{i}.glist = T_table{i}.Properties.RowNames;
end

for i = 1 : 6 %DEGs
    % T4 TFs-DEGs-time-Activation4-up-Interp-Dscrtz.csv
    csv = sprintf('%s%s',myDir_Dsc,sprintf('DEGs-time-Activation%d-Interp-Dscrtz.csv',i));
    
    if ~isInterp
        csv = sprintf('%s%s',myDir_Dsc,sprintf('DEGs-time-Activation%d-Dscrtz.csv',i));
    end
    
    T_table{i} = readtable(csv,'ReadRowNames',true,'ReadVariableNames',true);

    DEGs_at{i}.tp = i; 
    DEGs_at{i}.it = T_array(i)/.25 + 1; 
    DEGs_at{i}.dmtrx = table2array(T_table{i}); 
    DEGs_at{i}.table = T_table{i}; 
    DEGs_at{i}.ngenes = size( DEGs_at{i}.dmtrx , 1 );
    DEGs_at{i}.glist = T_table{i}.Properties.RowNames;
end

% Derive PPaMtrx
PoPaGList = [];

if ~isTPA
    ntp = 6;
    pRegtp_array = 1 : 6;
end

for i = 1 : ntp
    PoPaGList = [PoPaGList; TF_at{pRegtp_array(i)}.glist];
end
target = cell(1,DEGs_at{tartp}.ngenes);
npRtp = length(pRegtp_array);
PPaMtrx = [];
for i = 1 : npRtp
    PPaMtrx = [PPaMtrx; TF_at{pRegtp_array(i)}.dmtrx];
end

% if npRtp == 1
%     PPaMtrx.matrix = TF_at{pRegtp_array}.dmtrx;
%     PPaMtrx.at = TF_at{pRegtp_array}.it;
% end

rmax = 1;
% ESS = 1E-10;
ESS = 1E-5;
for i = 1 : length(target) % for each potential target
    % alignment potential regulator
    Atarget = DEGs_at{tartp}.dmtrx(i,:);
%     B.matrix = Atarget;
%     B.at = TF_at{tartp}.it;
%     ADM = f_alignment(PPaMtrx, B, 0);

    ADM = [PPaMtrx; Atarget];
    % get topology
    [ topologyi, BDeu_Memory ] = f_getTopology( ADM, n_levels, rmax, ESS ); 
    target{i}.name = DEGs_at{tartp}.glist{i};
    target{i}.indx = topologyi.indx;
    target{i}.BDeu = topologyi.BDeu;
    target{i}.Pa = PoPaGList(target{i}.indx);
    target{i}.PoPaGList = PoPaGList;
    target{i}.BDeu_Memory = BDeu_Memory;
end


% Print BDeu Score




% 
% %% Analysis
% BDeu = zeros(DEGs_at{tartp}.ngenes,1);
% for i = 1 : DEGs_at{tartp}.ngenes % for each potential target
%     BDeu(i) = target{i}.BDeu;
% end
% 
% % BDeu_max = max(BDeu);
% % indx = find(abs(BDeu-BDeu_max) < .1 ); % corresponds to strongest causal relationship;
% % nl = length(indx); % # of structure s.t. has highest score.
% % DBN_target = cell(1,nl);
% % for i = 1 : nl
% %    DBN_target{1,i} = target{1,indx(i)};
% % end
% % fprintf('Totally have %d sub-structure\n',length(DBN_target));
% 
% % for i = 1 : length(DBN_target)
% %    fprintf('-----------------------\n');
% %    fprintf('structure: %d\n',i);
% %    fprintf('Regulators:\n');
% %    disp(DBN_target{i}.Pa);
% %    fprintf('Target: %s\n',DBN_target{i}.name);
% %    fprintf('BDeu score: %G\n',DBN_target{i}.BDeu);  
% % end
% 
% %% grep structure.
% structure = cell(1,length(target));
% for i = 1 : length(target)
%     structure{i} = [ target{i}.Pa ; target{i}.name ]; % last one is target;
% end
% 
% % output expression profiles.
% T = readtable('Profiles-ANan-DEGs-Filtered.csv','ReadVariableNames',true,'ReadRowNames',true);
% for i = 1 : length(target)
%     sub_T = T(structure{i},:);
%     
%     writetable(sub_T,sprintf('./%s/Sub_Table%d.csv',outputDir,i)...
%         ,'WriteRowNames',true,'WriteVariableNames',true);
% end
% 
% % Visulization
% % myFiles = dir(fullfile(outputDir,'/Sub_Table*.csv')); %gets all wav files in struct
% % for i = 1 : length(myFiles)
% %     csv = sprintf('%s/%s',outputDir,myFiles(i).name);
% %     f_plotTable3( csv, [], 'Normalized' );
% % end
% 
% %% output file for cytoscape
% % iteration 30 targets, output file for cytoscape
% ntargets = length(target);
% % T_cytoscape = table('VariableNames',{'Source' 'Target' 'BDeu'});
% Source = [];
% Target = [];
% BDeu = [];
% for i = 1 : ntargets
%     tSource = target{i}.Pa;
%     nPa = length(tSource);
%     tTarget = [];
%     for j = 1 : nPa
%         tTarget = [tTarget;cellstr(target{i}.name)];
%     end             
%     tBDeu = target{i}.BDeu * ones(nPa,1);    
%     Source = [Source ; tSource];
%     Target = [Target ; tTarget];
%     BDeu = [BDeu; tBDeu];
% end
% %%
% [TSource,~] = f_tranlate(Source,'./ORF_Translation/table_mapping.csv');
% [TTarget,~] = f_tranlate(Target,'./ORF_Translation/table_mapping.csv');
% T_cytoscape = table(Source,Target,BDeu,TSource,TTarget);
% T_cytoscape_sort = sortrows(T_cytoscape,'BDeu','descend'); % highest score at the top
% 
% writetable(T_cytoscape_sort,sprintf('./%s/Network_Cytoscape.csv',outputDir)...
%         ,'WriteRowNames',true,'WriteVariableNames',true);
% 

