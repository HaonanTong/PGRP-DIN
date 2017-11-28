function [] = f_trgt_analysis( trgt, Target_Folder )
%F_TRGT_ANALYSIS gives analysis on the output of DBN
% Output: For each target output 3 Tables
% ATxGxxxxx_#Parent.csv where # = 1, 2, 3
% | Target | Source | BDeu  |
% |   G1   |   G2   | score |
% |   ..   |   ..   |  ...  |
% |   ..   |   ..   |  ...  |

nTarget = length(trgt);
    
if ~exist(Target_Folder,'dir')
    mkdir(Target_Folder);
end 
% T_OneParent = table('VariableNames',{'Source' 'Target' 'BDeu'});
% T_TwoParent = table('VariableNames',{'Source1' 'Source2' 'Target' 'BDeu'});
% T_ThreeParent = table('VariableNames',{'Source1' 'Source2' 'Source3' 'Target' 'BDeu'});

for t = 1 : nTarget
    target = trgt{t};
    % Generate T_OneParent
    Target = [];
    Source = [];
    BDeu = [];
    for comb = 1 : size(target.BDeu_Memory,2)
        if isempty(target.BDeu_Memory{1,comb})
                continue;
        end
        tSource = target.PoPaGList{target.BDeu_Memory{1,comb}.Indx};
        tBDeu = target.BDeu_Memory{1,comb}.BDeu;
        tTarget = target.name;
        Source = [Source ; tSource];
        Target = [Target ; tTarget];
        BDeu = [BDeu; tBDeu]; 
    end
    T_OneParent = table(Target,Source,BDeu);
    T_OneParent = sortrows(T_OneParent,'BDeu','descend'); % highest score at the top

    writetable(T_OneParent,sprintf('./%s/%s_1Parent.csv',Target_Folder,target.name)...
        ,'WriteRowNames',true,'WriteVariableNames',true);
    
%     % Generate T_TwoParent
%     nPa = 2;
%     Target = [];
%     Source1 = [];
%     Source2 = [];
%     BDeu = [];
%     for comb = 1 : size(target.BDeu_Memory,2)
%         if isempty(target.BDeu_Memory{nPa,comb})
%                 continue;
%         end
%         
%         indx = target.BDeu_Memory{nPa,comb}.Indx;
%         tSource1 = target.PoPaGList{indx(nPa-1)};
%         tSource2 = target.PoPaGList{indx(nPa)};
%         tBDeu = target.BDeu_Memory{nPa,comb}.BDeu;
%         tTarget = target.name;
%         Source1 = [Source1 ; tSource1];
%         Source2 = [Source2 ; tSource2];
%         Target = [Target ; tTarget];
%         BDeu = [BDeu; tBDeu]; 
%     end
%     T_TwoParent = table(Target,Source1,Source2,BDeu);
%     T_TwoParent = sortrows(T_TwoParent,'BDeu','descend'); % highest score at the top
%     
%     writetable(T_TwoParent,sprintf('./%s/%s_2Parent.csv',Target_Folder,target.name)...
%             ,'WriteRowNames',true,'WriteVariableNames',true);
%      
%     % Generate T_ThreeParent
%     nPa=3;
%     Target = [];
%     Source1 = [];
%     Source2 = [];
%     Source3 = [];
%     BDeu = [];
%     
%     for comb = 1 : size(target.BDeu_Memory,2)
%         if isempty(target.BDeu_Memory{nPa,comb})
%                 continue;
%         end
% %         TSource = target.PoPaGList(target.BDeu_Memory{2,comb}.Indx);
%         
%         indx = target.BDeu_Memory{nPa,comb}.Indx;
%         tSource1 = target.PoPaGList{indx(nPa-2)};
%         tSource2 = target.PoPaGList{indx(nPa-1)};
%         tSource3 = target.PoPaGList{indx(nPa)};
%         tBDeu = target.BDeu_Memory{nPa,comb}.BDeu;
%         tTarget = target.name;
%         Source1 = [Source1 ; tSource1];
%         Source2 = [Source2 ; tSource2];
%         Source3 = [Source3 ; tSource3];
%         Target = [Target ; tTarget];
%         BDeu = [BDeu; tBDeu]; 
%     end
%     T_ThreeParent = table(Target,Source1,Source2,Source3,BDeu);
%     T_ThreeParent = sortrows(T_ThreeParent,'BDeu','descend'); % highest score at the top
%     
%     writetable(T_ThreeParent,sprintf('./%s/%s_3Parent.csv',Target_Folder,target.name)...
%             ,'WriteRowNames',true,'WriteVariableNames',true);
% end

end

