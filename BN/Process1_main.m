% f_main([1,2]);clc;close all;
% f_main([1,2,3]);clc;close all;
% f_main([1,2,3,4]);clc;close all;
% f_main([1,2,3,4,5]);clc;close all;
% f_main([5,6]);clc;close all;
%f_main([2,3,4,5]);clc;close all;
%f_main([3,4,5]);clc;close all;
%f_main([4,5]);clc;close all;
% target_12 = f_main([1,2]);clc;close all;

% trgt = f_main([1,2]);clc;close all;
% trgt = f_main([1,2,3]);clc;close all;
% trgt = f_main([1,2,3,4]);clc;close all;
% trgt = f_main([1,2,3,4,5]);clc;close all;
% trgt = f_main([5,6]);clc;close all;

tp_array_List = {[1,2],[1,2,3],[1,2,3,4],[1,2,3,4,5],[1,2,3,4,5,6]};
% for n_levels = [2,3,5,7,10,15,20]
%     for i = 1 : length(tp_array_List)   
%         tp_array = tp_array_List{i};
%         ntp = length(tp_array);
% 
%     %     n_levels = 3; %discretize level
%         outputDir = sprintf('Target_Analysis/Target_Analysis_nlevels%d_',n_levels);
% 
%         for j = 1 : ntp
%             outputDir = strcat(outputDir,'%d');
%             outputDir = sprintf(outputDir,tp_array(j));
%         end
% 
%         % Reconstruct Network
%         trgt = f_main(tp_array, n_levels);
% 
%         % Analyze and Record Score
%         f_trgt_analysis(trgt,outputDir);
% 
%     end
% end

for n_levels = 2%[2,3]%[2,3,4,5]%,7,10,15,20]
    for i = 1 : length(tp_array_List)   
        tp_array = tp_array_List{i};
        ntp = length(tp_array);


        % Reconstruct Network
        trgt = f_main(tp_array, n_levels, 'isInterp',0,'isTPA',0);

        % Analyze and Record Score
        outputDir = sprintf('Target_Analysis/Target_Analysis_nlevels%d',n_levels);
        f_trgt_analysis(trgt,outputDir);

    end
end
