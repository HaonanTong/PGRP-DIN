function [ expr_n_mv ] = f_normalize( expr )
% zero mean standard deviation

[~, tp] = size(expr);

expr_bar = mean(expr,2);
expr_n_m = expr - repmat(expr_bar,1,tp);
expr_n_mv = expr_n_m./repmat(sqrt(var(expr_n_m,0,2)),1,tp);


end

