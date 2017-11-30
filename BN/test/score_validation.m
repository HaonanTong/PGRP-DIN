data = [ 0 1 0 1 0 1 0 1 0 1 0 0 0 0 0 1 1 1 1 1 1 0 1 0 1 0 1 0 1 0 ; ...
         1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 ]; ...
%          0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 
N = size(data,1);

% config = [0 1]';
% [ Nij ] = f_getNijk( config, data, 2 )

% Hyperparameters
ns = [2,2]; % node size for each node;
ESS = 3; % Equivalent Sample Size;

% BDeu Score for a specific structure;    
% dag = [ 0 0 1 ; 0 0 1; 0 0 0]; % struture to be evaluated;
dag = [ 0 1; 0 0];
params = cell(1,N);
for i=1:N
    params{i} = {'prior_type', 'dirichlet', 'dirichlet_weight', ESS};
end

bayes_score = score_dags(num2cell(data), ns, {dag}, 'params', params)
bayes_score2 = f_calculateScore(data,2,ESS)


% Angels
BDei = 0;
REG = data(1,:);
TAR = data(2,:);
n_levels = 2;
for r = 1:n_levels % # states of each target or regulator
    Nj = REG == r-1;
    BDei = BDei + log(gamma(ESS/n_levels)/gamma(sum(Nj)+ESS/n_levels));
    for d = 1:n_levels
        Njk = REG == r-1 & TAR == d-1;
        BDei = BDei + log(gamma(sum(Njk)+ESS/(n_levels*n_levels))/gamma(ESS/(n_levels*n_levels)));
    end
end


% BDeu Score for iteration of all structure;
% Illustration of Markov Equivalence;
dags = mk_all_dags(N);

params = cell(1,N);
for i=1:N
    params{i} = {'prior_type', 'dirichlet', 'dirichlet_weight', ESS};
end

bayes_score = score_dags(num2cell(data), ns, dags, 'params', params);