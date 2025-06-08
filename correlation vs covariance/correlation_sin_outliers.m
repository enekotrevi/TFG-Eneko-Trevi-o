clear all
close all

% Load Data
% -------------------------------------------------------------------------
DATA = xlsread('Five_times_Kinematics_HT_Query_final.xlsx');

MF = DATA(:,2);
X_lin = DATA(:,5:end)';

% -------------------------------------------------------------------------
% a) Generate quadratic kernel
% -------------------------------------------------------------------------
X = X_lin;
lin_dim = size(X,1);
for i = 1:lin_dim
    X(lin_dim+i, :) = X(i, :).^2;
end

% Add only the interactions between adjacent points of the same variable
quad_dim = size(X,1);
for i = 1:lin_dim
    if mod(i,5) == 0
        i = i + 1;
    end
    X(quad_dim+i, :) = X(i,:) .* X(i+1,:);
end

% Remove zero rows
rowsToKeep = [];
for i = 1:size(X,1)
    if sum(abs(X(i,:))) ~= 0
        rowsToKeep = [rowsToKeep; i];
    end
end
X = X(rowsToKeep, :);

% Also remove the last row (artifact from meaningless multiplication)
X_dim = size(X,1);
X(X_dim,:) = [];

% -------------------------------------------------------------------------
% Determine the number of significant components
% -------------------------------------------------------------------------

% PCA baseline (using correlation matrix instead of covariance matrix)
% -------------------------------------------------------------------------
S = corrcoef(X');
[V, D] = eig(S);
n = size(X, 2);
dim = size(X, 1);

sort_EVal = sort(diag(D), 'descend');
TotVar = sum(sort_EVal);
Var = 0;
% Identify number of PCs explaining at least 95% of total variance
for i = 1:dim
    Var = Var + sort_EVal(i);
    Var99 = Var/TotVar;
    if Var99 > 0.95
        break
    end
end
n_Dimensions = i;

% -------------------------------------------------------------------------
% Select an appropriate statistic
% -------------------------------------------------------------------------
p = n_Dimensions;
[nada, s_idxs] = sort(diag(D), 'descend');
Vp = V(:, s_idxs(1:p));
Y = Vp' * (X - repmat(mean(X',1)', [1 n]));
