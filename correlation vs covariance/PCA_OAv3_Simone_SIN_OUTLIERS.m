clear all
close all

% Load Data
% -------------------------------------------------------------------------
DATA = xlsread('Five_times_Kinematics_HT_Query_final.xlsx');

MF = DATA(:,2);
X_lin= DATA(:,5:end)';

% -------------------------------------------------------------------------
% a) Generate quadratic kernel
% Given the relevance and the interpretations of the data, in the quadratic
% kernel only interaction of subsequent values are considered.
% -------------------------------------------------------------------------
X=X_lin;
lin_dim=size( X,1);
for i=1:lin_dim
    X(lin_dim+i,:)=X(i,:).^2;
end

%We now add only the interactions between adiacent points of the same
%variable

quad_dim=size( X,1);
for i=1:lin_dim
    if mod(i, 5) == 0
        i=i+1;
    end
    X(quad_dim+i,:)=X(i,:).*X(i+1,:);
end

%remove zero rows
% Initialize an empty array to store indices of rows to keep
rowsToKeep = [];
for i = 1:size(X, 1)
        % Calculate the sum of the absolute values of the elements in the row
        if sum(abs(X(i, :))) ~= 0
            % If the sum is not zero, add the row index to rowsToKeep
            rowsToKeep = [rowsToKeep; i];
        end
end
    
% Create the modified matrix by selecting only the rows to keep
X= X(rowsToKeep, :);
%also remove the last row which was obtained moltipling miningless indexes
X_dim=size( X,1);
X(X_dim,:) = [];



% -------------------------------------------------------------------------
% a) Determine the number of significant components
% -------------------------------------------------------------------------

% PCA baseline
% ------------
S = cov( X' );
%S = corrcoef( X' );
[V,D] = eig( S );
n = size( X, 2 );
dim=size( X, 1 );

sort_EVal = sort( diag( D ), 'descend' );
TotVar=sum(sort_EVal);
Var=0;
%identify sum of PCs up to 99% of the variance
for i=1 : dim
    Var=Var+sort_EVal(i);
    Var99=Var/TotVar;
    if Var99>0.95
        break
    end
end
n_Dimensions=i;




% -------------------------------------------------------------------------
% d) Select an appropriate statistic
% -------------------------------------------------------------------------
%p = n_signif_components 
p=n_Dimensions;
% Compute the significant PCs
[nada, s_idxs] = sort( diag( D ), 'descend');
Vp = V(:, s_idxs(1:p));
Y = Vp' * (X - repmat( mean( X' )', [1 n]) );









