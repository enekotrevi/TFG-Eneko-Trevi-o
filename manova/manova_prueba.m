clear; close all; clc;

%% -------------------- CARGAR DATOS --------------------
DATA = xlsread('Five_times_Kinematics_HT_Query_MATLAB_SIN_OUTLIERS.xlsx');
MF = DATA(:, 2);                 % Etiquetas de grupo
X_lin = DATA(:, 3:end)';         % Datos de entrada (transpuestos)

%% -------------------- KERNEL CUADRÁTICO --------------------
X = X_lin;
lin_dim = size(X, 1);

% Añadir cuadrados
for i = 1:lin_dim
    X(lin_dim + i, :) = X(i, :).^2;
end

% Añadir interacciones adyacentes
quad_dim = size(X, 1);
for i = 1:lin_dim
    if mod(i, 5) == 0 || i == lin_dim
        continue;
    end
    X(quad_dim + i, :) = X(i, :) .* X(i + 1, :);
end

% Eliminar filas con suma cero y última fila redundante
X(sum(abs(X), 2) == 0, :) = [];
X(end, :) = [];

%% -------------------- PCA --------------------
S = cov(X');
[V, D] = eig(S);
evals = sort(diag(D), 'descend');
TotVar = sum(evals);

Var = 0;
for i = 1:length(evals)
    Var = Var + evals(i);
    if Var / TotVar > 0.95
        break
    end
end
n_Dim = i;
fprintf('Se necesitan %d componentes para explicar >95%% de la varianza.\n', n_Dim);

[~, s_idxs] = sort(diag(D), 'descend');
Vp = V(:, s_idxs(1:n_Dim));
Y = Vp' * (X - mean(X, 2));  % Y tiene dimensiones [n_Dim, N]
Y = Y';  % Para tener sujetos en filas

%% -------------------- MANOVA CON FITRM + MANOVA --------------------
% Convertir a tabla
T = array2table(Y);
for i = 1:n_Dim
    T.Properties.VariableNames{i} = ['PC' num2str(i)];
end
T.Grupo = categorical(MF);  % Agregamos la variable de grupo como categórica

% Ajustar modelo MANOVA multivariado
rm = fitrm(T, sprintf('PC1-PC%d ~ Grupo', n_Dim));
manovatbl = manova(rm);  % Resultado multivariante (como en tu imagen)

%% -------------------- RESULTADOS --------------------
disp('--- MANOVA multivariada (estilo tabla de resultados) ---');
disp(manovatbl);
