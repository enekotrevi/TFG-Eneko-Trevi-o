clear; close all; clc;

%% --- Cargar datos ---
data = xlsread('Five_times_Kinematics_HT_Query_Quadratic_SIN_OUTLIERS.xlsx');
Age     = data(:,2);
BMI     = data(:,3);
Therapy = data(:,4);
X_lin   = data(:,5:end)';

%% --- Reagrupar factores categóricos (0 y 1) ---
AgeGroup = zeros(size(Age));
AgeGroup(Age >= 68) = 1;

BMIGroup = zeros(size(BMI));
BMIGroup(BMI >= 30) = 1;

TherapyGroup = zeros(size(Therapy));
TherapyGroup(Therapy ~= 1) = 1;

%% --- Kernel cuadrático + PCA (usando tu código original) ----------------
% Generación del kernel cuadrático con interacciones adyacentes
X = X_lin;
lin_dim = size(X,1);
for i = 1:lin_dim
    X(lin_dim+i,:) = X(i,:).^2;
end

quad_dim = size(X,1);
for i = 1:lin_dim
    if mod(i,5) == 0
        i = i + 1;
    end
    if i < lin_dim
        X(quad_dim+i,:) = X(i,:) .* X(i+1,:);
    end
end

% Eliminar filas con solo ceros
rowsToKeep = [];
for i = 1:size(X, 1)
    if sum(abs(X(i, :))) ~= 0
        rowsToKeep = [rowsToKeep; i];
    end
end
X = X(rowsToKeep, :);

% Eliminar última fila irrelevante
X(end,:) = [];

% PCA baseline
S = cov(X');
[V, D] = eig(S);
n = size(X, 2);
dim = size(X, 1);

sort_EVal = sort(diag(D), 'descend');
TotVar = sum(sort_EVal);
Var = 0;
for i = 1:dim
    Var = Var + sort_EVal(i);
    Var99 = Var / TotVar;
    if Var99 > 0.95
        break
    end
end
n_Dimensions = i;

% Proyección sobre PCs significativas
[~, s_idxs] = sort(diag(D), 'descend');
Vp = V(:, s_idxs(1:n_Dimensions));
Y_all = Vp' * (X - repmat(mean(X')', [1 n]));

% Separar primeras 2 PCs (significativas por permutación) y PCs 3–6
Y_signif    = Y_all(1:2,:)';   % 84% varianza
Y_nonsignif = Y_all(3:6,:)';   % 11% varianza

%% --- Función auxiliar segura para MANOVA1 ---
function p = safe_manova1(Y, group_vector, factor_name)
    [groups, ~, idx] = unique(group_vector);
    valid_idx = true(size(group_vector));
    for i = 1:length(groups)
        if sum(idx == i) < 2
            fprintf('⚠️  Eliminado grupo "%s" del factor "%s" por tamaño insuficiente.\n', ...
                num2str(groups(i)), factor_name);
            valid_idx(idx == i) = false;
        end
    end
    Y_clean = Y(valid_idx,:);
    g_clean = group_vector(valid_idx);
    if numel(unique(g_clean)) < 2
        warning('No hay suficientes grupos válidos para el factor %s.', factor_name);
        p = NaN;
    else
        p = manova1(Y_clean, grp2idx(g_clean));
    end
end

%% --- MANOVA 1: PCs significativas (84%) ---
fprintf('\n--- MANOVA 1: PCs Significativas (84%% varianza) ---\n');
p_age1     = safe_manova1(Y_signif, AgeGroup, 'AgeGroup');
p_bmi1     = safe_manova1(Y_signif, BMIGroup, 'BMIGroup');
p_therapy1 = safe_manova1(Y_signif, TherapyGroup, 'TherapyGroup');

results1 = table([p_age1; p_bmi1; p_therapy1], ...
    'VariableNames', {'p_value'}, ...
    'RowNames', {'AgeGroup','BMIGroup','TherapyGroup'});
disp(results1);

%% --- MANOVA 2: PCs no significativas (11%) ---
fprintf('\n--- MANOVA 2: PCs No Significativas (11%% varianza) ---\n');
p_age2     = safe_manova1(Y_nonsignif, AgeGroup, 'AgeGroup');
p_bmi2     = safe_manova1(Y_nonsignif, BMIGroup, 'BMIGroup');
p_therapy2 = safe_manova1(Y_nonsignif, TherapyGroup, 'TherapyGroup');

results2 = table([p_age2; p_bmi2; p_therapy2], ...
    'VariableNames', {'p_value'}, ...
    'RowNames', {'AgeGroup','BMIGroup','TherapyGroup'});
disp(results2);
