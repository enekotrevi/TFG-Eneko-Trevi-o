%% ========================================================================
% ANÁLISIS COMPLETO: COVARIANCE‐PCA → MANOVA (AgeGroup, BMI_Value, GroupName)
% Fichero: Five_times_Kinematics_HT_Query_SIN_OUTLIERS_con tres_variables.xlsx
% ========================================================================

%% 0) Carga de datos y definición de factores
% -------------------------------------------
Tdata = readtable('Five_times_Kinematics_HT_Query_SIN_OUTLIERS_con_tres_variables.xlsx');

% Factores
AgeGroup  = categorical(Tdata.AgeGroup);     % p.ej. '60-67','68-75'
BMI_Value = categorical(Tdata.BMI_Value);    % p.ej. 'High','Low'
GroupName = categorical(Tdata.GroupName);    % p.ej. 'Conservador','Artroplástia'

% Construcción de X_lin con todas las columnas cinemáticas
numericVars = Tdata.Properties.VariableNames(5:end);
X_lin = Tdata{:, numericVars}';   % → [#features × #observaciones]

%% 1) Generación del “quadratic kernel”
% --------------------------------------
X = X_lin;
lin_dim = size(X,1);

% a) Añade cuadrados de cada fila
for i = 1:lin_dim
    X(lin_dim + i, :) = X(i, :).^2;
end

% b) Añade interacciones solo entre puntos adyacentes
quad_dim = size(X,1);
for i = 1:lin_dim
    if mod(i,5)==0
        continue;  % saltar índices múltiplos de 5
    end
    X(quad_dim + i, :) = X(i, :) .* X(i+1, :);
end

% c) Elimina filas todas‐cero
rowsToKeep = find(any(X~=0, 2));
X = X(rowsToKeep, :);

% d) Elimina la última fila sobrante
X(end, :) = [];

%% 2) PCA sobre la matriz de covarianza
% -------------------------------------
% a) Matriz de covarianza
S = cov(X');

% b) Descomposición en vectores y valores propios
[V, D] = eig(S);

% c) Orden y varianza acumulada
eVals  = sort(diag(D), 'descend');
cumVar = cumsum(eVals) / sum(eVals);

% d) Selección de componentes (≥95%)
p = find(cumVar >= 0.95, 1, 'first');

% e) Vectores propios para PC’s
[~, idx] = sort(diag(D), 'descend');
Vp = V(:, idx(1:p));

% f) Cálculo de scores: Y = Vp' * (X – media)
n = size(X, 2);
Xc = X - repmat(mean(X,2), 1, n);
Y  = Vp' * Xc;    % → [p × n]

%% 3) Preparación de tabla para MANOVA
% ------------------------------------
Y_scores = Y';    % → [n × p]
varNames  = arrayfun(@(k) sprintf('PC%d',k), 1:p, 'Uni', false);

Tpc = array2table(Y_scores, 'VariableNames', varNames);
Tpc.Age       = AgeGroup;
Tpc.BMI       = BMI_Value;
Tpc.Treatment = GroupName;

%% 4) Ajuste MANOVA
% -----------------
fm = sprintf('PC1-PC%d ~ Age*BMI*Treatment', p);
rm = fitrm(Tpc, fm);

%% 5) Resultados MANOVA global
% ----------------------------
disp('=== MANOVA global (Wilks, Pillai, etc.) ===');
manovatbl = manova(rm);
disp(manovatbl);

%% 6) ANOVA univariados por PC
% ----------------------------
disp('=== ANOVA univariados (ranova) ===');
uniTbl = ranova(rm);
disp(uniTbl);

%% 7) Post-hoc: comparaciones múltiples
% -------------------------------------
groupNames = {'AgeGroup','BMI_Value','GroupName'};
factors   = {AgeGroup, BMI_Value, GroupName};

for k = 1:p
    fprintf('\n===== PC %d =====\n', k);
    dataPC = Y(k, :)';
    [~, ~, stats] = anovan(dataPC, factors, ...
        'model','interaction', ...
        'varnames', groupNames, ...
        'display','off');

    for f = 1:numel(groupNames)
        idx = f;  % 1→AgeGroup, 2→BMI_Value, 3→GroupName
        fprintf('\n-- Post-hoc %s --\n', groupNames{f});

        % Capturamos la salida y la convertimos en tabla
        [c, ~] = multcompare(stats, ...
            'Dimension', idx, ...
            'CType',     'bonferroni', ...
            'Display',   'off');        % evita la figura
        Tmc = array2table(c, ...
            'VariableNames', {'Group1','Group2','Lower','MeanDiff','Upper','pValue'});
        disp(Tmc);
    end
end


%% Opcional: tamaños de efecto univariados
%eta2 = etaSquared(rm);
%disp('=== Tamaño del efecto (eta²) ===');
%disp(eta2);
