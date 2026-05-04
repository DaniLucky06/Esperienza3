% MASSA IN KILOGRAMMI
% PERIODI IN SECONDI

dati = readmatrix('dati_pendolo_semplice.csv');

%misure
masse = dati(:,1) ./1000;
periodi = dati(:,2);
media_periodo = 1.86120245020182;
periodo_teorico = 1.86016115911165;

%errore di risoluzione
errore_ris_masse = 0.1/1000;

%sigma msiure
sgm_masse = errore_ris_masse/sqrt(12);
sgm_periodo = dati(:,3);





%% calcolo fit lineare
x_err = sgm_masse;
y_err = sgm_periodo;
x = masse;
y = periodi;

% Varie somme
W = 1 ./ (y_err .^ 2); % "Pesi" - Weights
S_W   = sum(W);
S_XW  = sum(x .* W);
S_YW  = sum(y .* W);
S_XXW = sum((x .^ 2) .* W);
S_XYW = sum(x .* y .* W);

% Delta
D_W = S_W * S_XXW - S_XW^2;

% Fit iniziale
b = (1 / D_W) * (S_W * S_XYW - S_XW * S_YW);

% migliorare
y_err_i = sqrt(y_err .^ 2 + (b .* x_err) .^ 2); % errore 2
b_err = sqrt(S_W / D_W);
b1 = b + 2 * b_err;

while abs(b - b1) > b_err % loop di miglioramento (per k non parte nemmeno) y_err_i = sqrt(y_err .^ 2 + (b .* x_err) .^ 2);
    % Varie somme
    W = 1 ./ (y_err_i .^ 2); % "Pesi" - Weights
    S_W   = sum(W);
    S_XW  = sum(x .* W);
    S_YW  = sum(y .* W);
    S_XXW = sum((x .^ 2) .* W);
    S_XYW = sum(x .* y .* W);
    
    % Delta
    D_W = S_W * S_XXW - S_XW^2;

    b_err = sqrt(S_W / D_W);
    b1 = b;
    b = (1 / D_W) * (S_W * S_XYW - S_XW * S_YW);
end

% Intercetta
a = (1 / D_W) * (S_XXW * S_YW - S_XW * S_XYW);

% Errori su fit
a_err = sqrt(S_XXW / D_W);
b_err = sqrt(S_W / D_W); % errore finale


%% grafici
figure(1);
hold on;

errorbar(masse,periodi,sgm_periodo,sgm_periodo,sgm_masse,sgm_masse,  ...
       'o', ...
        ...
       'Color', 'k', 'LineWidth', 0.1, ...
       ...
       'MarkerFaceColor', 'k' , 'MarkerSize', 2, ...
       'HandleVisibility', 'off');

% valori vari per comparare: media periodo e periodo teorico

yline(media_periodo, '-' , ...
    num2str(media_periodo), 'Color', 'r', 'LabelHorizontalAlignment', 'left', ...
    'DisplayName','Periodo medio');

yline(periodo_teorico, '-' , ...
    num2str(periodo_teorico), 'Color', 'g', 'LabelHorizontalAlignment', 'left', ...
    'DisplayName','Periodo teorico');


%lunghezza asse x y
xlim([0, 1]);

%grafico fit lineare
x = linspace(min(masse)-0.2,max(masse)+0.2);
y = a+b*x;
plot(x,y,'DisplayName','Fit lineare')

legend('Location','best');
<<<<<<< HEAD
=======



>>>>>>> origin/Morla
