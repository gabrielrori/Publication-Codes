clc
clear all
close all

%Opciones de selección
opt.dict = 'wave';      %Diccionario 'rand' 'wave'
opt.stepT = 'cte';      %Nombre de tamaño de paso 'cte' 'BB1' 'YnB'
opt.wtype = 'cont';     %Estrategia de Warm Start 'cont' 'line' 'cuad' 
opt.wname = 'db4';      %Tipo de Transformada Wavelet 'db1' 'db4' 'bior4.4'
opt.sample_img = 'baboon.gif'; % 'lena.bmp' 'barbara.bmp' 'goldhill.gif' 'house.gif' 'baboon.gif'

%Constantes lógicas
opt.steplim = 0;        %Activación de límites para el tamaño de paso
opt.warms = 0;          %Activación de Warm Start
opt.screening = 0;      %Activación de Screening

%Parámetros
L = 20;                 %número de iteraciones
opt.alpha = 0.05;       %tamaño de paso
opt.lambda = 80*10^-3;  %parámetro de regularización
opt.sigma = 1.5*10^-1;  %ruido aditivo
opt.ws_inc = 2;         %incremento para Warm Start
opt.warms_iter = 5;     %iteración hasta la cual se evalua Warm Start
opt.screeningfeatures = 1; %intensidad de Screening
opt.n = 3;              %nivel de evaluación de la Transformada Wavelets

%Generacion de datos y variables iniciales
[x_orig, b, x_init, var_init] = data_gen(opt);

%Casos

%ISTA
x = x_init;
var = var_init;
[stats_ista, x_ista, var_ista] = general_ista(opt, var, b, x, L);

%FISTA
x = x_init;
var = var_init;
[stats_fista, x_fista, var_fista] = general_fista(opt, var, b, x, L);

%FISTA TPAL (tamaño de paso automático con límites)
opt.stepT = 'BB1';
opt.steplim = 1;
opt.stepCH = 0.2;       %valor máximo de tamaño de paso de Cauchy
opt.stepCL = 0.1;       %valor mínimo de tamaño de paso de Caichy
opt.stepH = 0.3 ;       %valor máximo de tamaño de paso automático
opt.stepL = 0.2;        %valor mínimo de tamaño de paso automático
opt.step0_mult = 0.1;   %coeficiente que multiplica al t. de paso de Cauchy
opt.step_mult = 0.3;    %coeficiente que multiplica al t. de paso automático
x = x_init;
var = var_init;
[stats_fista_auto, x_fista_auto, var_fista_auto] = general_fista(opt, var, b, x, L);

%FISTA WS suave
opt.warms = 1;
opt.stepT = 'cte';
opt.ws_inc = 2;
x = x_init;
var = var_init;
[stats_fista_wsL, x_fista_wsL, var_fista_wsL] = general_fista(opt, var, b, x, L);

%FISTA WS fuerte
opt.ws_inc = 10;
x = x_init;
var = var_init;
[stats_fista_wsM, x_fista_wsM, var_fista_wsM] = general_fista(opt, var, b, x, L);

%FISTA WS excesivo
opt.ws_inc = 20;
x = x_init;
var = var_init;
[stats_fista_wsH, x_fista_wsH, var_fista_wsH] = general_fista(opt, var, b, x, L);
        
%FISTA SC suave
opt.warms = 0;
opt.screening = 1;
opt.screeningfeatures = 0.3;
x = x_init;
var = var_init;
[stats_fista_scL, x_fista_sL, var_fista_scL] = general_fista(opt, var, b, x, L);

%%FISTA SC medio
opt.screeningfeatures = 0.5;
x = x_init;
var = var_init;
[stats_fista_scM, x_fista_scM, var_fista_scM] = general_fista(opt, var, b, x, L);

%FISTA SC fuerte
opt.screeningfeatures = 1;
x = x_init;
var = var_init;
[stats_fista_scH, x_fista_scH, var_fista_scH] = general_fista(opt, var, b, x, L);

%FISTA OPT
opt.warms = 1;
opt.ws_inc = 1.5;
opt.screening = 1;
opt.screeningfeatures = 1;
opt.stepT = 'BB1';
opt.steplim = 1;
opt.step0_mult = 0.1;
opt.step_mult = 0.3;
opt.stepCH =   0.2;
opt.stepCL =   0.1;
opt.stepH =  0.3;
opt.stepL =   0.2;
x = x_init;
var = var_init;
[stats_fista_all, x_fista_all, var_fista_all] = general_fista(opt, var, b, x, L);

%Graficos

%Comparación entre FISTA Regular, TPAL, WSO, SCO y OPT
figure(1)
subplot(1,3,1), plot(stats_fista.fx,'LineWidth',3)
hold on
plot(stats_fista_auto.fx,'LineWidth',3)
plot(stats_fista_wsM.fx,'LineWidth',3)
plot(stats_fista_scH.fx,'LineWidth',3)
plot(stats_fista_all.fx,'LineWidth',3)
set(gca,'FontSize',24)
grid;
xlabel('Iteraciones'), ylabel('Funcion de error')
xlim([1 20])
subplot(1,3,2), plot(stats_fista.l1,'LineWidth',3)
hold on
plot(stats_fista_auto.l1,'LineWidth',3)
plot(stats_fista_wsM.l1,'LineWidth',3)
plot(stats_fista_scH.l1,'LineWidth',3)
plot(stats_fista_all.l1,'LineWidth',3)
set(gca,'FontSize',24)
grid;
xlabel('Iteraciones'), ylabel('Norma l1 de x')
xlim([1 20])
subplot(1,3,3), plot(stats_fista.F,'LineWidth',3)
hold on
plot(stats_fista_auto.F,'LineWidth',3)
plot(stats_fista_wsM.F,'LineWidth',3)
plot(stats_fista_scH.F,'LineWidth',3)
plot(stats_fista_all.F,'LineWidth',3)
set(gca,'FontSize',24)
legend('FISTA', 'FISTA TPAL','FISTA WSO',  'FISTA SCO', 'FISTA OPT',  'FontSize',  26)
grid;
xlabel('Iteraciones'), ylabel('Funcion de costo F(x)')
xlim([1 20])

%Comparación de tiempo para FISTA Regular, TPAL, WSO, SCO y OPT
figure(2)
plot(var_fista.time(1:L), stats_fista.F(1:L),'LineWidth',3)
hold on
plot(var_fista_auto.time(1:L), stats_fista_auto.F(1:L),'LineWidth',3)
plot(var_fista_wsM.time(1:L), stats_fista_wsM.F(1:L),'LineWidth',3)
plot(var_fista_scH.time(1:L), stats_fista_scH.F(1:L),'LineWidth',3)
plot(var_fista_all.time(1:L), stats_fista_all.F(1:L),'LineWidth',3)
set(gca,'FontSize',24)
legend('FISTA', 'FISTA TPAL','FISTA WSO',  'FISTA SCO', 'FISTA OPT',  'FontSize',  26)
grid;
xlabel('Tiempo (s)'), ylabel('Funcion de Costo F(x)')
xlim([0 2.7])

figure(3)
subplot(1,3,1), plot(stats_ista.fx,'LineWidth',3)
hold on
plot(stats_fista.fx,'LineWidth',3)
set(gca,'FontSize',24)
grid;
xlabel('Iteraciones'), ylabel('Funcion de error')
xlim([1 20])
subplot(1,3,2), plot(stats_ista.l1,'LineWidth',3)
hold on
plot(stats_fista.l1,'LineWidth',3)
set(gca,'FontSize',24)
grid;
xlabel('Iteraciones'), ylabel('Norma l1 de x')
xlim([1 20])
subplot(1,3,3), plot(stats_ista.F,'LineWidth',3)
hold on
plot(stats_fista.F,'LineWidth',3)
set(gca,'FontSize',24)
legend('ISTA','FISTA', 'FontSize', 26)
grid;
xlabel('Iteraciones'), ylabel('Funcion de costo F(x)')
xlim([1 20])

figure(4)
subplot(1,3,1), plot(stats_fista.fx,'LineWidth',3)
hold on
plot(stats_fista_auto.fx,'LineWidth',3)
set(gca,'FontSize',24)
grid;
xlabel('Iteraciones'), ylabel('Funcion de error')
% xlim([0 100])
xlim([1 20])
subplot(1,3,2), plot(stats_fista.l1,'LineWidth',3)
hold on
plot(stats_fista_auto.l1,'LineWidth',3)
set(gca,'FontSize',24)
grid;
xlabel('Iteraciones'), ylabel('Norma l1 de x')
% xlim([0 100])
xlim([1 20])
subplot(1,3,3), plot(stats_fista.F,'LineWidth',3)
hold on
plot(stats_fista_auto.F,'LineWidth',3)
set(gca,'FontSize',24)
legend('FISTA','FISTA TPAL', 'FontSize', 26)
grid;
xlabel('Iteraciones'), ylabel('Funcion de costo F(x)')
xlim([1 20])

figure(5)
subplot(1,3,1), plot(stats_fista.fx,'LineWidth',3)
hold on
plot(stats_fista_wsL.fx,'LineWidth',3)
plot(stats_fista_wsM.fx,'LineWidth',3)
plot(stats_fista_wsH.fx,'LineWidth',3)
set(gca,'FontSize',24)
grid;
xlabel('Iteraciones'), ylabel('Funcion de error')
xlim([1 20])
subplot(1,3,2), plot(stats_fista.l1,'LineWidth',3)
hold on
plot(stats_fista_wsL.l1,'LineWidth',3)
plot(stats_fista_wsM.l1,'LineWidth',3)
plot(stats_fista_wsH.l1,'LineWidth',3)
set(gca,'FontSize',24)
grid;
xlabel('Iteraciones'), ylabel('Norma l1 de x')
xlim([1 20])
subplot(1,3,3), plot(stats_fista.F,'LineWidth',3)
hold on
plot(stats_fista_wsL.F,'LineWidth',3)
plot(stats_fista_wsM.F,'LineWidth',3)
plot(stats_fista_wsH.F,'LineWidth',3)
set(gca,'FontSize',24)
legend('FISTA', 'FISTA WS suave','FISTA WS fuerte',  'FISTA WS excesivo',  'FontSize',  26)
grid;
xlabel('Iteraciones'), ylabel('Funcion de costo F(x)')
xlim([1 20])

figure(6)
subplot(1,3,1), plot(stats_fista.fx,'LineWidth',3)
hold on
plot(stats_fista_wsM.fx,'LineWidth',3)
set(gca,'FontSize',24)
grid;
xlabel('Iteraciones'), ylabel('Funcion de error')
xlim([1 20])
subplot(1,3,2), plot(stats_fista.l1,'LineWidth',3)
hold on
plot(stats_fista_wsM.l1,'LineWidth',3)
set(gca,'FontSize',24)
grid;
xlabel('Iteraciones'), ylabel('Norma l1 de x')
xlim([1 20])
subplot(1,3,3), plot(stats_fista.F,'LineWidth',3)
hold on
plot(stats_fista_wsM.F,'LineWidth',3)
set(gca,'FontSize',24)
legend('FISTA', 'FISTA WSO', 'FontSize',  26)
grid;
xlabel('Iteraciones'), ylabel('Funcion de costo F(x)')
xlim([1 20])

figure(7)
subplot(1,3,1), plot(stats_fista.fx,'LineWidth',3)
hold on
plot(stats_fista_scL.fx,'LineWidth',3)
plot(stats_fista_scM.fx,'LineWidth',3)
plot(stats_fista_scH.fx,'LineWidth',3)
set(gca,'FontSize',24)
grid;
xlabel('Iteraciones'), ylabel('Funcion de error')
xlim([1 20])
subplot(1,3,2), plot(stats_fista.l1,'LineWidth',3)
hold on
plot(stats_fista_scL.l1,'LineWidth',3)
plot(stats_fista_scM.l1,'LineWidth',3)
plot(stats_fista_scH.l1,'LineWidth',3)
set(gca,'FontSize',24)
grid;
xlabel('Iteraciones'), ylabel('Norma l1 de x')
xlim([1 20])
subplot(1,3,3), plot(stats_fista.F,'LineWidth',3)
hold on
plot(stats_fista_scL.F,'LineWidth',3)
plot(stats_fista_scM.F,'LineWidth',3)
plot(stats_fista_scH.F,'LineWidth',3)
set(gca,'FontSize',24)
legend('FISTA', 'FISTA SC suave','FISTA SC medio',  'FISTA SC fuerte',  'FontSize',  26)
grid;
xlabel('Iteraciones'), ylabel('Funcion de costo F(x)')
xlim([1 20])

figure(8)
subplot(1,3,1), plot(stats_fista.fx,'LineWidth',3)
hold on
plot(stats_fista_scH.fx,'LineWidth',3)
set(gca,'FontSize',24)
grid;
xlabel('Iteraciones'), ylabel('Funcion de error')
xlim([1 20])
subplot(1,3,2), plot(stats_fista.l1,'LineWidth',3)
hold on
plot(stats_fista_scH.l1,'LineWidth',3)
set(gca,'FontSize',24)
grid;
xlabel('Iteraciones'), ylabel('Norma l1 de x')
xlim([1 20])
subplot(1,3,3), plot(stats_fista.F,'LineWidth',3)
hold on
plot(stats_fista_scH.F,'LineWidth',3)
set(gca,'FontSize',24)
legend('FISTA', 'FISTA SCO',  'FontSize',  26)
grid;
xlabel('Iteraciones'), ylabel('Funcion de costo F(x)')
xlim([1 20])