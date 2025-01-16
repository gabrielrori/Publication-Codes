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



T = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BUCLE

for i = 1:T
    i
    
    
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
    opt.ws_inc = 10;
    x = x_init;
    var = var_init;
    [stats_fista_wsM, x_fista_wsM, var_fista_wsM] = general_fista(opt, var, b, x, L);

    %FISTA SC suave
    opt.warms = 0;
    opt.screening = 1;
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

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stats = stats_ista.F; %23  0.0777
    time = var_ista.time;
    conv = zeros(1,L);
    convergencia = 0;

    for(j=1:L)
        dif(j) = stats(j)-stats(j+1);
        if(abs(dif(j))<100) && (j> 1)
    %     if (stats(i) < 900)
            convergencia=1;
            conv(j)=1;
    %        break
        else
            conv(j)=0;
        end
    end

    conv_iter_ista(i) =  L;
    conv_time_ista(i) = time(conv_iter_ista(i)+1);


    %%%%%%%%%%%%%%%%%%%%%%%%%%
    stats = stats_fista.F; %23  0.0777
    time = var_fista.time;
    conv = zeros(1,L);
    convergencia = 0;
    %0.02

    for(j=1:L)
        dif(j) = stats(j)-stats(j+1);
        if(abs(dif(j))<100) && (j> 1)
    %     if (stats(i) < 900)
            convergencia=1;
            conv(j)=1;
    %        break
        else
            conv(j)=0;
        end
    end

    conv_iter_fista(i) = find(conv,1,'first');
    conv_time_fista(i) = time(conv_iter_fista(i)+1);



    %%%%%%%%%%%%%%%%%%%%%%%%%


    stats = stats_fista_auto.F; %16 0.0461
    time = var_fista_auto.time;
    conv = zeros(1,L);
    convergencia = 0;
    %0.02

    for(j=1:L)
        dif(j) = stats(j)-stats(j+1);
        if(abs(dif(j))<100) && (j> 1)
    %     if (stats(i) < 900)
            convergencia=1;
            conv(j)=1;
    %        break
        else
            conv(j)=0;
        end
    end

    conv_iter_fista_auto(i) =  find(conv,1,'first');
    conv_time_fista_auto(i) = time(conv_iter_fista_auto(i)+1);



    %%%%%%%%%%%%%%%%%%%%%%%%%%

    stats = stats_fista_wsM.F; %19  0.0499
    time = var_fista_wsM.time;
    conv = zeros(1,L);
    convergencia = 0;
    %0.02

    for(j=1:L)
        dif(j) = stats(j)-stats(j+1);
        if(abs(dif(j))<100) && (j> 1)
    %     if (stats(i) < 900)
            convergencia=1;
            conv(j)=1;
    %        break
        else
            conv(j)=0;
        end
    end

    conv_iter_fista_ws(i) = find(conv,1,'first');
    conv_time_fista_ws(i) = time(conv_iter_fista_ws(i)+1);



    %%%%%%%%%%%%%%%%%%%%%%%%%%



    stats = stats_fista_scH.F; %15  0.0347
    time = var_fista_scH.time;
    conv = zeros(1,L);

    convergencia = 0;
    %0.02

    for(j=1:L)
        dif(j) = stats(j)-stats(j+1);
        if(abs(dif(j))<100) && (j> 1)
    %     if (stats(i) < 900)
            convergencia=1;
            conv(j)=1;
    %        break
        else
            conv(j)=0;
        end
    end

    conv_iter_fista_sc(i) =  find(conv,1,'first');
    conv_time_fista_sc(i) = time(conv_iter_fista_sc(i)+1);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    stats = stats_fista_all.F; %11 0.0301
    time = var_fista_all.time;


    conv = zeros(1,L);

    convergencia = 0;
    %0.02

    for(j=1:L)
        dif(j) = stats(j)-stats(j+1);
        if(abs(dif(j))<100) && (j> 1)
    %     if (stats(i) < 900)
            convergencia=1;
            conv(j)=1;
    %        break
        else
            conv(j)=0;
        end
    end

    conv_iter_fista_all(i) = find(conv,1,'first');
    conv_time_fista_all(i) = time(conv_iter_fista_all(i)+1);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    stats_ista_fx_1(i) = stats_ista.fx(1);
    stats_ista_l1_1(i) = stats_ista.l1(1);
    stats_ista_F_1(i) = stats_ista.F(1);

    
    stats_ista_fx_5(i) = stats_ista.fx(5);
    stats_ista_l1_5(i) = stats_ista.l1(5);
    stats_ista_F_5(i) = stats_ista.F(5);
    
    
    stats_ista_fx_10(i) = stats_ista.fx(10);
    stats_ista_l1_10(i) = stats_ista.l1(10);
    stats_ista_F_10(i) = stats_ista.F(10);
    
    
    stats_ista_fx_15(i) = stats_ista.fx(15);
    stats_ista_l1_15(i) = stats_ista.l1(15);
    stats_ista_F_15(i) = stats_ista.F(15);


    stats_ista_fx_20(i) = stats_ista.fx(20);
    stats_ista_l1_20(i) = stats_ista.l1(20);
    stats_ista_F_20(i) = stats_ista.F(20);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%


    stats_fista_fx_1(i) = stats_fista.fx(1);
    stats_fista_l1_1(i) = stats_fista.l1(1);
    stats_fista_F_1(i) = stats_fista.F(1);
    
    stats_fista_fx_5(i) = stats_fista.fx(5);
    stats_fista_l1_5(i) = stats_fista.l1(5);
    stats_fista_F_5(i) = stats_fista.F(5);

    stats_fista_fx_10(i) = stats_fista.fx(10);
    stats_fista_l1_10(i) = stats_fista.l1(10);
    stats_fista_F_10(i) = stats_fista.F(10);

    
    stats_fista_fx_15(i) = stats_fista.fx(15);
    stats_fista_l1_15(i) = stats_fista.l1(15);
    stats_fista_F_15(i) = stats_fista.F(15);

    stats_fista_fx_20(i) = stats_fista.fx(20);
    stats_fista_l1_20(i) = stats_fista.l1(20);
    stats_fista_F_20(i) = stats_fista.F(20);


    %%%%%%%%%%%%%%%%%%%%%%%%%%


    stats_fista_auto_fx_1(i) = stats_fista_auto.fx(1);
    stats_fista_auto_l1_1(i) = stats_fista_auto.l1(1);
    stats_fista_auto_F_1(i) = stats_fista_auto.F(1);
    
    stats_fista_auto_fx_5(i) = stats_fista_auto.fx(5);
    stats_fista_auto_l1_5(i) = stats_fista_auto.l1(5);
    stats_fista_auto_F_5(i) = stats_fista_auto.F(5);


    stats_fista_auto_fx_10(i) = stats_fista_auto.fx(10);
    stats_fista_auto_l1_10(i) = stats_fista_auto.l1(10);
    stats_fista_auto_F_10(i) = stats_fista_auto.F(10);

    
    stats_fista_auto_fx_15(i) = stats_fista_auto.fx(15);
    stats_fista_auto_l1_15(i) = stats_fista_auto.l1(15);
    stats_fista_auto_F_15(i) = stats_fista_auto.F(15);

    stats_fista_auto_fx_20(i) = stats_fista_auto.fx(20);
    stats_fista_auto_l1_20(i) = stats_fista_auto.l1(20);
    stats_fista_auto_F_20(i) = stats_fista_auto.F(20);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    stats_fista_ws_fx_1(i) = stats_fista_wsM.fx(1);
    stats_fista_ws_l1_1(i) = stats_fista_wsM.l1(1);
    stats_fista_ws_F_1(i) = stats_fista_wsM.F(1);
    
    stats_fista_ws_fx_5(i) = stats_fista_wsM.fx(5);
    stats_fista_ws_l1_5(i) = stats_fista_wsM.l1(5);
    stats_fista_ws_F_5(i) = stats_fista_wsM.F(5);

    stats_fista_ws_fx_10(i) = stats_fista.fx(10);
    stats_fista_ws_l1_10(i) = stats_fista_wsM.l1(10);
    stats_fista_ws_F_10(i) = stats_fista_wsM.F(10);
    
    stats_fista_ws_fx_15(i) = stats_fista_wsM.fx(15);
    stats_fista_ws_l1_15(i) = stats_fista_wsM.l1(15);
    stats_fista_ws_F_15(i) = stats_fista_wsM.F(15);


    stats_fista_ws_fx_20(i) = stats_fista_wsM.fx(20);
    stats_fista_ws_l1_20(i) = stats_fista_wsM.l1(20);
    stats_fista_ws_F_20(i) = stats_fista_wsM.F(20);
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    stats_fista_sc_fx_1(i) = stats_fista_scH.fx(1);
    stats_fista_sc_l1_1(i) = stats_fista_scH.l1(1);
    stats_fista_sc_F_1(i) = stats_fista_scH.F(1);
    
    
    stats_fista_sc_fx_5(i) = stats_fista_scH.fx(5);
    stats_fista_sc_l1_5(i) = stats_fista_scH.l1(5);
    stats_fista_sc_F_5(i) = stats_fista_scH.F(5);


    stats_fista_sc_fx_10(i) = stats_fista_scH.fx(10);
    stats_fista_sc_l1_10(i) = stats_fista_scH.l1(10);
    stats_fista_sc_F_10(i) = stats_fista_scH.F(10);
    
    stats_fista_sc_fx_15(i) = stats_fista_scH.fx(15);
    stats_fista_sc_l1_15(i) = stats_fista_scH.l1(15);
    stats_fista_sc_F_15(i) = stats_fista_scH.F(15);


    stats_fista_sc_fx_20(i) = stats_fista_scH.fx(20);
    stats_fista_sc_l1_20(i) = stats_fista_scH.l1(20);
    stats_fista_sc_F_20(i) = stats_fista_scH.F(20);


    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    stats_fista_all_fx_1(i) = stats_fista_all.fx(1);
    stats_fista_all_l1_1(i) = stats_fista_all.l1(1);
    stats_fista_all_F_1(i) = stats_fista_all.F(1);
    
    
    stats_fista_all_fx_5(i) = stats_fista_all.fx(5);
    stats_fista_all_l1_5(i) = stats_fista_all.l1(5);
    stats_fista_all_F_5(i) = stats_fista_all.F(5);

    stats_fista_all_fx_10(i) = stats_fista_all.fx(10);
    stats_fista_all_l1_10(i) = stats_fista_all.l1(10);
    stats_fista_all_F_10(i) = stats_fista_all.F(10);
    
    stats_fista_all_fx_15(i) = stats_fista_all.fx(15);
    stats_fista_all_l1_15(i) = stats_fista_all.l1(15);
    stats_fista_all_F_15(i) = stats_fista_all.F(15);


    stats_fista_all_fx_20(i) = stats_fista_all.fx(20);
    stats_fista_all_l1_20(i) = stats_fista_all.l1(20);
    stats_fista_all_F_20(i) = stats_fista_all.F(20);

 

    
    
    
    x_all = waverec2(x_fista_all.val,x_fista_all.size,'db4');
    x_scH = waverec2(x_fista_scH.val,x_fista_scH.size,'db4');
    x_ws = waverec2(x_fista_wsM.val,x_fista_wsM.size,'db4');
    x_f = waverec2(x_fista.val,x_fista.size,'db4');
    x_i = waverec2(x_ista.val,x_ista.size,'db4');

    x_auto = waverec2(x_fista_auto.val,x_fista_auto.size,'db4');


    snr_all(i) = psnr(x_all,x_orig);
    snr_scH(i) = psnr(x_scH, x_orig);
    snr_ws(i) = psnr(x_ws,x_orig);
    snr_f(i) = psnr(x_f, x_orig);
    snr_i(i) = psnr(x_i, x_orig);
    snr_auto(i) = psnr(x_auto,x_orig);
    
    ssim_all(i) = ssim(x_all,x_orig);
    ssim_scH(i) = ssim(x_scH, x_orig);
    ssim_ws(i) = ssim(x_ws,x_orig);
    ssim_f(i) = ssim(x_f, x_orig);
    ssim_i(i) = ssim(x_i, x_orig);
    ssim_auto(i) = ssim(x_auto,x_orig);
    
    
end






%FUERA DEL BUCLE


%CONVERGENCIA Y TIEMPO DE CONVERGENCIA

conv_iter_mean_ista =  mean(conv_iter_ista)
conv_iter_std_ista =  std(conv_iter_ista)
conv_time_mean_ista = mean(conv_time_ista)
conv_time_std_ista = std(conv_time_ista)

conv_iter_mean_fista =  mean(conv_iter_fista)
conv_iter_std_fista =  std(conv_iter_fista)
conv_time_mean_fista = mean(conv_time_fista)
conv_time_std_fista = std(conv_time_fista)

conv_iter_mean_fista_auto =  mean(conv_iter_fista_auto)
conv_iter_std_fista_auto =  std(conv_iter_fista_auto)
conv_time_mean_fista_auto = mean(conv_time_fista_auto)
conv_time_std_fista_auto = std(conv_time_fista_auto)


conv_iter_mean_fista_ws =  mean(conv_iter_fista_ws)
conv_iter_std_fista_ws =  std(conv_iter_fista_ws)
conv_time_mean_fista_ws = mean(conv_time_fista_ws)
conv_time_std_fista_ws = std(conv_time_fista_ws)


conv_iter_mean_fista_sc =  mean(conv_iter_fista_sc)
conv_iter_std_fista_sc =  std(conv_iter_fista_sc)
conv_time_mean_fista_sc = mean(conv_time_fista_sc)
conv_time_std_fista_sc = std(conv_time_fista_sc)


conv_iter_mean_fista_all =  mean(conv_iter_fista_all)
conv_iter_std_fista_all =  std(conv_iter_fista_all)
conv_time_mean_fista_all = mean(conv_time_fista_all)
conv_time_std_fista_all = std(conv_time_fista_all)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ISTA

stats_ista_fx_1_mean =  mean(stats_ista_fx_1)
stats_ista_fx_1_std =  std(stats_ista_fx_1)
stats_ista_l1_1_mean = mean(stats_ista_l1_1)
stats_ista_l1_1_std = std(stats_ista_l1_1)
stats_ista_F_1_mean = mean(stats_ista_F_1)
stats_ista_F_1_std = std(stats_ista_F_1)

stats_ista_fx_5_mean =  mean(stats_ista_fx_5)
stats_ista_fx_5_std =  std(stats_ista_fx_5)
stats_ista_l1_5_mean = mean(stats_ista_l1_5)
stats_ista_l1_5_std = std(stats_ista_l1_5)
stats_ista_F_5_mean = mean(stats_ista_F_5)
stats_ista_F_5_std = std(stats_ista_F_5)

stats_ista_fx_10_mean =  mean(stats_ista_fx_10)
stats_ista_fx_10_std =  std(stats_ista_fx_10)
stats_ista_l1_10_mean = mean(stats_ista_l1_10)
stats_ista_l1_10_std = std(stats_ista_l1_10)
stats_ista_F_10_mean = mean(stats_ista_F_10)
stats_ista_F_10_std = std(stats_ista_F_10)

stats_ista_fx_15_mean =  mean(stats_ista_fx_15)
stats_ista_fx_15_std =  std(stats_ista_fx_15)
stats_ista_l1_15_mean = mean(stats_ista_l1_15)
stats_ista_l1_15_std = std(stats_ista_l1_15)
stats_ista_F_15_mean = mean(stats_ista_F_15)
stats_ista_F_15_std = std(stats_ista_F_15)

stats_ista_fx_20_mean =  mean(stats_ista_fx_20)
stats_ista_fx_20_std =  std(stats_ista_fx_20)
stats_ista_l1_20_mean = mean(stats_ista_l1_20)
stats_ista_l1_20_std = std(stats_ista_l1_20)
stats_ista_F_20_mean = mean(stats_ista_F_20)
stats_ista_F_20_std = std(stats_ista_F_20)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FISTA

stats_fista_fx_1_mean =  mean(stats_fista_fx_1)
stats_fista_fx_1_std =  std(stats_fista_fx_1)
stats_fista_l1_1_mean = mean(stats_fista_l1_1)
stats_fista_l1_1_std = std(stats_fista_l1_1)
stats_fista_F_1_mean = mean(stats_fista_F_1)
stats_fista_F_1_std = std(stats_fista_F_1)


stats_fista_fx_5_mean =  mean(stats_fista_fx_5)
stats_fista_fx_5_std =  std(stats_fista_fx_5)
stats_fista_l1_5_mean = mean(stats_fista_l1_5)
stats_fista_l1_5_std = std(stats_fista_l1_5)
stats_fista_F_5_mean = mean(stats_fista_F_5)
stats_fista_F_5_std = std(stats_fista_F_5)

stats_fista_fx_10_mean =  mean(stats_fista_fx_10)
stats_fista_fx_10_std =  std(stats_fista_fx_10)
stats_fista_l1_10_mean = mean(stats_fista_l1_10)
stats_fista_l1_10_std = std(stats_fista_l1_10)
stats_fista_F_10_mean = mean(stats_fista_F_10)
stats_fista_F_10_std = std(stats_fista_F_10)

stats_fista_fx_15_mean =  mean(stats_fista_fx_15)
stats_fista_fx_15_std =  std(stats_fista_fx_15)
stats_fista_l1_15_mean = mean(stats_fista_l1_15)
stats_fista_l1_15_std = std(stats_fista_l1_15)
stats_fista_F_15_mean = mean(stats_fista_F_15)
stats_fista_F_15_std = std(stats_fista_F_15)

stats_fista_fx_20_mean =  mean(stats_fista_fx_20)
stats_fista_fx_20_std =  std(stats_fista_fx_20)
stats_fista_l1_20_mean = mean(stats_fista_l1_20)
stats_fista_l1_20_std = std(stats_fista_l1_20)
stats_fista_F_20_mean = mean(stats_fista_F_20)
stats_fista_F_20_std = std(stats_fista_F_20)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FISTA AUTO

stats_fista_auto_fx_1_mean =  mean(stats_fista_auto_fx_1)
stats_fista_auto_fx_1_std =  std(stats_fista_auto_fx_1)
stats_fista_auto_l1_1_mean = mean(stats_fista_auto_l1_1)
stats_fista_auto_l1_1_std = std(stats_fista_auto_l1_1)
stats_fista_auto_F_1_mean = mean(stats_fista_auto_F_1)
stats_fista_auto_F_1_std = std(stats_fista_auto_F_1)

stats_fista_auto_fx_5_mean =  mean(stats_fista_auto_fx_5)
stats_fista_auto_fx_5_std =  std(stats_fista_auto_fx_5)
stats_fista_auto_l1_5_mean = mean(stats_fista_auto_l1_5)
stats_fista_auto_l1_5_std = std(stats_fista_auto_l1_5)
stats_fista_auto_F_5_mean = mean(stats_fista_auto_F_5)
stats_fista_auto_F_5_std = std(stats_fista_auto_F_5)

stats_fista_auto_fx_10_mean =  mean(stats_fista_auto_fx_10)
stats_fista_auto_fx_10_std =  std(stats_fista_auto_fx_10)
stats_fista_auto_l1_10_mean = mean(stats_fista_auto_l1_10)
stats_fista_auto_l1_10_std = std(stats_fista_auto_l1_10)
stats_fista_auto_F_10_mean = mean(stats_fista_auto_F_10)
stats_fista_auto_F_10_std = std(stats_fista_auto_F_10)


stats_fista_auto_fx_15_mean =  mean(stats_fista_auto_fx_15)
stats_fista_auto_fx_15_std =  std(stats_fista_auto_fx_15)
stats_fista_auto_l1_15_mean = mean(stats_fista_auto_l1_15)
stats_fista_auto_l1_15_std = std(stats_fista_auto_l1_15)
stats_fista_auto_F_15_mean = mean(stats_fista_auto_F_15)
stats_fista_auto_F_15_std = std(stats_fista_auto_F_15)


stats_fista_auto_fx_20_mean =  mean(stats_fista_auto_fx_20)
stats_fista_auto_fx_20_std =  std(stats_fista_auto_fx_20)
stats_fista_auto_l1_20_mean = mean(stats_fista_auto_l1_20)
stats_fista_auto_l1_20_std = std(stats_fista_auto_l1_20)
stats_fista_auto_F_20_mean = mean(stats_fista_auto_F_20)
stats_fista_auto_F_20_std = std(stats_fista_auto_F_20)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FISTA WS


stats_fista_ws_fx_1_mean =  mean(stats_fista_ws_fx_1)
stats_fista_ws_fx_1_std =  std(stats_fista_ws_fx_1)
stats_fista_ws_l1_1_mean = mean(stats_fista_ws_l1_1)
stats_fista_ws_l1_1_std = std(stats_fista_ws_l1_1)
stats_fista_ws_F_1_mean = mean(stats_fista_ws_F_1)
stats_fista_ws_F_1_std = std(stats_fista_ws_F_1)


stats_fista_ws_fx_5_mean =  mean(stats_fista_ws_fx_5)
stats_fista_ws_fx_5_std =  std(stats_fista_ws_fx_5)
stats_fista_ws_l1_5_mean = mean(stats_fista_ws_l1_5)
stats_fista_ws_l1_5_std = std(stats_fista_ws_l1_5)
stats_fista_ws_F_5_mean = mean(stats_fista_ws_F_5)
stats_fista_ws_F_5_std = std(stats_fista_ws_F_5)


stats_fista_ws_fx_10_mean =  mean(stats_fista_ws_fx_10)
stats_fista_ws_fx_10_std =  std(stats_fista_ws_fx_10)
stats_fista_ws_l1_10_mean = mean(stats_fista_ws_l1_10)
stats_fista_ws_l1_10_std = std(stats_fista_ws_l1_10)
stats_fista_ws_F_10_mean = mean(stats_fista_ws_F_10)
stats_fista_ws_F_10_std = std(stats_fista_ws_F_10)

stats_fista_ws_fx_15_mean =  mean(stats_fista_ws_fx_15)
stats_fista_ws_fx_15_std =  std(stats_fista_ws_fx_15)
stats_fista_ws_l1_15_mean = mean(stats_fista_ws_l1_15)
stats_fista_ws_l1_15_std = std(stats_fista_ws_l1_15)
stats_fista_ws_F_15_mean = mean(stats_fista_ws_F_15)
stats_fista_ws_F_15_std = std(stats_fista_ws_F_15)


stats_fista_ws_fx_20_mean =  mean(stats_fista_ws_fx_20)
stats_fista_ws_fx_20_std =  std(stats_fista_ws_fx_20)
stats_fista_ws_l1_20_mean = mean(stats_fista_ws_l1_20)
stats_fista_ws_l1_20_std = std(stats_fista_ws_l1_20)
stats_fista_ws_F_20_mean = mean(stats_fista_ws_F_20)
stats_fista_ws_F_20_std = std(stats_fista_ws_F_20)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FISTA SC

stats_fista_sc_fx_1_mean =  mean(stats_fista_sc_fx_1)
stats_fista_sc_fx_1_std =  std(stats_fista_sc_fx_1)
stats_fista_sc_l1_1_mean = mean(stats_fista_sc_l1_1)
stats_fista_sc_l1_1_std = std(stats_fista_sc_l1_1)
stats_fista_sc_F_1_mean = mean(stats_fista_sc_F_1)
stats_fista_sc_F_1_std = std(stats_fista_sc_F_1)

stats_fista_sc_fx_5_mean =  mean(stats_fista_sc_fx_5)
stats_fista_sc_fx_5_std =  std(stats_fista_sc_fx_5)
stats_fista_sc_l1_5_mean = mean(stats_fista_sc_l1_5)
stats_fista_sc_l1_5_std = std(stats_fista_sc_l1_5)
stats_fista_sc_F_5_mean = mean(stats_fista_sc_F_5)
stats_fista_sc_F_5_std = std(stats_fista_sc_F_5)

stats_fista_sc_fx_10_mean =  mean(stats_fista_sc_fx_10)
stats_fista_sc_fx_10_std =  std(stats_fista_sc_fx_10)
stats_fista_sc_l1_10_mean = mean(stats_fista_sc_l1_10)
stats_fista_sc_l1_10_std = std(stats_fista_sc_l1_10)
stats_fista_sc_F_10_mean = mean(stats_fista_sc_F_10)
stats_fista_sc_F_10_std = std(stats_fista_sc_F_10)


stats_fista_sc_fx_15_mean =  mean(stats_fista_sc_fx_15)
stats_fista_sc_fx_15_std =  std(stats_fista_sc_fx_15)
stats_fista_sc_l1_15_mean = mean(stats_fista_sc_l1_15)
stats_fista_sc_l1_15_std = std(stats_fista_sc_l1_15)
stats_fista_sc_F_15_mean = mean(stats_fista_sc_F_15)
stats_fista_sc_F_15_std = std(stats_fista_sc_F_15)


stats_fista_sc_fx_20_mean =  mean(stats_fista_sc_fx_20)
stats_fista_sc_fx_20_std =  std(stats_fista_sc_fx_20)
stats_fista_sc_l1_20_mean = mean(stats_fista_sc_l1_20)
stats_fista_sc_l1_20_std = std(stats_fista_sc_l1_20)
stats_fista_sc_F_20_mean = mean(stats_fista_sc_F_20)
stats_fista_sc_F_20_std = std(stats_fista_sc_F_20)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FISTA ALL

stats_fista_all_fx_1_mean =  mean(stats_fista_all_fx_1)
stats_fista_all_fx_1_std =  std(stats_fista_all_fx_1)
stats_fista_all_l1_1_mean = mean(stats_fista_all_l1_1)
stats_fista_all_l1_1_std = std(stats_fista_all_l1_1)
stats_fista_all_F_1_mean = mean(stats_fista_all_F_1)
stats_fista_all_F_1_std = std(stats_fista_all_F_1)


stats_fista_all_fx_5_mean =  mean(stats_fista_all_fx_5)
stats_fista_all_fx_5_std =  std(stats_fista_all_fx_5)
stats_fista_all_l1_5_mean = mean(stats_fista_all_l1_5)
stats_fista_all_l1_5_std = std(stats_fista_all_l1_5)
stats_fista_all_F_5_mean = mean(stats_fista_all_F_5)
stats_fista_all_F_5_std = std(stats_fista_all_F_5)



stats_fista_all_fx_10_mean =  mean(stats_fista_all_fx_10)
stats_fista_all_fx_10_std =  std(stats_fista_all_fx_10)
stats_fista_all_l1_10_mean = mean(stats_fista_all_l1_10)
stats_fista_all_l1_10_std = std(stats_fista_all_l1_10)
stats_fista_all_F_10_mean = mean(stats_fista_all_F_10)
stats_fista_all_F_10_std = std(stats_fista_all_F_10)


stats_fista_all_fx_15_mean =  mean(stats_fista_all_fx_15)
stats_fista_all_fx_15_std =  std(stats_fista_all_fx_15)
stats_fista_all_l1_15_mean = mean(stats_fista_all_l1_15)
stats_fista_all_l1_15_std = std(stats_fista_all_l1_15)
stats_fista_all_F_15_mean = mean(stats_fista_all_F_15)
stats_fista_all_F_15_std = std(stats_fista_all_F_15)

stats_fista_all_fx_20_mean =  mean(stats_fista_all_fx_20)
stats_fista_all_fx_20_std =  std(stats_fista_all_fx_20)
stats_fista_all_l1_20_mean = mean(stats_fista_all_l1_20)
stats_fista_all_l1_20_std = std(stats_fista_all_l1_20)
stats_fista_all_F_20_mean = mean(stats_fista_all_F_20)
stats_fista_all_F_20_std = std(stats_fista_all_F_20)

%ISTA
    mean_snr_i =  mean(snr_i)
    std_snr_i = std(snr_i)
    mean_ssim_i =  mean(ssim_i)
    std_ssim_i = std(ssim_i)

%FISTA
    mean_snr_f =  mean(snr_f)
    std_snr_f = std(snr_f)
    mean_ssim_f =  mean(ssim_f)
    std_ssim_f = std(ssim_f)
      
%TPAL
    mean_snr_auto =  mean(snr_auto)
    std_snr_auto = std(snr_auto)
    mean_ssim_auto =  mean(ssim_auto)
    std_ssim_auto = std(ssim_auto)

%WS
    mean_snr_ws =  mean(snr_ws)
    std_snr_ws = std(snr_ws)
    mean_ssim_ws =  mean(ssim_ws)
    std_ssim_ws = std(ssim_ws)
    
%SCH
    mean_snr_scH =  mean(snr_scH)
    std_snr_scH = std(snr_scH)
    mean_ssim_scH =  mean(ssim_scH)
    std_ssim_scH = std(ssim_scH)
    
%All
    mean_snr_all = mean(snr_all)
    std_snr_all = std(snr_all)
    mean_ssim_all =  mean(ssim_all)
    std_ssim_all = std(ssim_all)
    
    
    
