%Función general de ISTA
function [stats, x, var] = general_ista(opt, var, b, x, L)
    
    %Cálculo de la primera estadística previo a computar el algoritmo
    Ax = Ax_fun(x, opt, var);
    z = Ax- b;
    stats.fx(1) = 0.5*sum(z(:).*z(:)); %Función de Error
    stats.l1(1) = norm(x.val.*var.c_mask,1); %Función de la norma l1 (tamaño)
    stats.F(1) = stats.fx(1) + opt.lambda.*stats.l1(1); %Función de Costo
    
    %Inicio de conteo de tiempo
    var.time(1) = 0;
    tic
    
    %Creación del vector de lambdas para Warm Start
    if opt.warms
        var.lambda_vec_init = var.lambda_vec;
        var.lambda_vec_ws = opt.ws_inc*var.lambda_vec;
    end
    maskS = ones(length(x.val),1);
    
    %Reducción de tamaño por Screening
    if opt.screening
        var.mask_SC = screening_fun(b, opt, var, x);
        x.val = x.val.*var.mask_SC;
    end
        
    %Primer paso de Cauchy
    i = 0;
    ATz = ATx_fun(z, opt, var);
    g = ATz;
    var.alpha = step_fun(x, g, opt, var, i);

    %Inicio del algoritmo iterativo
    for i = 1:L

        %Warm Start
        if opt.warms
            var.lambda_vec =  warms_fun(opt, var,i);
        end
        
        %Guardar valores anteriores
        var.x_ant = x.val;
        var.g_ant = g.val;

        %Shrinkage
        arg_p = x.val - var.alpha.*g.val;
        v = abs(arg_p);
        mask = v > var.alpha.*var.lambda_vec;
        mask = mask.*maskS;
        
        %Actualización de valores
        x.val = sign(arg_p).*(mask.*v - var.alpha.*var.lambda_vec);

        %Gradiente
        Ax = Ax_fun(x, opt, var);
        z = Ax-b;
        ATz  = ATx_fun(z, opt, var);
        g = ATz;

        %Tamano de paso BB1
        var.alpha = step_fun(x, g, opt, var, i);    

        %Estadisticas
        stats.fx(i+1) = 0.5*sum(z(:).*z(:)); %Función de Error
        stats.l1(i+1) = norm(x.val.*var.c_mask,1);  %Función de la norma l1
        stats.F(i+1) = stats.fx(i+1) + opt.lambda.*stats.l1(i+1); %F. de Costo
        maskS = ones(length(x.val),1);
        var.time(i+1) = toc;

    end
end