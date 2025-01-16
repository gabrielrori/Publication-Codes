%Función general de FISTA
function [stats, x, var] = general_fista(opt, var, b, x, L)
    
    %Cálculo de la primera estadística previo a computar el algoritmo
    Ax = Ax_fun(x, opt, var);
    z = Ax- b;
    stats.fx(1) = 0.5*sum(z(:).*z(:)); %Función de Error
    stats.l1(1) = norm(x.val.*var.c_mask,1); %Función de la norma l1 (tamaño)
    stats.F(1) = stats.fx(1) + opt.lambda.*stats.l1(1); %Función de Costo

    %Inicio de conteo de tiempo
    var.time(1) = 0;
    tic
    y = x; 
    t = 1;
        
    %Creación del vector de lambdas para Warm Start
    if opt.warms
        var.lambda_vec_init = var.lambda_vec;
        var.lambda_vec_ws = opt.ws_inc*var.lambda_vec;
    end
    maskS = ones(length(x.val),1);
    
    %Reducción de tamaño por Screening
    if opt.screening
        var.mask_SC = screening_fun(b, opt, var, y);
        x.val = x.val.*var.mask_SC;
        y.val = y.val.*var.mask_SC;
    end
        
    %Primer paso de Cauchy
    i = 0;
    ATz = ATx_fun(z, opt, var);
    g = ATz;
    var.g(1,:) = g.val;
    var.alpha(1) = step_fun(y, g, opt, var, i);
    alphaCi = var.alpha(1);
    
    %Inicio del algoritmo iterativo
    for i = 1:L

        %Warm Start
        if opt.warms
            var.lambda_vec =  warms_fun(opt, var,i);
        end
        
        %Cálculo de tamaño de paso de Cauchy de iter. previa para Yuan B
        if opt.stepT == 'YnB' 
            ATz = ATx_fun(z, opt, var);
            g = ATz;
            alphaCi = 0;
            alphaCi = step_fun(z, y, g, opt, var, i, alphaCi);
        end

        %Guardar valores anteriores
        var.x_ant = x.val;
        var.g_ant = g.val;
        var.y_ant = y.val;

        %Shrinkage
        arg_p = y.val - var.alpha(i).*g.val;
        v = abs(arg_p);
        mask = v > var.alpha(i).*var.lambda_vec;
        mask = mask.*maskS;
        x.val = sign(arg_p).*(mask.*v - var.alpha(i).*var.lambda_vec);
        
        %Actualización de valores
        t_next = (1+(1+4*t^2)^(1/2))/2;
        y.val = x.val + ((t-1)/t_next).*(x.val-var.x_ant);
        t = t_next;
          
        %Gradiente
        Ax = Ax_fun(y, opt, var);
        z = Ax-b;
        ATz  = ATx_fun(z, opt, var);
        g = ATz;
        var.g(i+1,:) = g.val;

        %Tamano de paso BB1
        var.alpha(i+1) = step_fun(y, g, opt, var, i , alphaCi);    

        %Cálculo de estadisticas
        stats.fx(i+1) = 0.5*sum(z(:).*z(:)); %Función de Error
        stats.l1(i+1) = norm(y.val.*var.c_mask,1); %Función de la norma l1
        stats.F(i+1) = stats.fx(i+1) + opt.lambda.*stats.l1(i+1); %F. de Costo
        maskS = ones(length(y.val),1);
        var.time(i+1) = toc;
    end
end