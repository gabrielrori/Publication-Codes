function alpha = step_fun(x, g, opt, var, i, alphaCi)

    %Si se encuentra en la iteracion 0, se calcula el tamaño de paso de Cauchy 
    if ((isequal(opt.stepT,'BB1') ) || (isequal(opt.stepT,'YnB') )) && i == 0
        Ag = Ax_fun(g, opt, var); %Cálculo de A*g
        alpha0 = norm(g.val)^2/norm(Ag)^2; %Tamaño de paso de Cauchy
        alpha0 = opt.step0_mult*alpha0; %Multiplicación por un coeficiente
        %Limites para el tamaño de paso  (opcional)
        if opt.steplim
            if alpha0 > opt.stepCH 
                alpha0 = opt.stepCH ;
            end
            if alpha0 < opt.stepCL 
                alpha0 = opt.stepCL;
            end
        end    
    end
    
	switch opt.stepT
        
        %Tamaño de paso constante
        case 'cte'
        alpha = opt.alpha;
        
        %Tamaño de paso de Barzilai-Borwein 
        case 'BB1'
        if i == 0
            alpha = alpha0; %Si es iteración 0, usar Cauchy
        else
            %Diferencia entre gradiente y valores de iteración actual y
            %pasada
            s = x.val - var.x_ant;
            y = g.val - var.g_ant;
            alpha = s'*y/norm(y)^2; %Tamaño de paso BB1
            alpha = opt.step_mult*alpha; %Multiplicación por un coeficiente
            %Limites para el tamaño de paso  (opcional)
            if opt.steplim
                if alpha > opt.stepH 
                    alpha = opt.stepH ;
                end
                if alpha < opt.stepL
                    alpha = opt.stepL;
                end
            end
        end
        
        %Tamaño de paso de Yuan B
        case 'YnB'
        if i == 0
            alpha = alpha0;  %Si es iteración 0, usar Cauchy
        else
            %Intercalamiento  entre tamaño de paso de Cauchy y Yuan
            f = mod(i,3);
            if i == 3 ||  f
                %Tamaño de paso de Cauchy
                alpha = alpha0;
            else
                %Tamaño de paso de Yuan B
                s = x.val - var.x_ant;
                alpha = 2/((((1/alphaCi) -1/(alpha0))^2+(4*norm(g.val)^2)/(norm(s)^2))^(1/2) + 1/alphaCi+1/alpha0);
            end  
            alpha = opt.step_mult*alpha; %Multiplicación por un coeficiente
            %Limites para el tamaño de paso  (opcional)
            if opt.steplim
                if alpha >  opt.stepH 
                    alpha =  opt.stepH ;
                end
                if alpha <  opt.stepL
                    alpha =  opt.stepL;
                end
            end
        end
    end
end