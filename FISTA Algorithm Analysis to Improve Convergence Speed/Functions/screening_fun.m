%Función de Screening
function mask_SC =  screening_fun(b, opt, var, x)
    switch opt.dict
        
        %Caso de datos sintéticos aleatorios
        case 'rand'
            lambda_max = max(b'*var.A);%Calculo de lambda_max
            lambda = opt.screeningfeatures*lambda_max; %Porcentaje de
            %lambda_max a ser evaluado (diferentes intensidades de screening)
            normb = norm(b); %norma l_2 de b
            normAi = 0.0394*vecnorm(var.A,1)'; %aproximacion de la norma l2
            %de las columnas de A  en base a la norma l1 (menor tiempo de
            %computo)
            absbA = abs(b'*var.A)'; 
            pk = (normb.*normAi+absbA)./(normb.*normAi+lambda_max); %Cálculo
            %del vector pk
            pklambdaMax = pk*lambda_max; %Vector pk*lambda_max
            %Comparacion de lambda y pk*lambda_max
            mask_SC = lambda < pklambdaMax ; %máscara lógica de 0's y 1's
            mask_SC = double(mask_SC); %máscara de doubles
            
        %Caso de imágenes    
        case 'wave'
            c = x.val; %valores de la transformada wavelet de la imágen
            s = x.size; %tamaño de los niveles de la transformada
            c_size = length(c); %tamaño de todos los coeficientes
            cA_size = s(1,1)*s(1,2); %tamaño de los coeficientes aproximados
            cD_size = c_size - cA_size; %tamaño de los coeficientes de detalle
            c_mask = [zeros(cA_size,1); ones(cD_size,1)]; %máscara de 0's
            %(coeficientes aproximados) y 1's (coeficientes de detalle)
            c = c_mask.*c; %se evaluan solo los coeficiente de detalle
            lambda_max = max(c); %maximo valor de los coeficientes de detalle
            lambda = opt.screeningfeatures*lambda_max;  %Porcentaje de
            %lambda_max a ser evaluado (diferentes intensidades de screening)
            xk2 = (norm(c)/norm(b))/length(c); %aproximación del cálculo de
            %la norma de las columnas de A
            y2 = norm(b) ; 
            pk = (y2*xk2+ abs(c))/(y2*xk2+ lambda_max); %Cálculo del vector pk
            pklambdaMax = pk*lambda_max; %Vector pk*lambda_max
            %Comparacion de lambda y pk*lambda_max
            mask_SC = lambda < pklambdaMax ;  %máscara lógica de 0's y 1's
            mask_SC = double(mask_SC); %máscara de doubles
            mask_SC = mask_SC .*c_mask; %máscara que solo toma los
            %coeficientes de detalle)
            mask_SC(1:cA_size) = ones(cA_size, 1); %máscara decoefientes
            %aproximados sin cambiar (1's)
    end
end