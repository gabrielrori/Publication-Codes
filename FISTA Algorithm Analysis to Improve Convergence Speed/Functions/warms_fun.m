%Función de Warm Start
function lambda_vec =  warms_fun(opt, var,i)

    if i > (opt.warms_iter)
        lambda_vec = var.lambda_vec_init; %vector inicial de lambdas
    else
         switch opt.wtype
             
             %Decrecimiento lineal hallado en base a dos puntos
             case 'line'
                 lambda_vec = var.lambda_vec_ws*(((1-opt.ws_inc)*i+(opt.warms_iter*opt.ws_inc+opt.ws_inc-1))/(opt.warms_iter*opt.ws_inc));
             
             %Decrecimiento cuadrático hallado en base al vertice
             %(derivada = 0) y otro punto
             case 'cuad'
                 lambda_vec = var.lambda_vec_ws*(((1-opt.ws_inc)*i*i+(opt.ws_inc-1)*2*i+(opt.warms_iter*opt.warms_iter*opt.ws_inc+1-opt.ws_inc))/(opt.warms_iter*opt.warms_iter*opt.ws_inc));
             
             %Caso con lambda continuo, sin decrecimiento
             case 'cont'
                 lambda_vec = var.lambda_vec_ws;
         end
    end
end