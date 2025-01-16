%Función para la generación de datos iniciales y variables necesrias
function [x_orig, b, x, var] = data_gen(opt)
    switch opt.dict
        
        %Caso de datos sintéticos aleatorios
        case 'rand'
            %Matriz A aleatoria normalizada
            var.A = sqrt(1/opt.N)*randn(opt.N, opt.N);
            %Datos originales (vector de datos aleatorios)
            x_orig = randn(opt.N,1);
            %Generar datos originales que sean sparse 
            ind = randperm(opt.N); 
            R = 0.25*(opt.N); 
            x_orig(ind(1:R)) = 0; 
            %Datos observados b
            b = var.A*x_orig + opt.sigma*randn(opt.N,1);
            %Datos iniciales aleatorios
            x.val = randn(opt.N,1);            
            %Vector de lambdas (necesario para shrinkage y Warm Start)
            var.c_mask = ones(opt.N,1);
            var.lambda_vec = opt.lambda*var.c_mask;
        
        %Caso de imágenes
        case 'wave'
            %Datos originales: lectura de imagen
            x_orig = imread(opt.sample_img);
            S = length(x_orig);
            x_orig = double(x_orig)./255.0;
            %Imagen observada: imágene original y suma de ruido Gaussiano
            b = x_orig;
            b = b + opt.sigma*randn(S,S);
            %Datos iniciales de las mismas dimensiones que la transformada
            %directa 2D de los datos observados
            [x.val,x.size] = wavedec2(b,opt.n,opt.wname);
            opt.N = length(x.val);
            x.val = randn(opt.N, 1);
            %Vector de lambdas necesario para shrinkage y Warm Start:
            %cantidad de coeficientes de aproximación
            cA_size = x.size(1,1)*x.size(1,2);
            %cantidad de coeficientes de detalle
            cD_size = sum(x.size(2:opt.n+1,1).*x.size(2:opt.n+1,2))*3; 
            var.c_mask = [zeros(cA_size,1); ones(cD_size,1)];
            %lambda solo influye en los coeficientes de detalle
            var.lambda_vec = opt.lambda*var.c_mask;
    end
end