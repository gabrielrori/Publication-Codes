%Función A'*x, empleada para el cálculo de la gradiente (A'*(A*x-b))
function ATx = ATx_fun(z, opt, var)
    switch opt.dict
        
        %En el caso de datos sintéticos aleatorios, debido a que la matriz A 
        %es conocida, se tiene multiplicación de la matriz A' por el vector
        %z = A*x-b
        case 'rand'
            ATx.val = var.A'*z; 
        
        %En el caso de imágenes, se toma la Transformada Discreta Wavelets 2D
        case 'wave'
            [ATx.val,ATx.size] = wavedec2(z,opt.n,opt.wname);
            ATx.val = ATx.val';
    end

end