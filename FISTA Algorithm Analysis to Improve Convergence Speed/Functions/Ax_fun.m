%Función A*x, empleada para la gradiente y las funciones de Error y Costo
function Ax = Ax_fun(z, opt, var)
    switch opt.dict
        
        %En el caso de datos sintéticos aleatorios, debido a que la matriz A
        %es conocida, se tiene multiplicación de la matriz A por el vector x
        case 'rand'
            Ax = var.A*z.val;
        
        %En el caso de imágenes, se toma la Transformada Discreta Wavelets
        %2D Inversa
        case 'wave'
            Ax = waverec2(z.val,z.size,opt.wname);
    end

end