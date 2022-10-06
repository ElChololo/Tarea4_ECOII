classdef TestHausmanAlternativo
    %Objeto para realizar test de Hausman
    
    properties
        % X deben ser definidas como las variables que se desee probar
        % exogeneidad
        X
        %% Z variables exógenas a X
        Z
        %% Variable que queremos demostar que depende de Y
        Y
    end
    
    methods
        function obj = TestHausmanAlternativo(X,Z,Y)
            %Se definen los vectores de variables que definirán las
            %propiedades X,Z,Y
            obj.X = X;
            obj.Z = Z;
            obj.Y=Y;
        end
        
        function v_est = reg_xz(obj)
            %Función que devuelve un vector de errores v estimados.
            %Provienen de los residuos de una regresión de X en z
            y=obj.X;
            x=[ones(length(obj.Z),1) obj.Z];
            coef_mco=(x'*x)\(x'*y);
            v_est=y-coef_mco*x;
            
        end
        function u_est = reg_yx(obj)
            y=obj.Y;
            x=[ones(length(obj.X,1)) obj.X];
            coef_mco=(x'*x)\(x'*y);
            u_est = y-coef_mco*coef_mco;
        end
        
        function bool_test_t=reg_u_xv(obj,u_est,v_est)
            y=u_est;
            x=obj.X;
            x=[ones(length(v_est),1) x v_est];
            coef_mco = (x'*x)\(x'*y);
            resid = y-coef_mco*x;
            resid_2 = resid'*resid;
            sigma_2= resid_2/(size(x,1)-size(x,2));
            %Ahora hay que ver que TODOS los coeficientes asociados a los
            %v_est sean significativos
            var_betas = sigma_2*inv((x'*x));
            %Identificar las posiciones de los estimadores v_est
            posicion_v_est=size(x,2)-size(obj.X,2)-1;
            significancia_v_est=zeros(size(v_est,2),1);
            %% CREO QUE NO ESTÁ AGARRANDO BIEN LOS COEF MCO RELEVANTES
            for ii=1:posicion_v_est
                est_t= (coef_mco(ii))/(var_betas(ii,ii))^(1/2);
                if abs(est_t)>1.96
                    significancia_v_est(ii)=0;
                else
                    significancia_v_est(ii)=1;
                end
            end
            %% Ahora vemos la significancia de los test, si al menos 1 estadistico no era significativo. 
            %%la función devolverá un valor 1 que significa que un v_est
            %%fue significativo en la regresión u_t. Por tanto, se rechaza
            %%la nula de exogeneidad de X
            if sum(significancia_v_est) ~=0
                bool_test_t=1;
            else
                bool_test_t=0;
            end
            
        end
    end
end

