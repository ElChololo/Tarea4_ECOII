classdef Test_White_homocedasticidad
    %Objeto diseñado para testear homocedasticidad utilizando el test de
    %White. Entregará el valor 1 si pasa el test de homocedasticidad, 0 si
    %se rechaza la hipótesis de homocedasticidad
    
    properties
        Y
        X
        u_est
    end
    
    methods
        function obj = Test_White_homocedasticidad(Y,X)
            % Constructor del objeto
            obj.Y = Y;
            obj.X=X;
        end
        
        function u_est = Regresion_Auxiliar(obj,Y,X)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            x=[ones(size(X,1)) X];
            coef_mco = (x'*x)\(x'*Y);
            u_est=Y-coef_mco*x;
            obj.u_est=u_est;
            
        end
        
        function Bool_Test_pass = Test_White(obj,X,u_est)
            %Se regresiona los residuos al cuadrado
            var_depe=u_est.^2;
            %Armar las variables al cuadrado
            regresores_cuadrado=zeros(size(X,1),size(X,2));
            for ii=1:size(X,2)
                regresores_cuadrado(ii)=X(1:end,ii).^2;
            end
            %Armar los productos cruzados
            regresores_cruzados=zeros(size(X,1),nchoosek(size(X,2),2));
            combinaciones=nchoosek(1:size(X,2),2);
            for ii=1:size(combinaciones,1)
               regresores_cruzados=X(1:end,combinaciones(ii,1)).* X(1:end,combinaciones(ii,2)); 
            end
            % Armar los regresores
            x=[ones(size(X,1),1) X regresores_cruzados regresores_cuadrado];
            coef_mco=(x'*x)/(x'*var_depe);
            resid=var_depe-coef_mco*x;
            RSS=resid'*resid;
            TSS_aux= var_depe-mean(var_depe);
            TSS= TSS_aux'*TSS_aux;
            R_2= 1- (RSS/TSS);
            LM= size(var_depe,1)* R_2;
            val_crit=chi2inv(0.95,size(x,2)-1);
            if LM > val_crit
                %se rechaza la hipotesis nula de homocedasticidad
                Bool_Test_pass = 0;
            else
                Bool_Test_pass= 1;
            end
        end
    end
end

