classdef Test_Chow_Sequencial
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        obs
        Y
        X
    end
    
    methods
        function obj = Test_Chow_Sequencial(Y,X)
           obj.obs=length(Y);
           obj.Y=Y;
           obj.X=X;
        end
        
        function Bool_test = Test_Chow(obj)
            % Retorna el valor 1 si pasó el test de Chow, es decir
            % No se rechaza la nula de que no hubieron quiebres
            % estructurales en algún período
            Bool_test_Mat=zeros(len(obj.Y)-2,1);
            for ii=2:obj.obs-1
               Bool_test_Mat(ii)=obj.test_chow(ii,obj.Y,obj.X); 
            end
            if sum(Bool_test_Mat)==0
                Bool_test = 1;
            else
                Bool_test = 0;
            end
                
        end
        
        function bol_est_chow=test_chow(obj,periodo_quiebre,Y,X)
            x_prequiebre=X(1:periodo_quiebre,1:end);
            x_postquiebre=X(periodo_quiebre+1:end,1:end);
            zeros_aux=zeros(periodo_quiebre,size(X,2));
            mat_sup=[x_prequiebre zeros_aux];
            mat_inf=[zeros_aux x_postquiebre];
            reg= [mat_sup;mat_inf];
            st_h1= Y'*( eye(size(reg,1))-( reg/(reg'*reg)*reg' ))*Y;
            st_h0=Y'*( eye(size(X,1))-( X/(X'*X)*X' ))*Y;
            aux_est_chow=(st_h1-st_h0)/st_h1;
            est_chow= ((size(X,1)-2*size(X,2)) / size(X,2) )* aux_est_chow;
            valor_crit_f= finv(0.95,5,size(x,1)-2*size(x,2));
            %Los test F se hacen a una cola
            if est_chow > valor_crit_f
                %se rechaza la nula
                bol_est_chow=1;
            else
                bol_est_chow=0;
            end
            
        end
    end
end

