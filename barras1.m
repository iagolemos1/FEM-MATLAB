%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programa SIST_BARRAS.M para análise estática por EF de sistemas
% constituídos por barras unidirecionais.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENTRADA DE DADOS
clc; clear all; close all;
%              E(N/m^2)   A(m^2)   L(m)
propriedades=[ 2.1e11     5e-4     0.2;
               2*2.1e11   5e-4     0.2;
               2.1e11     5e-4/2   2*0.2;
               3*2.1e11   5e-4     1.4*0.2];
           
F=3000;

% número de elementos
[nb_ele, dummy]=size(propriedades);

% número de nós
nb_node = nb_ele + 1;

% matriz de conectividade
mat_conect=[1 2
            2 3
            3 4
            4 5];
        
% condições de contorno nos gdl impostos [elemento deslocamento]
cond_cont=[1 0
           5 0];
       
% forças externas aplicadas nos gld livres [elemento força]
forcas_aplic= [2 F
               3 F
               4 3*F]; % valores em Newtons
           
% CONSTRUÇÃO DAS MATRIZES ELEMENTARES E MONTAGEM DA MATRIZ DE RIGIDEZ GLOBAL
K_global = zeros(nb_node); %declara matriz de rigidez global

for ii=1:nb_ele
   coef=propriedades(ii,1)*propriedades(ii,2)/propriedades(ii,3); %calculando EA/L para cada nó
   K_elementar=coef*[1 -1 ; -1 1]; %escreve a matriz de rigidez elementar [EA/L -EA/L; -EA/L EA/L]
   mat_ident=eye(nb_node);  %declara matriz identidade
   mat_transf=[mat_ident(mat_conect(ii,1),:); mat_ident(mat_conect(ii,2),:)]; %monta a matriz de transformação
   K_global= K_global + mat_transf'*K_elementar*mat_transf; %somatório para montar a matriz de rigidez global
end

% IMPOSIÇÃO DAS CONDIÇÕES DE CONTORNO PELO MÉTODO DO PARTICIONAMENTO DA
% MATRIZ DE RIGIDEZ

% identificação dos gdl livres e gdl impostos
gdl_livres = forcas_aplic(:,1); %[2 3 4]
gdl_impostos = cond_cont(:,1);  %[1 5]

% construção das submatrizes de rigidez
K_ll=K_global(gdl_livres,gdl_livres);
K_li=K_global(gdl_livres,gdl_impostos);
K_ii=K_global(gdl_impostos,gdl_impostos);

% construção dos vetores de forças nos gdl livres e de deslocamentos nos
% gdl impostos
f_liv=forcas_aplic(:,2); %[F F 3F]
d_imp=cond_cont(:,2); %[0 0]

% CÁLCULO DOS DESLOCAMENTOS NOS GDL LIVRES E FORÇAS DE REAÇÃO NOS GDL
% IMPOSTOS
d_liv = inv(K_ll)*(f_liv-K_li*d_imp);
f_imp = K_li'*d_liv+K_ii*d_imp;
d_completo=zeros(nb_node,1);
d_completo(gdl_livres)=d_liv;
delta_d=diff(d_completo);
for ii=1:nb_ele
sigma(ii)=delta_d(ii)*propriedades(ii,1)/propriedades(ii,3);
end

% APRESENTAÇÃO DOS RESULTADOS
disp('Deslocamentos(m)')
d_completo
disp('Reações(N)')
f_imp
disp('Tensões(N/m^2)')
sigma