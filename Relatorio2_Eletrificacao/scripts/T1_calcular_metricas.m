function [coef,E,R,M] = T1_calcular_metricas(coef,E,R)
% Erros sempre entre teoria e FFT da simulacao; THD exclui a componente CC.
coef.erro_par_V = hypot(coef.a_fft-coef.a_teo,coef.b_fft-coef.b_teo);
coef.erro_C_V = abs(coef.C_fft-coef.C_teo);
E.erro_total_pct = 100*abs(E.E_total_sim-E.E_total_angular)./E.E_total_angular;
E.erro_E1_pct = 100*abs(E.E1_sim-E.E1)./E.E1;
E.erro_E1_7_pct = 100*abs(E.E1_7_sim-E.E1_7)./E.E1_7;
E.THD7_teo_pct = 100*sqrt(max(0,E.E1_7-E.E1)./E.E1);
E.THD7_sim_pct = 100*sqrt(max(0,E.E1_7_sim-E.E1_sim)./E.E1_sim);
E.THD_total_teo_pct = 100*sqrt(max(0,E.E_total_angular-E.E_CC-E.E1)./E.E1);
E.THD_total_sim_pct = 100*sqrt(max(0,E.E_total_sim-E.E_CC_sim-E.E1_sim)./E.E1_sim);
R.erro_a0_V = abs(R.a0_sim-R.a0_teo);
R.erro_C1_pct = 100*abs(R.C1_sim-R.C1_teo)./R.C1_teo;
M = table();
M.Caso = E.Caso;
M.erro_a0_V = R.erro_a0_V;
M.erro_C1_pct = R.erro_C1_pct;
M.erro_coef_max_V = zeros(height(E),1);
for j=1:height(E)
    f=coef(coef.Caso==E.Caso(j),:);
    assert(height(f)==7,'Esperados sete pares de coeficientes por caso.');
    M.erro_coef_max_V(j)=max(f.erro_par_V);
end
M.erro_Etotal_pct=E.erro_total_pct;
M.erro_E1_pct=E.erro_E1_pct;
M.erro_E1_7_pct=E.erro_E1_7_pct;
M.THD7_teo_pct=E.THD7_teo_pct;
M.THD7_sim_pct=E.THD7_sim_pct;
M.THD_total_teo_pct=E.THD_total_teo_pct;
M.THD_total_sim_pct=E.THD_total_sim_pct;
M.diferenca_THD7_pp=abs(M.THD7_sim_pct-M.THD7_teo_pct);
M.diferenca_THD_total_pp=abs(M.THD_total_sim_pct-M.THD_total_teo_pct);
end
