function T1_verificar_resultados
% Auditoria: quadratura independente, FFT das amostras e numeros dos arquivos TeX.
raiz=fileparts(fileparts(mfilename('fullpath'))); pasta=fullfile(raiz,'resultados');
S=load(fullfile(pasta,'validacao_completa.mat'));
C=readtable(fullfile(pasta,'coeficientes_fft.csv'),'TextType','string');
E=readtable(fullfile(pasta,'energias.csv'),'TextType','string');
R=readtable(fullfile(pasta,'validacao_fft.csv'),'TextType','string');
M=readtable(fullfile(pasta,'metricas_erro_thd.csv'),'TextType','string');
assert(numel(S.dadosOndas)==24 && height(E)==24 && height(C)==168);
maxIntegral=0; maxFFT=0; maxParseval=0;
for j=1:24
    q=S.dadosOndas{j}; c=q.caso; N=numel(q.y);
    assert(N==S.N && all(isfinite(q.y)));
    f=C(C.Caso==string(c.nome),:);
    assert(isequal(f.n,(1:7)') && E.Caso(j)==string(c.nome));
    if c.tipo<=6
        pontos=unique([0,c.L,mod(c.phi,c.L),mod(c.phi+c.L/2,c.L)]);
    else
        pontos=[0,c.phi,c.phi+c.d,c.L];
        if c.alternado, pontos=[pontos,c.P+c.phi2,c.P+c.phi2+c.d2]; end
        pontos=unique(pontos);
    end
    fun=@(x) valorOnda(x,c);
    a0=quadratura(fun,pontos)/c.L;
    Et=quadratura(@(x)fun(x).^2,pontos);
    a=zeros(7,1); b=a;
    for n=1:7
        k=2*pi*n/c.L;
        a(n)=2/c.L*quadratura(@(x)fun(x).*cos(k*x),pontos);
        b(n)=2/c.L*quadratura(@(x)fun(x).*sin(k*x),pontos);
    end
    erro=max(abs([a0-R.a0_teo(j);a-f.a_teo;b-f.b_teo;Et-E.E_total_angular(j)]));
    maxIntegral=max(maxIntegral,erro);
    assert(erro<1e-10,'Integral teorica inconsistente: %s',c.nome);
    Y=fft(q.y)/N; cn=Y(2:8).*exp(-1i*pi*(1:7)'/N);
    af=2*real(cn); bf=-2*imag(cn);
    Ecc=c.L*a0^2; E1=c.L/2*(a(1)^2+b(1)^2);
    E7=c.L/2*sum(a.^2+b.^2);
    Es=c.L*mean(q.y.^2); Eccs=c.L*real(Y(1))^2;
    E1s=c.L/2*(af(1)^2+bf(1)^2); E7s=c.L/2*sum(af.^2+bf.^2);
    esperado=[Et,Ecc,E1,E7,Ecc+E7,100*E1/Et,100*E7/Et,100*(Ecc+E7)/Et,Et/S.omega, ...
        Es,Eccs,E1s,E7s,Eccs+E7s,100*E1s/Es,100*E7s/Es,100*(Eccs+E7s)/Es,Es/S.omega];
    publicado=E{j,{'E_total_angular','E_CC','E1','E1_7','E_CC_7','eta1','eta1_7', ...
        'etaCC_7','E_total_temporal','E_total_sim','E_CC_sim','E1_sim','E1_7_sim', ...
        'E_CC_7_sim','eta1_sim','eta1_7_sim','etaCC_7_sim','E_total_temporal_sim'}};
    erro=max(abs([af-f.a_fft;bf-f.b_fft;real(Y(1))-R.a0_sim(j); ...
        esperado(:)-publicado(:)]));
    maxFFT=max(maxFFT,erro);
    assert(erro<1e-9,'FFT ou energia inconsistente: %s',c.nome);
    assert(max(abs(f.C_teo-hypot(a,b)))<1e-11);
    assert(max(abs(f.C_fft-hypot(af,bf)))<1e-11);
    assert(abs(R.C1_sim(j)-hypot(af(1),bf(1)))<1e-11);
    assert(abs(R.Pq_sim(j)-Es/c.L)<1e-11);
    assert(abs(E.f0_Hz(j)-S.omega/c.L)<1e-11);
    % A soma direta do espectro inclui o bin de Nyquist uma unica vez.
    Eespectro=c.L*sum(abs(Y).^2);
    maxParseval=max(maxParseval,abs(Eespectro-Es));
    assert(abs(Eespectro-Es)<1e-10,'Falha de Parseval discreto.');
    thdFFT=100*sqrt((2*sum(abs(Y(3:N/2)).^2)+abs(Y(N/2+1))^2)/(2*abs(Y(2))^2));
    assert(abs(thdFFT-M.THD_total_sim_pct(j))<1e-9,'THD total nao confere.');
    assert(abs(100*sqrt(sum(f.C_fft(2:7).^2))/f.C_fft(1)-M.THD7_sim_pct(j))<1e-9);
    assert(E.E_CC_7(j)<=Et+1e-10 && E.E_CC_7_sim(j)<=Es+1e-10);
    % Os dados das figuras devem usar as mesmas amostras e reconstrucoes.
    W=readtable(fullfile(pasta,sprintf('onda_%02d.csv',j)));
    indices=unique([1:128:N,N]); u=q.x(indices)/c.L;
    teoria7=a0+zeros(size(u)); fft7=real(Y(1))+zeros(size(u));
    for n=1:7
        teoria7=teoria7+a(n)*cos(2*pi*n*u)+b(n)*sin(2*pi*n*u);
        fft7=fft7+af(n)*cos(2*pi*n*u)+bf(n)*sin(2*pi*n*u);
    end
    assert(max(abs(W{:,:}-[u,q.y(indices),teoria7,fft7]),[],'all')<1e-10, ...
        'Dados da figura da onda nao conferem: %s',c.nome);
    V=readtable(fullfile(pasta,sprintf('espectro_caso_%02d.csv',j)));
    assert(max(abs(V{:,:}-[(1:7)',hypot(a,b),hypot(af,bf)]),[],'all')<1e-10);
end
[Ci,Ei,Ri,Mi]=T1_calcular_metricas(C,E,R);
assert(max(abs(Ci.erro_par_V-C.erro_par_V))<1e-12);
assert(max(abs(Ei.THD_total_teo_pct-E.THD_total_teo_pct))<1e-10);
assert(max(abs(Ri.erro_C1_pct-R.erro_C1_pct))<1e-10);
assert(isequal(Mi.Properties.VariableNames,M.Properties.VariableNames));
assert(max(abs(Mi{:,2:end}-M{:,2:end}),[],'all')<1e-10);
% Conferir invariancias de amplitudes, percentuais e mudanca de periodo.
for j=1:6
    assert(max(abs([E.eta1(j)-E.eta1(j+6),E.eta1(j)-E.eta1(j+12)]))<1e-10);
    assert(abs(E.E_total_angular(j+12)/E.E_total_angular(j)-0.8)<1e-12);
end
assert(abs((3-sqrt(6))-(1/4+(sqrt(3)/2-sqrt(2))^2))<1e-14);
assert(abs(E.eta1(21)-120*(3-sqrt(6))/pi^2)<1e-10);
assert(abs(E.eta1(24)-200*(3-sqrt(6))/pi^2)<1e-10);
% Verificar todos os numeros tabulados, inclusive seus arredondamentos.
quantidade=0;
quantidade=quantidade+conferirTex(pasta,'tabela_validacao.tex', ...
    [R.a0_teo,R.a0_sim,R.C1_teo,R.C1_sim]);
quantidade=quantidade+conferirTex(pasta,'tabela_coeficientes_comparacao.tex', ...
    [C.a_teo,C.a_fft,C.b_teo,C.b_fft,C.erro_par_V]);
esperado=[];
for g=1:4
    ids=(g-1)*6+(1:6);
    for j=ids
        esperado=[esperado;E{j,{'E_total_angular','E1','E1_7'}}; ...
            E{j,{'E_total_sim','E1_sim','E1_7_sim'}}]; %#ok<AGROW>
    end
    for j=ids
        esperado=[esperado;E{j,{'eta1','eta1_7','etaCC_7'}}; ...
            E{j,{'eta1_sim','eta1_7_sim','etaCC_7_sim'}}]; %#ok<AGROW>
    end
end
quantidade=quantidade+conferirTex(pasta,'tabela_energias.tex',esperado);
quantidade=quantidade+conferirTex(pasta,'tabela_erro_thd.tex', ...
    M{:,{'erro_Etotal_pct','THD7_teo_pct','THD7_sim_pct','THD_total_teo_pct','THD_total_sim_pct'}});
T0=[S.dadosOndas{1}.caso.L;S.dadosOndas{13}.caso.L;S.dadosOndas{21}.caso.L]/S.omega;
quantidade=quantidade+conferirTex(fullfile(raiz,'Capitulos'), ...
    'resultados_implementacao.tex',[1e3*T0,1e9*T0/S.N]);
for j=1:3
    ids=[19 20 21]; f=C(C.Caso==E.Caso(ids(j)),:);
    quantidade=quantidade+conferirTex(pasta,sprintf('retangular_coef_%d.tex',j), ...
        [f.a_teo,f.b_teo]/S.dadosOndas{ids(j)}.caso.A);
end
config=readtable(fullfile(pasta,'configuracao_simulacao.csv'),'TextType','string');
assert(height(config)==24 && all(config.N==S.N) && all(config.Ciclos_FFT==1));
assert(max(abs(config.Ts_s-E.L/(S.omega*S.N)))<1e-16);
assert(all(config.Simulink==S.temSimulink));
if S.temSimulink
    modelo='T1_validacao_sinais';
    assert(~bdIsLoaded(modelo),'Feche o modelo antes de verificar os parametros salvos.');
    load_system(fullfile(pasta,[modelo,'.slx']));
    limpeza=onCleanup(@()close_system(modelo,0));
    assert(strcmp(get_param(modelo,'Solver'),'FixedStepDiscrete'));
    assert(strcmp(get_param(modelo,'SolverType'),'Fixed-step'));
    assert(abs(str2double(get_param(modelo,'FixedStep'))-config.Ts_s(end))<1e-16);
    assert(abs(str2double(get_param([modelo '/Tempo'],'SampleTime'))-config.Ts_s(end))<1e-16);
    assert(abs(str2double(get_param(modelo,'StopTime'))-config.StopTime_s(end))<1e-14);
    assert(strcmp(get_param([modelo '/Amostras'],'SaveFormat'),'Timeseries'));
end
fid=fopen(fullfile(pasta,'conferencia_numerica.txt'),'w');
fprintf(fid,'PASSOU: 24 cenarios; 168 pares, ordens 1 a 7.\n');
fprintf(fid,'Origem registrada das amostras: Simulink=%d.\n',S.temSimulink);
fprintf(fid,'Quadratura por trechos vs teoria: diferenca maxima %.12g.\n',maxIntegral);
fprintf(fid,'FFT recalculada e energias vs CSV: diferenca maxima %.12g.\n',maxFFT);
fprintf(fid,'Parseval discreto: diferenca maxima %.12g V^2 rad.\n',maxParseval);
fprintf(fid,'Conferidos %d valores nas 44 tabelas numericas TeX, com arredondamento.\n',quantidade);
fprintf(fid,'Conferidos os CSV das 24 ondas e dos 24 espectros apresentados nas figuras.\n');
fprintf(fid,'THD parcial e total verificados por soma direta do espectro, com Nyquist.\n');
fprintf(fid,'G=3-sqrt(6); porcentagens das alternadas verificadas.\n');
fprintf(fid,'Erro maximo por par de coeficientes: %.12g V.\n',max(M.erro_coef_max_V));
fprintf(fid,'Erro maximo da media: %.12g V.\n',max(M.erro_a0_V));
fprintf(fid,'Erro relativo maximo C1: %.12g %%.\n',max(M.erro_C1_pct));
fprintf(fid,'Erro relativo maximo Etotal: %.12g %%.\n',max(M.erro_Etotal_pct));
fprintf(fid,'Erro relativo maximo E1: %.12g %%.\n',max(M.erro_E1_pct));
fprintf(fid,'Erro relativo maximo E1:7: %.12g %%.\n',max(M.erro_E1_7_pct));
fprintf(fid,'Diferenca maxima do THD total: %.12g pontos percentuais.\n',max(M.diferenca_THD_total_pp));
fclose(fid);
disp(fileread(fullfile(pasta,'conferencia_numerica.txt')));
end

function y=valorOnda(x,c)
% Definicoes por trechos, sem usar coeficientes analiticos ou FFT.
if c.tipo<=6
    z=mod(x-c.phi,c.L);
    switch c.tipo
        case 1, y=c.A*z/c.L;
        case 2, y=c.A-c.A*z/c.L;
        case 3
            y=2*c.A*z/c.L; ids=z>=c.L/2;
            y(ids)=2*c.A*(c.L-z(ids))/c.L;
        case 4, y=-c.A+2*c.A*z/c.L;
        case 5, y=c.A-2*c.A*z/c.L;
        case 6
            y=-c.A+4*c.A*z/c.L; ids=z>=c.L/2;
            y(ids)=3*c.A-4*c.A*z(ids)/c.L;
    end
else
    pulso=x>=c.phi & x<c.phi+c.d;
    if c.alternado, pulso=pulso | (x>=c.P+c.phi2 & x<c.P+c.phi2+c.d2); end
    y=c.A*double(pulso);
    if c.tipo==8, y=2*y-c.A; end
end
end

function v=quadratura(fun,pontos)
v=0;
for j=1:numel(pontos)-1
    v=v+integral(fun,pontos(j),pontos(j+1),'AbsTol',1e-12,'RelTol',1e-12);
end
end

function quantidade=conferirTex(pasta,nome,esperados)
s=fileread(fullfile(pasta,nome));
tabelas=regexp(s,'(?s)\\begin\{table\}.*?\\end\{table\}','match');
s=strjoin(tabelas,newline);
trechos=regexp(s,'\$([^$]*)\$','tokens'); valores={};
for j=1:numel(trechos)
    t=strrep(trechos{j}{1},'{,}','.');
    t=strrep(t,'\!\times\!10^{','e'); t=strrep(t,'}','');
    if ~isempty(regexp(t,'^-?[0-9]+(\.[0-9]+)?(e-?[0-9]+)?$','once'))
        valores{end+1}=t; %#ok<AGROW>
    end
end
esperados=reshape(esperados.',[],1); quantidade=numel(valores);
assert(quantidade==numel(esperados),'Quantidade de valores nao confere: %s',nome);
for j=1:quantidade
    partes=strsplit(valores{j},'e');
    expoente=0; if numel(partes)>1, expoente=str2double(partes{2}); end
    v=str2double(valores{j});
    ponto=strfind(partes{1},'.');
    tolerancia=1e-13;
    if ~isempty(ponto), tolerancia=tolerancia+0.50001*10^(expoente-(numel(partes{1})-ponto)); end
    assert(abs(v-esperados(j))<=tolerancia, ...
        'Valor %d incorreto em %s: lido %.16g esperado %.16g tolerancia %.16g.', ...
        j,nome,v,esperados(j),tolerancia);
end
end
