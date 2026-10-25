function T1_Fourier_MATLAB(modo)
% Reproduz os 24 cenarios, gera dados da FFT e executa o Simulink.
% Executar na pasta do projeto: addpath('scripts'); T1_Fourier_MATLAB
raiz = fileparts(fileparts(mfilename('fullpath')));
saida = fullfile(raiz,'resultados');
if ~exist(saida,'dir'), mkdir(saida); end
diary(fullfile(saida,'execucao_matlab.txt'));
if nargin>0 && strcmp(modo,'artefatos')
    % Recria somente tabelas e graficos a partir de uma execucao ja salva.
    S=load(fullfile(saida,'validacao_completa.mat'));
    assert(isfield(S,'dadosOndas'),'Execute o modo completo para salvar as amostras simuladas.');
    T1_gerar_apresentacao(saida,S.coef,S.energias,S.resumo,S.dadosOndas,S.temSimulink);
    registrarAmbiente(saida,S.N,S.omega,S.temSimulink);
    fprintf('Artefatos recriados a partir dos 24 cenarios salvos.\n'); diary off; return
end
A = 1; P = 2*pi; omega = 120*pi; N = 2^16;
casos = struct('nome',{},'tipo',{},'A',{},'P',{},'L',{},'phi',{}, ...
    'd',{},'phi2',{},'d2',{},'alternado',{});
for i = 1:6
    casos(end+1) = caso(sprintf('Fig1_sinal%d',i),i,A,P,P,0,0,0,0,false);
end
for fator = [1,0.8]
    for i = 1:6
        casos(end+1) = caso(sprintf('Fig2_%.1fP_sinal%d',fator,i), ...
            i,A,P,fator*P,P/10,0,0,0,false);
    end
end
for tipo = 7:8
    casos(end+1) = caso(sprintf('Ret%d_simetrica',tipo),tipo,A,P,P,P/3,P/3,0,0,false);
    casos(end+1) = caso(sprintf('Ret%d_assimetrica',tipo),tipo,A,P,P,P/6,P/3,0,0,false);
    casos(end+1) = caso(sprintf('Ret%d_alternada',tipo),tipo,A,P,2*P,P/6,P/3,P/4,P/2,true);
end

temSimulink = license('test','Simulink') && ~isempty(ver('simulink'));
modelo = 'T1_validacao_sinais';
if temSimulink
    Simulink.fileGenControl('set','CacheFolder',fullfile(tempdir,'T1_simulink_cache'), ...
        'CodeGenFolder',fullfile(tempdir,'T1_simulink_codegen'),'createDir',true);
    load_system('simulink');
    if bdIsLoaded(modelo), close_system(modelo,0); end
    new_system(modelo);
    add_block('simulink/Sources/Digital Clock',[modelo '/Tempo'], ...
        'Position',[35 55 85 85]);
    add_block('simulink/User-Defined Functions/MATLAB Function',[modelo '/Onda'], ...
        'Position',[145 40 280 100]);
    add_block('simulink/Sinks/To Workspace',[modelo '/Amostras'], ...
        'VariableName','onda_simulada','SaveFormat','Timeseries', ...
        'MaxDataPoints','inf','Position',[345 50 455 90]);
    add_line(modelo,'Tempo/1','Onda/1');
    add_line(modelo,'Onda/1','Amostras/1');
    set_param(modelo,'SolverType','Fixed-step','Solver','FixedStepDiscrete', ...
        'ReturnWorkspaceOutputs','on','SaveTime','off','SaveOutput','off');
    arvore = sfroot;
    bloco = arvore.find('-isa','Stateflow.EMChart','Path',[modelo '/Onda']);
end

coef = []; energias = []; resumo = []; dadosOndas = cell(numel(casos),1);
for ic = 1:numel(casos)
    c = casos(ic); Ts = c.L/(omega*N);
    x = c.L*((0:N-1)'+0.5)/N;
    y = onda(x,c);
    [a0,a,b,E] = teoria(c);
    ysim = y;
    if temSimulink
        set_param([modelo '/Tempo'],'SampleTime',sprintf('%.17g',Ts));
        set_param(modelo,'FixedStep',sprintf('%.17g',Ts), ...
            'StopTime',sprintf('%.17g',(N-1)*Ts));
        bloco.Script = codigoBloco(c,omega,Ts);
        simOut = sim(modelo);
        dados = simOut.get('onda_simulada'); ys = dados.Data(:);
        assert(numel(ys)==N,'Numero de amostras diferente do previsto.');
        ysim = ys;
    end
    % A FFT e as energias numericas usam diretamente as amostras simuladas.
    Y = fft(ysim)/N;
    cn = Y(2:8).*exp(-1i*pi*(1:7)'/N);
    af = 2*real(cn); bf = -2*imag(cn); a0f = real(Y(1));
    Pnum = mean(ysim.^2);
    Ecc = c.L*a0^2; E1 = c.L/2*(a(1)^2+b(1)^2);
    E7 = c.L/2*sum(a.^2+b.^2); Ecc7 = Ecc+E7;
    Esim = c.L*Pnum; Eccsim = c.L*a0f^2;
    E1sim = c.L/2*(af(1)^2+bf(1)^2);
    E7sim = c.L/2*sum(af.^2+bf.^2); Ecc7sim = Eccsim+E7sim;
    coef = [coef; table(repmat(string(c.nome),7,1),(1:7)',a,b,af,bf, ...
        hypot(a,b),hypot(af,bf), ...
        'VariableNames',{'Caso','n','a_teo','b_teo','a_fft','b_fft', ...
        'C_teo','C_fft'})]; %#ok<AGROW>
    energias = [energias; table(string(c.nome),c.L,omega/c.L,a0,E,Ecc,E1,E7,Ecc7, ...
        100*E1/E,100*E7/E,100*Ecc7/E,E/omega, ...
        Esim,Eccsim,E1sim,E7sim,Ecc7sim,100*E1sim/Esim, ...
        100*E7sim/Esim,100*Ecc7sim/Esim,Esim/omega, ...
        'VariableNames',{'Caso','L','f0_Hz','a0','E_total_angular','E_CC', ...
        'E1','E1_7','E_CC_7','eta1','eta1_7','etaCC_7','E_total_temporal', ...
        'E_total_sim','E_CC_sim','E1_sim','E1_7_sim','E_CC_7_sim', ...
        'eta1_sim','eta1_7_sim','etaCC_7_sim','E_total_temporal_sim'})]; %#ok<AGROW>
    resumo = [resumo; table(string(c.nome),a0,a0f,hypot(a(1),b(1)),hypot(af(1),bf(1)),E/c.L,Pnum, ...
        'VariableNames',{'Caso','a0_teo','a0_sim','C1_teo','C1_sim','Pq_teo','Pq_sim'})]; %#ok<AGROW>
    dadosOndas{ic} = struct('caso',c,'x',x,'y',ysim,'a0',a0,'a',a,'b',b, ...
        'a0_fft',a0f,'a_fft',af,'b_fft',bf);
    assert(all(abs([a0f-a0;af-a;bf-b])<2e-4),'Coeficientes fora da tolerancia de amostragem.');
    assert(abs(Pnum-E/c.L)<2e-4,'Valor quadratico fora da tolerancia de amostragem.');
    fprintf('%2d/24 %s: a0 teoria/sim = %.6f / %.6f; C1 = %.6f / %.6f\n', ...
        ic,c.nome,a0,a0f,hypot(a(1),b(1)),hypot(af(1),bf(1)));
end
assert(all(energias.etaCC_7<=100+1e-10),'Energia truncada superior a energia total.');
writetable(coef,fullfile(saida,'coeficientes_fft.csv'));
writetable(energias,fullfile(saida,'energias.csv'));
writetable(resumo,fullfile(saida,'validacao_fft.csv'));
save(fullfile(saida,'validacao_completa.mat'),'casos','coef','energias','resumo','dadosOndas','N','omega','temSimulink');
if temSimulink
    save_system(modelo,fullfile(saida,[modelo '.slx']));
    close_system(modelo,0);
end
T1_gerar_apresentacao(saida,coef,energias,resumo,dadosOndas,temSimulink);
registrarAmbiente(saida,N,omega,temSimulink);
fprintf('Concluido: %d cenarios. Simulink executado: %d.\n',numel(casos),temSimulink);
diary off;
end

function registrarAmbiente(saida,N,omega,temSimulink)
fid = fopen(fullfile(saida,'ambiente.txt'),'w');
fprintf(fid,'MATLAB: %s\nN = %d\nomega = %.17g rad/s\nSimulink executado: %d\n', ...
    version,N,omega,temSimulink);
fprintf(fid,'Geracao Simulink: Digital Clock -> MATLAB Function -> To Workspace.\n');
fprintf(fid,'Amostras no ponto medio; FFT corrigida por exp(-i*pi*n/N).\n');
fclose(fid);
end

function c = caso(nome,tipo,A,P,L,phi,d,phi2,d2,alternado)
c = struct('nome',nome,'tipo',tipo,'A',A,'P',P,'L',L,'phi',phi, ...
    'd',d,'phi2',phi2,'d2',d2,'alternado',alternado);
end

function y = onda(x,c)
if c.tipo<=6
    u = mod(x-c.phi,c.L)/c.L;
    switch c.tipo
        case 1, y = c.A*u;
        case 2, y = c.A*(1-u);
        case 3, y = c.A*(1-2*abs(u-0.5));
        case 4, y = c.A*(2*u-1);
        case 5, y = c.A*(1-2*u);
        case 6, y = c.A*(1-4*abs(u-0.5));
    end
else
    z = mod(x,c.L);
    ligado = z>=c.phi & z<c.phi+c.d;
    if c.alternado
        ligado = ligado | (z>=c.P+c.phi2 & z<c.P+c.phi2+c.d2);
    end
    y = c.A*double(ligado);
    if c.tipo==8, y = 2*y-c.A; end
end
end

function [a0,a,b,E] = teoria(c)
n = (1:7)'; a = zeros(7,1); b = a;
if c.tipo<=6
    a0 = c.A/2*double(c.tipo<=3);
    switch c.tipo
        case 1, b = -c.A./(pi*n);
        case 2, b = c.A./(pi*n);
        case 3, a = 2*c.A*((-1).^n-1)./(pi^2*n.^2);
        case 4, b = -2*c.A./(pi*n);
        case 5, b = 2*c.A./(pi*n);
        case 6, a = 4*c.A*((-1).^n-1)./(pi^2*n.^2);
    end
    psi = 2*pi*c.phi/c.L;
    aa = a.*cos(n*psi)-b.*sin(n*psi);
    b = a.*sin(n*psi)+b.*cos(n*psi); a = aa;
    E = c.L*c.A^2/3;
else
    k = 2*pi*n/c.L;
    esq = c.phi; dir = c.phi+c.d;
    if c.alternado
        esq = [esq,c.P+c.phi2]; dir = [dir,c.P+c.phi2+c.d2];
    end
    a0 = c.A*sum(dir-esq)/c.L;
    a = 2*c.A./(c.L*k).*sum(sin(k*dir)-sin(k*esq),2);
    b = 2*c.A./(c.L*k).*sum(cos(k*esq)-cos(k*dir),2);
    E = c.A^2*sum(dir-esq);
    if c.tipo==8
        a0 = 2*a0-c.A; a = 2*a; b = 2*b; E = c.L*c.A^2;
    end
end
end

function s = codigoBloco(c,omega,Ts)
% Equacoes por trechos; nao usa coeficientes teoricos nem sintese harmonica.
s = sprintf(['function y = f(t)\n' ...
    'x=%.17g*(t+%.17g/2);\nA=%.17g; L=%.17g; phi=%.17g;\n'],omega,Ts,c.A,c.L,c.phi);
if c.tipo<=6
    expressoes = {'A*u','A*(1-u)','A*(1-2*abs(u-0.5))', ...
        'A*(2*u-1)','A*(1-2*u)','A*(1-4*abs(u-0.5))'};
    s = [s,sprintf('u=mod(x-phi,L)/L;\ny=%s;\nend\n',expressoes{c.tipo})];
else
    s = [s,sprintf('z=mod(x,L);\nligado=(z>=phi)&&(z<phi+%.17g);\n',c.d)];
    if c.alternado
        s = [s,sprintf('ligado=ligado||((z>=%.17g)&&(z<%.17g));\n', ...
            c.P+c.phi2,c.P+c.phi2+c.d2)];
    end
    s = [s,sprintf('y=A*double(ligado);\n')];
    if c.tipo==8, s = [s,sprintf('y=2*y-A;\n')]; end
    s = [s,sprintf('end\n')];
end
end

