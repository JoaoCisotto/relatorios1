function [coeficientes, energias] = T1_simulacoes_imagens(varargin)
% T1_SIMULACOES_IMAGENS Simulacoes e figuras produzidas diretamente no MATLAB.
%
% Basta abrir este arquivo no MATLAB e clicar em Run, ou executar:
%   addpath('scripts');
%   T1_simulacoes_imagens
%
% Exemplos:
%   T1_simulacoes_imagens('Modo','salvos')
%     Usa as amostras do Simulink ja salvas em validacao_completa.mat.
%
%   T1_simulacoes_imagens('Cenarios',[1 3 6 19 24],'MostrarFiguras',true)
%     Exporta apenas esses casos e mantem as figuras abertas para edicao.
%
%   T1_simulacoes_imagens('A',2,'N',2^17)
%     Refaz as simulacoes no MATLAB com outra amplitude e amostragem.
%
%   T1_simulacoes_imagens('GerarImagens',false)
%     Executa todos os calculos e exporta os dados, sem acessar graficos.
%
% O arquivo e independente de Simulink, LaTeX e dos outros scripts no modo
% 'matlab'. Nao e necessario instalar toolboxes para os calculos e as figuras.
% Desenvolvido no MATLAB R2026a; utiliza tiledlayout e exportgraphics.
%
% Cada caso exporta uma imagem com tres paineis:
%   1. Onda simulada e reconstrucoes teorica/FFT ate a ordem 7;
%   2. Amplitudes das ordens 1 a 7: teoria versus FFT;
%   3. Energias total, da fundamental e das ordens 1 a 7: teoria/simulacao.
% Sao contemplados todos os 24 cenarios do relatorio.

pastaScript = fileparts(mfilename('fullpath'));
[~,nomePasta] = fileparts(pastaScript);
if strcmpi(nomePasta,'scripts')
    raiz = fileparts(pastaScript);
else
    raiz = pastaScript;
end
p = inputParser;
positivo = @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0;
logico = @(x) islogical(x) && isscalar(x);
addParameter(p,'Modo','matlab',@(x) any(strcmpi(string(x),["matlab","salvos"])));
addParameter(p,'A',1,positivo);
addParameter(p,'P',2*pi,positivo);
addParameter(p,'Omega',120*pi,positivo);
addParameter(p,'N',2^16,@(x) positivo(x) && x>=32 && x==floor(x));
addParameter(p,'Cenarios',1:24,@(x) isnumeric(x) && isvector(x) && ...
    ~isempty(x) && all(isfinite(x)) && all(x==floor(x)) && all(x>=1 & x<=24));
addParameter(p,'Periodos',2,@(x) positivo(x) && x==floor(x));
addParameter(p,'GerarImagens',true,logico);
addParameter(p,'MostrarFiguras',false,logico);
addParameter(p,'SalvarFIG',true,logico);
addParameter(p,'Resolucao',300,positivo);
addParameter(p,'PastaSaida',fullfile(raiz,'resultados','imagens_matlab'), ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p,'ArquivoSalvo',fullfile(raiz,'resultados','validacao_completa.mat'), ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
parse(p,varargin{:}); op = p.Results;
selecionados = unique(op.Cenarios(:)','stable');
saida = char(op.PastaSaida);
if ~exist(saida,'dir'), mkdir(saida); end

if strcmpi(op.Modo,'salvos')
    assert(isfile(op.ArquivoSalvo), ...
        'Arquivo de amostras ausente. Execute no modo matlab ou indique ArquivoSalvo.');
    S = load(op.ArquivoSalvo,'dadosOndas','omega','temSimulink');
    assert(isfield(S,'dadosOndas') && numel(S.dadosOndas)==24, ...
        'O arquivo salvo deve conter as amostras dos 24 cenarios.');
    dados = S.dadosOndas; omega = S.omega;
    fonte = 'dados salvos do MATLAB';
    if isfield(S,'temSimulink') && S.temSimulink
        fonte = 'amostras salvas do Simulink';
    end
else
    omega = op.Omega; casos = definirCasos(op.A,op.P);
    dados = cell(24,1);
    for ic=1:24
        c=casos(ic); x=c.L*((0:op.N-1)'+0.5)/op.N;
        % A onda e gerada por trechos, independentemente da serie teorica.
        dados{ic}=struct('caso',c,'x',x,'y',simularOnda(x,c));
    end
    fonte = 'simulacao numerica no MATLAB';
end

coeficientes = table(); energias = table(); imagens = table();
dadosExportados = cell(numel(selecionados),1);
for j=1:numel(selecionados)
    ic=selecionados(j); q=dados{ic}; c=q.caso;
    x=q.x(:); y=q.y(:); N=numel(y); n=(1:7)';
    assert(N==numel(x) && N>=32,'Amostras invalidas no cenario %d.',ic);
    esperado=c.L*((0:N-1)'+0.5)/N;
    assert(max(abs(x-esperado))<1e-10*max(1,c.L), ...
        'As amostras devem formar um periodo completo em pontos medios.');

    [a0,a,b,Etotal] = coeficientesTeoricos(c);
    Y=fft(y)/N;
    cn=Y(2:8).*exp(-1i*pi*n/N); % correcao de meia amostra
    a0fft=real(Y(1)); afft=2*real(cn); bfft=-2*imag(cn);
    C=hypot(a,b); Cfft=hypot(afft,bfft);
    En=c.L/2*C.^2; Enfft=c.L/2*Cfft.^2;
    Etotalfft=c.L*mean(y.^2);
    Ecc=c.L*a0^2; Eccfft=c.L*a0fft^2;
    E7=sum(En); E7fft=sum(Enfft);
    eta=100*[En(1),E7,Ecc+E7]/Etotal;
    etafft=100*[Enfft(1),E7fft,Eccfft+E7fft]/Etotalfft;

    linhaCoef=table(repmat(ic,7,1),repmat(string(c.nome),7,1),n,a,b,afft,bfft, ...
        C,Cfft,En,Enfft,'VariableNames',{'Cenario','Caso','n','a_teoria', ...
        'b_teoria','a_fft','b_fft','C_teoria','C_fft','E_n_teoria','E_n_simulacao'});
    linhaEnergia=table(ic,string(c.nome),string(fonte),omega/c.L,a0,a0fft, ...
        Etotal,Etotalfft,En(1),Enfft(1),E7,E7fft, ...
        eta(1),etafft(1),eta(2),etafft(2),eta(3),etafft(3), ...
        'VariableNames',{'Cenario','Caso','Fonte','f0_Hz','a0_teoria','a0_simulacao', ...
        'E_total_teoria','E_total_simulacao','E1_teoria','E1_simulacao', ...
        'E1_7_teoria','E1_7_simulacao','eta1_teoria','eta1_simulacao', ...
        'eta1_7_teoria','eta1_7_simulacao','etaCC_7_teoria','etaCC_7_simulacao'});
    linhaCoef.erro_par_V=hypot(afft-a,bfft-b);
    linhaCoef.erro_C_V=abs(Cfft-C);
    linhaEnergia.erro_Etotal_pct=100*abs(Etotalfft-Etotal)/Etotal;
    linhaEnergia.erro_C1_pct=100*abs(Cfft(1)-C(1))/C(1);
    linhaEnergia.THD7_teoria_pct=100*sqrt(sum(C(2:7).^2))/C(1);
    linhaEnergia.THD7_simulacao_pct=100*sqrt(sum(Cfft(2:7).^2))/Cfft(1);
    linhaEnergia.THD_total_teoria_pct=100*sqrt(max(0,Etotal-Ecc-En(1))/En(1));
    linhaEnergia.THD_total_simulacao_pct=100*sqrt(max(0,Etotalfft-Eccfft-Enfft(1))/Enfft(1));
    coeficientes=[coeficientes;linhaCoef]; %#ok<AGROW>
    energias=[energias;linhaEnergia]; %#ok<AGROW>
    q.a0=a0; q.a=a; q.b=b; q.a0_fft=a0fft; q.a_fft=afft; q.b_fft=bfft;
    dadosExportados{j}=q;

    if op.GerarImagens
        base=sprintf('%02d_%s',ic,regexprep(c.nome,'[^a-zA-Z0-9_]','_'));
        fig=desenharFigura(q,omega,En,Enfft,Etotal,Etotalfft,eta,etafft,ic,fonte,op);
        try
            exportgraphics(fig,fullfile(saida,[base,'.png']),'Resolution',op.Resolucao);
            exportgraphics(fig,fullfile(saida,[base,'.pdf']),'ContentType','vector');
            if op.SalvarFIG, savefig(fig,fullfile(saida,[base,'.fig'])); end
        catch falha
            if ~op.MostrarFiguras && isgraphics(fig), close(fig); end
            rethrow(falha);
        end
        if ~op.MostrarFiguras, close(fig); end
        imagens=[imagens;table(ic,string(c.nome),string(base), ...
            'VariableNames',{'Cenario','Caso','Arquivo'})]; %#ok<AGROW>
    end
    fprintf('%02d/24 %-26s f0=%6.2f Hz; eta1 T/S=%8.4f / %8.4f %%\n', ...
        ic,c.nome,omega/c.L,eta(1),etafft(1));
end

% V^2 rad nas tabelas; dividir por omega para V^2 s.
writetable(coeficientes,fullfile(saida,'comparacao_coeficientes.csv'));
writetable(energias,fullfile(saida,'comparacao_energias.csv'));
save(fullfile(saida,'dados_simulacoes.mat'),'dadosExportados','coeficientes', ...
    'energias','omega','fonte','selecionados');
if op.GerarImagens
    writetable(imagens,fullfile(saida,'indice_imagens.csv'));
    escreverInclusoes(saida,imagens,fonte);
end
fprintf('\nConcluido: %d cenarios; %d pares de coeficientes.\n', ...
    numel(selecionados),height(coeficientes));
fprintf('Saida: %s\n',saida);
if ~op.GerarImagens
    fprintf('Somente calculos e dados exportados. GerarImagens=false.\n');
end
end

function casos=definirCasos(A,P)
casos=struct('nome',{},'tipo',{},'A',{},'P',{},'L',{},'phi',{}, ...
    'd',{},'phi2',{},'d2',{},'alternado',{});
for i=1:6
    casos(end+1)=caso(sprintf('Fig1_sinal%d',i),i,A,P,P,0,0,0,0,false);
end
for fator=[1,0.8]
    for i=1:6
        casos(end+1)=caso(sprintf('Fig2_%.1fP_sinal%d',fator,i), ...
            i,A,P,fator*P,P/10,0,0,0,false);
    end
end
for tipo=7:8
    casos(end+1)=caso(sprintf('Ret%d_simetrica',tipo),tipo,A,P,P,P/3,P/3,0,0,false);
    casos(end+1)=caso(sprintf('Ret%d_assimetrica',tipo),tipo,A,P,P,P/6,P/3,0,0,false);
    casos(end+1)=caso(sprintf('Ret%d_alternada',tipo),tipo,A,P,2*P,P/6,P/3,P/4,P/2,true);
end
end

function c=caso(nome,tipo,A,P,L,phi,d,phi2,d2,alternado)
c=struct('nome',nome,'tipo',tipo,'A',A,'P',P,'L',L,'phi',phi, ...
    'd',d,'phi2',phi2,'d2',d2,'alternado',alternado);
end

function y=simularOnda(x,c)
if c.tipo<=6
    u=mod(x-c.phi,c.L)/c.L;
    switch c.tipo
        case 1, y=c.A*u;
        case 2, y=c.A*(1-u);
        case 3, y=c.A*(1-2*abs(u-0.5));
        case 4, y=c.A*(2*u-1);
        case 5, y=c.A*(1-2*u);
        case 6, y=c.A*(1-4*abs(u-0.5));
    end
else
    z=mod(x,c.L); ligado=z>=c.phi & z<c.phi+c.d;
    if c.alternado
        ligado=ligado | (z>=c.P+c.phi2 & z<c.P+c.phi2+c.d2);
    end
    y=c.A*double(ligado);
    if c.tipo==8, y=2*y-c.A; end
end
end

function [a0,a,b,E]=coeficientesTeoricos(c)
n=(1:7)'; a=zeros(7,1); b=a;
if c.tipo<=6
    a0=c.A/2*double(c.tipo<=3);
    switch c.tipo
        case 1, b=-c.A./(pi*n);
        case 2, b=c.A./(pi*n);
        case 3, a=2*c.A*((-1).^n-1)./(pi^2*n.^2);
        case 4, b=-2*c.A./(pi*n);
        case 5, b=2*c.A./(pi*n);
        case 6, a=4*c.A*((-1).^n-1)./(pi^2*n.^2);
    end
    psi=2*pi*c.phi/c.L;
    deslocado=a.*cos(n*psi)-b.*sin(n*psi);
    b=a.*sin(n*psi)+b.*cos(n*psi); a=deslocado;
    E=c.L*c.A^2/3;
else
    k=2*pi*n/c.L; esq=c.phi; dir=c.phi+c.d;
    if c.alternado
        esq=[esq,c.P+c.phi2]; dir=[dir,c.P+c.phi2+c.d2];
    end
    a0=c.A*sum(dir-esq)/c.L;
    a=2*c.A./(c.L*k).*sum(sin(k*dir)-sin(k*esq),2);
    b=2*c.A./(c.L*k).*sum(cos(k*esq)-cos(k*dir),2);
    E=c.A^2*sum(dir-esq);
    if c.tipo==8, a0=2*a0-c.A; a=2*a; b=2*b; E=c.L*c.A^2; end
end
end

function fig=desenharFigura(q,omega,En,Enfft,Etotal,Etotalfft,eta,etafft,ic,fonte,op)
c=q.caso; N=numel(q.y);
visibilidade='off'; if op.MostrarFiguras, visibilidade='on'; end
fig=figure('Name',sprintf('Cenario %02d - %s',ic,c.nome), ...
    'Visible',visibilidade,'Color','w','Units','centimeters','Position',[2 2 9.5 20]);
layout=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
partes=strsplit(tituloCaso(ic),' - ');
cabecalho=[partes,{sprintf('A=%.3g V; f_0=%.3g Hz',c.A,omega/c.L),fonte}];
sgtitle(layout,cabecalho,'FontSize',10,'FontWeight','bold','Interpreter','tex');

% Seleciona pontos reais das amostras somente para a visualizacao.
indices=unique(round(linspace(1,N,min(N,1500))));
x=[]; y=[];
for periodo=0:op.Periodos-1
    x=[x;q.x(indices)+periodo*c.L]; y=[y;q.y(indices)]; %#ok<AGROW>
end
teoria=q.a0+zeros(size(x)); reconstruida=q.a0_fft+zeros(size(x));
for n=1:7
    teoria=teoria+q.a(n)*cos(2*pi*n*x/c.L)+q.b(n)*sin(2*pi*n*x/c.L);
    reconstruida=reconstruida+q.a_fft(n)*cos(2*pi*n*x/c.L)+q.b_fft(n)*sin(2*pi*n*x/c.L);
end
ax=nexttile(layout);
plot(ax,1000*x/omega,y,'.','Color',[0.25 0.25 0.25],'MarkerSize',3); hold(ax,'on');
plot(ax,1000*x/omega,teoria,'b-','LineWidth',1.25);
plot(ax,1000*x/omega,reconstruida,'r--','LineWidth',1.1);
grid(ax,'on'); box(ax,'on'); ylabel(ax,'Tensao (V)'); xlabel(ax,'Tempo (ms)');
title(ax,'Forma de onda e reconstrucoes ate a ordem 7');
legend(ax,'Amostras simuladas','Serie teorica','Serie obtida por FFT', ...
    'Location','best','FontSize',8);
xlim(ax,[0,1000*op.Periodos*c.L/omega]);

ax=nexttile(layout); f=(1:7)'*omega/c.L;
stem(ax,f,hypot(q.a,q.b),'b','filled','LineWidth',1); hold(ax,'on');
plot(ax,f,hypot(q.a_fft,q.b_fft),'ro','MarkerSize',6,'LineWidth',1.2);
grid(ax,'on'); box(ax,'on'); xticks(ax,f);
xlabel(ax,'Frequencia (Hz)'); ylabel(ax,'Amplitude de pico (V)');
title(ax,'FFT: fundamental e harmonicas de ordens 2 a 7');
legend(ax,'Teoria','FFT da simulacao','Location','best','FontSize',8);

ax=nexttile(layout);
valores=[Etotal,Etotalfft;En(1),Enfft(1);sum(En),sum(Enfft)]/omega;
bar(ax,valores,'grouped'); grid(ax,'on'); box(ax,'on');
xticks(ax,1:3); xticklabels(ax,{'Sinal completo','Fundamental','Ordens 1 a 7'});
ylabel(ax,'Energia quadratica (V^2 s)');
title(ax,{sprintf('\\eta_1 (T/S): %.4f%% / %.4f%%',eta(1),etafft(1)), ...
    sprintf('\\eta_{1:7} (T/S): %.4f%% / %.4f%%',eta(2),etafft(2))}, ...
    'Interpreter','tex','FontSize',8);
legend(ax,'Teoria','Simulacao','Location','best','FontSize',8);
set(findall(fig,'Type','axes'),'FontName','Arial','FontSize',9);
end

function t=tituloCaso(ic)
tipos={'Rampa crescente unipolar','Rampa decrescente unipolar','Triangular unipolar', ...
    'Rampa crescente bipolar','Rampa decrescente bipolar','Triangular bipolar'};
if ic<=18
    t=tipos{mod(ic-1,6)+1};
    if ic>6 && ic<=12, t=[t,' - \delta=P, \phi=P/10']; end
    if ic>12, t=[t,' - \delta=0.8P, \phi=P/10']; end
else
    tipos={'Retangular unipolar simetrica','Retangular unipolar assimetrica', ...
        'Retangular unipolar alternada','Retangular bipolar simetrica', ...
        'Retangular bipolar assimetrica','Retangular bipolar alternada'};
    t=tipos{ic-18};
end
t=sprintf('Caso %02d: %s',ic,t);
end

function escreverInclusoes(saida,imagens,fonte)
fid=fopen(fullfile(saida,'figuras_para_relatorio.tex'),'w');
limpa=onCleanup(@()fclose(fid));
fprintf(fid,'%% Copie os arquivos PDF para a pasta images/matlab do seu relatorio.\n');
fprintf(fid,'%% Figuras geradas a partir de %s.\n',fonte);
for j=1:height(imagens)
    base=char(imagens.Arquivo(j)); ic=imagens.Cenario(j);
    fprintf(fid,'\\begin{figure}[htbp]\n\\centering\n');
    fprintf(fid,'\\includegraphics[width=\\linewidth]{images/matlab/%s.pdf}\n',base);
    caption=strrep(tituloCaso(ic),'\delta=P','$\delta=P$');
    caption=strrep(caption,'\delta=0.8P','$\delta=0{,}8P$');
    caption=strrep(caption,'\phi=P/10','$\phi=P/10$');
    fprintf(fid,'\\caption{%s. Comparacao entre teoria e simulacao: onda, FFT e energias.}\n',caption);
    fprintf(fid,'\\label{fig:matlab_caso_%02d}\n\\end{figure}\n\n',ic);
end
end
