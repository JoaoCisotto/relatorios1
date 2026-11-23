function res = plotar_3etapas(out, fs, opc)
% Equivalência em 3 etapas a partir da saída 'out' da simulação (Simulink).
%   Etapa 1: out.x_circuito  (circuito)
%   Etapa 2: out.x_chaveado  (equações chaveadas, entrada = PWM)
%   Etapa 3: out.x_medio     (modelo médio, entrada = d_medio)
% Colunas: [ib ic ipv vdc vc]
%
% Uso:   out = sim('SEU_MODELO');
%        res = plotar_3etapas(out, 20e3);
%        res = plotar_3etapas(out, 20e3, struct('t_degrau',0.1,'t_zoom',0.05));
%
% Opções (campos de opc, todos opcionais):
%   arquivo  nome base dos arquivos          ('fig_equiv')
%   pasta    pasta de saída                  ('images')
%   largura, altura   tamanho em cm          (8.8, 11  -> coluna IEEE)
%   modo     'auto' | 'bruto' | 'envelope'   ('auto': envelope se a simulação > 0,1 s)
%   t_degrau instante de um degrau [s], desenha linha vertical   ([])
%   t_zoom   instante inicial [s] da figura de 3 períodos (fig_equiv_zoom). []: não gera
%   nsp      amostras por período na grade uniforme (100)
% Salva PDF vetorial (e PNG) e imprime os erros para a Tabela do relatório.
if nargin < 1 || isempty(out), out = evalin('base','out'); end
if nargin < 2 || isempty(fs),  fs  = 20e3; end
if nargin < 3, opc = struct; end
opc = padrao(opc);
Ts = 1/fs;  nsp = opc.nsp;  dt = Ts/nsp;

[t1,X1] = ler_sinal(out,'x_circuito',5);
[t2,X2] = ler_sinal(out,'x_chaveado',5);
[t3,X3] = ler_sinal(out,'x_medio',5);
T  = min([t1(end) t2(end) t3(end)]);
np = floor(T/Ts + 1e-9);
tu = (0:dt:np*Ts)';
x1 = interp1(t1,X1,tu,'linear','extrap');  clear X1 t1
x2 = interp1(t2,X2,tu,'linear','extrap');  clear X2 t2
x3 = interp1(t3,X3,tu,'linear','extrap');  clear X3 t3
n = np*nsp;
bl = @(x) reshape(x(1:n,:), nsp, np, 5);                  % blocos de 1 período

% médias por período e envelopes (mín/máx) de cada período
xm  = reshape(mean(bl(x1),1), np, 5);
mx1 = reshape(max(bl(x1),[],1), np, 5);  mn1 = reshape(min(bl(x1),[],1), np, 5);
mx2 = reshape(max(bl(x2),[],1), np, 5);  mn2 = reshape(min(bl(x2),[],1), np, 5);
tm  = ((1:np)' - 0.5)*Ts;
x3m = interp1(tu, x3, tm);

% erros (o 1º período é ignorado na etapa 3: d_medio ainda está se formando)
res.max12 = max(abs(x1 - x2), [], 1);
e3 = x3m(2:end,:) - xm(2:end,:);
res.rms3 = sqrt(mean(e3.^2,1));   res.max3 = max(abs(e3),[],1);
nomes = {'ib [A]','ic [A]','ipv [A]','vdc [V]','vc [V]'};
fprintf('\n%-9s %-16s %-16s %-16s\n','Grandeza','Max etapas 1-2','RMS media-modelo','Max media-modelo');
for j = 1:5
    fprintf('%-9s %-16.3g %-16.3g %-16.3g\n',nomes{j},res.max12(j),res.rms3(j),res.max3(j));
end
try   % diagnóstico: d_medio deve ficar em [Db Dc Dpv]
    [td,Dm] = ler_sinal(out,'d_medio',3);  i = td >= td(end) - Ts;
    fprintf('\nd_medio (ultimo periodo): [%.4f %.4f %.4f]\n', mean(Dm(i,:),1));
catch
end

% ---------------- cores ----------------
cAz = [0.00 0.35 0.75];  cLar = [0.95 0.55 0.10];  cLarC = [1.00 0.87 0.70];
cVd = [0.10 0.55 0.20];  cVm  = [0.80 0.10 0.10];

% ---------------- figura principal ----------------
usa_env = strcmp(opc.modo,'envelope') || (strcmp(opc.modo,'auto') && T > 0.1);
ids = [1 2 4];  rot = {'$i_b$ [A]','$i_c$ [A]','$v_{dc}$ [V]'};
f = figure('Color','w','Units','centimeters'); f.Position = [2 2 opc.largura opc.altura];
ax = gobjects(3,1);
for q = 1:3
    ax(q) = subplot(3,1,q); hold on; j = ids(q);
    if usa_env
        % etapa 2 ao fundo (faixa laranja clara), depois circuito, média e modelo por cima
        h2 = fill([tm;flipud(tm)]*1e3, [mx2(:,j);flipud(mn2(:,j))], cLarC, 'EdgeColor',cLar, 'LineWidth',0.3);
        h1 = plot(tm*1e3, mx1(:,j), 'Color',cAz, 'LineWidth',0.4);
        plot(tm*1e3, mn1(:,j), 'Color',cAz, 'LineWidth',0.4, 'HandleVisibility','off');
    else
        d = 2;  k = 1:d:numel(tu);
        h2 = plot(tu(k)*1e3, x2(k,j), 'Color',cLar, 'LineWidth',0.9);
        h1 = plot(tu(k)*1e3, x1(k,j), 'Color',cAz,  'LineWidth',0.5);
    end
    hm = plot(tm*1e3, xm(:,j), 'Color',cVd, 'LineWidth',1.2);
    h3 = plot(tu(1:max(1,round(nsp/10)):end)*1e3, x3(1:max(1,round(nsp/10)):end,j), '-.', 'Color',cVm, 'LineWidth',1.0);
    if q == 1, Hleg = [h1 h2 hm h3]; end
    ylabel(rot{q}, 'Interpreter','latex');
    xlim([0 np*Ts*1e3]);
    if ~isempty(opc.t_degrau)
        yl = ylim; plot([1 1]*opc.t_degrau*1e3, yl, 'k:', 'LineWidth',0.6, 'HandleVisibility','off'); ylim(yl);
    end
end
xlabel('Tempo [ms]', 'Interpreter','latex');
formatar(ax);
lg = legend(ax(1), Hleg, {'Etapa 1: circuito','Etapa 2: equações','Média por período','Etapa 3: modelo médio'}, ...
            'Location','northoutside','Orientation','horizontal','NumColumns',2,'FontSize',6.5);
lg.Box = 'off';
salvar(f, opc.arquivo, opc.pasta);

% ---------------- figura de detalhe (3 períodos) ----------------
if ~isempty(opc.t_zoom)
    t0 = floor(opc.t_zoom/Ts)*Ts;  m = tu >= t0 - 1e-12 & tu <= t0 + 3*Ts + 1e-12;
    tz = (tu(m) - t0)*1e3;
    g = figure('Color','w','Units','centimeters'); g.Position = [2 2 opc.largura 9];
    az = gobjects(3,1);
    for q = 1:3
        az(q) = subplot(3,1,q); hold on; j = ids(q);
        hb = plot(tz, x2(m,j), 'Color',cLar, 'LineWidth',1.4);          % etapa 2 ao fundo
        hf = plot(tz, x1(m,j), 'Color',cAz,  'LineWidth',0.7);          % circuito por cima
        if q == 1, Hz = [hf hb]; end
        ylabel(rot{q}, 'Interpreter','latex'); xlim([0 3*Ts*1e3]);
    end
    xlabel('Tempo [ms]', 'Interpreter','latex');
    formatar(az);
    lz = legend(az(1), Hz, {'Etapa 1: circuito','Etapa 2: equações'}, ...
                'Location','northoutside','Orientation','horizontal','FontSize',6.5);
    lz.Box = 'off';
    salvar(g, [opc.arquivo '_zoom'], opc.pasta);
end
end

% ======================= funções locais =======================
function opc = padrao(opc)
d = struct('arquivo','fig_equiv','pasta','images','largura',8.8,'altura',11, ...
           'modo','auto','t_degrau',[],'t_zoom',[],'nsp',100);
nm = fieldnames(d);
for i = 1:numel(nm)
    if ~isfield(opc,nm{i}) || isempty(opc.(nm{i})) && ~ismember(nm{i},{'t_degrau','t_zoom'})
        opc.(nm{i}) = d.(nm{i});
    end
end
end

function formatar(ax)
for i = 1:numel(ax)
    set(ax(i), 'FontSize',7, 'TickLabelInterpreter','latex', 'Box','on', ...
        'XGrid','on','YGrid','on', 'GridAlpha',0.25, 'LineWidth',0.5, 'Layer','top');
    if i < numel(ax), set(ax(i),'XTickLabel',[]); end
end
end

function salvar(f, nome, pasta)
if ~exist(pasta,'dir'), mkdir(pasta); end
arq = fullfile(pasta,[nome '.pdf']);
try
    exportgraphics(f, arq, 'ContentType','vector');
catch                                   % MATLAB anterior à R2020a
    f.PaperUnits = 'centimeters'; f.PaperSize = f.Position(3:4); f.PaperPosition = [0 0 f.Position(3:4)];
    print(f, arq, '-dpdf', '-painters');
end
try, exportgraphics(f, fullfile(pasta,[nome '.png']), 'Resolution',300); catch, end
fprintf('Figura salva em %s\n', arq);
end

function [t,X] = ler_sinal(out,nome,ncol)
nomes = out.getElementNames;
if ~ismember(nome,nomes)
    error('A saída não tem ''%s''. Disponíveis: %s', nome, strjoin(nomes,', '));
end
S = out.get(nome);
if isa(S,'timeseries'),  t = S.Time(:); X = S.Data;
elseif isstruct(S),      t = S.time(:); X = S.signals.values;
else,                    t = out.tout(:); X = S;
end
X = squeeze(X);
if size(X,1) ~= numel(t) && size(X,2) == numel(t), X = X.'; end
if size(X,1) ~= numel(t)
    error('''%s'': %d amostras de dados contra %d instantes de tempo.', nome, size(X,1), numel(t));
end
if size(X,2) ~= ncol
    error('''%s'' deveria ter %d colunas, tem %d.', nome, ncol, size(X,2));
end
[t,iu] = unique(t);  X = X(iu,:);       % remove instantes repetidos (passo variável)
end
