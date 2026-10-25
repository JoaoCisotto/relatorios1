function T1_gerar_apresentacao(pasta,coef,E,R,dados,temSimulink)
% Tabelas teoria/simulacao e imagens vetoriais dos dados calculados no MATLAB.
% A renderizacao em LaTeX dispensa o motor grafico interativo do MATLAB.
fid=fopen(fullfile(pasta,'resumo_validacao.tex'),'w');
fprintf(fid,['Foram simulados 24 cenários, com $N=65536$ amostras por período real. ' ...
    'As tabelas apresentam lado a lado os valores da teoria (T) e da simulação (S). ']);
if temSimulink
    fprintf(fid,'Os valores S foram obtidos aplicando a FFT do MATLAB às amostras geradas no Simulink.\n');
else
    fprintf(fid,'Os valores S foram obtidos aplicando a FFT às ondas geradas numericamente no MATLAB.\n');
end
fclose(fid);

% Media e fundamental de todos os 24 cenarios.
fid=fopen(fullfile(pasta,'tabela_validacao.tex'),'w');
grupos={1:6,7:12,13:18,19:24};
legendas={'Ondas básicas: média e amplitude da fundamental em V.', ...
    'Ondas deslocadas com $\delta=P$: média e fundamental em V.', ...
    'Ondas deslocadas com $\delta=0{,}8P$: média e fundamental em V.', ...
    'Retangulares: média e amplitude da fundamental em V.'};
labels={'tab:fft','tab:fft_P','tab:fft_08P','tab:fft_ret'};
for j=1:4
    inicioTabela(fid,legendas{j},labels{j},'lrrrr');
    fprintf(fid,'Caso & $a_{0,\\mathrm T}$ & $a_{0,\\mathrm S}$ & $C_{1,\\mathrm T}$ & $C_{1,\\mathrm S}$ \\\\ \\hline\n');
    for r=grupos{j}
        fprintf(fid,'%s & %s & %s & %s & %s \\\\\n',rotulo(r), ...
            decimal(R.a0_teo(r),6),decimal(R.a0_sim(r),6), ...
            decimal(R.C1_teo(r),6),decimal(R.C1_sim(r),6));
    end
    fimTabela(fid);
end
fclose(fid);

% Todos os sete pares de coeficientes de todos os 24 cenarios.
fid=fopen(fullfile(pasta,'tabela_coeficientes_comparacao.tex'),'w');
for j=1:height(R)
    f=coef(coef.Caso==R.Caso(j),:);
    assert(height(f)==7 && isequal(f.n,(1:7)'),'Cobertura harmonica incompleta.');
    legenda=[nomeGrafico(j),': coeficientes em V, teoria (T) e FFT da simulação (S).'];
    inicioTabela(fid,legenda,sprintf('tab:coef_comp%d',j),'crrrr');
    fprintf(fid,'$n$ & $a_{n,\\mathrm T}$ & $a_{n,\\mathrm S}$ & $b_{n,\\mathrm T}$ & $b_{n,\\mathrm S}$ \\\\ \\hline\n');
    for n=1:7
        fprintf(fid,'%d & %s & %s & %s & %s \\\\\n',n, ...
            decimal(f.a_teo(n),6),decimal(f.a_fft(n),6), ...
            decimal(f.b_teo(n),6),decimal(f.b_fft(n),6));
    end
    fimTabela(fid);
end
fclose(fid);

% Energia absoluta e participacao percentual: resultados teoricos E e
% estimativas da simulacao calculadas das amostras e dos coeficientes FFT.
fid=fopen(fullfile(pasta,'tabela_energias.tex'),'w');
idsEnergia={1:6,7:12,13:18,19:24};
for grupo=1:4
    textos={'ondas básicas','ondas deslocadas com $\delta=P$', ...
        'ondas deslocadas com $\delta=0{,}8P$','retangulares'};
    inicioTabela(fid,sprintf('Energia angular em V$^2$rad: %s; T = teoria, S = simulação.',textos{grupo}), ...
        sprintf('tab:energia_abs%d',grupo),'llrrr');
    fprintf(fid,'Caso & Origem & $E_L^{(x)}$ & $E_1$ & $E_{1:7}$ \\\\ \\hline\n');
    for r=idsEnergia{grupo}
        fprintf(fid,'%s & T & %s & %s & %s \\\\\n',rotulo(r), ...
            decimal(E.E_total_angular(r),6),decimal(E.E1(r),6),decimal(E.E1_7(r),6));
        fprintf(fid,' & S & %s & %s & %s \\\\\n',decimal(E.E_total_sim(r),6), ...
            decimal(E.E1_sim(r),6),decimal(E.E1_7_sim(r),6));
    end
    fimTabela(fid);
    labelsEnergia={'tab:energias','tab:energias_P','tab:energias_08P','tab:energias_ret'};
    lab=labelsEnergia{grupo};
    inicioTabela(fid,sprintf('Percentuais de energia (\\%%): %s, teoria e simulação.',textos{grupo}),lab,'llrrr');
    fprintf(fid,'Caso & Origem & $\\eta_1$ & $\\eta_{1:7}$ & $\\eta_{\\mathrm{CC}:7}$ \\\\ \\hline\n');
    for r=idsEnergia{grupo}
        fprintf(fid,'%s & T & %s & %s & %s \\\\\n',rotulo(r), ...
            decimal(E.eta1(r),4),decimal(E.eta1_7(r),4),decimal(E.etaCC_7(r),4));
        fprintf(fid,' & S & %s & %s & %s \\\\\n',decimal(E.eta1_sim(r),4), ...
            decimal(E.eta1_7_sim(r),4),decimal(E.etaCC_7_sim(r),4));
    end
    fimTabela(fid);
end
fclose(fid);

% Tabelas exclusivamente teoricas usadas no desenvolvimento.
for j=1:3
    idsRet=[19,20,21]; r=idsRet(j); f=coef(coef.Caso==R.Caso(r),:);
    fid=fopen(fullfile(pasta,sprintf('retangular_coef_%d.tex',j)),'w');
    titulos={'simétrico','assimétrico','com dois ciclos alternados'};
    inicioTabela(fid,sprintf('Coeficientes teóricos unipolares do exemplo %s, divididos por $A$.',titulos{j}), ...
        sprintf('tab:ret%d',j),'crr');
    fprintf(fid,'$n$ & $a_n/A$ & $b_n/A$ \\\\ \\hline\n');
    for n=1:7
        fprintf(fid,'%d & %s & %s \\\\\n',n,decimal(f.a_teo(n)/dados{r}.caso.A,6), ...
            decimal(f.b_teo(n)/dados{r}.caso.A,6));
    end
    fimTabela(fid); fclose(fid);
end
gerarImagens(pasta,dados);
end

function inicioTabela(fid,caption,label,cols)
fprintf(fid,'\\begin{table}[H]\n\\centering\\small\n\\setlength{\\tabcolsep}{3.5pt}\n');
fprintf(fid,'\\caption{%s}\\label{%s}\n\\begin{tabular}{%s}\\hline\n',caption,label,cols);
end

function fimTabela(fid)
fprintf(fid,'\\hline\\end{tabular}\n\\end{table}\n');
end

function t=decimal(x,n)
if abs(x)<0.5*10^(-n), x=0; end
t=['$',strrep(sprintf(['%.',num2str(n),'f'],x),'.','{,}'),'$'];
end

function t=rotulo(r)
if r<=18
    t=sprintf('S%d',mod(r-1,6)+1);
else
    nomes={'Uni. sim.','Uni. assim.','Uni. alt.','Bi. sim.','Bi. assim.','Bi. alt.'};
    t=nomes{r-18};
end
end

function gerarImagens(pasta,dados)
raiz=fileparts(pasta); antiga=pwd; restaura=onCleanup(@()cd(antiga)); cd(raiz);
imagens=fullfile(pasta,'figuras'); if ~exist(imagens,'dir'), mkdir(imagens); end
grupos={1:3,4:6,7:9,10:12,13:15,16:18,19:21,22:24};
arquivos={'comparacao_basicas_unipolares','comparacao_basicas_bipolares', ...
    'comparacao_deslocadas_P_unipolares','comparacao_deslocadas_P_bipolares', ...
    'comparacao_deslocadas_08P_unipolares','comparacao_deslocadas_08P_bipolares', ...
    'comparacao_retangulares_unipolares','comparacao_retangulares_bipolares'};
captions={'Sinais 1 a 3 da Figura 1: ondas unipolares.', ...
    'Sinais 4 a 6 da Figura 1: ondas bipolares.', ...
    'Ondas unipolares deslocadas, $\delta=P$.', ...
    'Ondas bipolares deslocadas, $\delta=P$.', ...
    'Ondas unipolares deslocadas, $\delta=0{,}8P$.', ...
    'Ondas bipolares deslocadas, $\delta=0{,}8P$.', ...
    'Retangulares unipolares: casos simétrico, assimétrico e alternado.', ...
    'Retangulares bipolares: casos simétrico, assimétrico e alternado.'};
report=fopen(fullfile(pasta,'graficos_fourier.tex'),'w');
cobertos=[];
for g=1:numel(grupos)
    ficheiro=fullfile(imagens,[arquivos{g},'.tex']); fid=fopen(ficheiro,'w');
    inicioImagem(fid);
    for r=grupos{g}
        cobertos(end+1)=r; %#ok<AGROW>
        q=dados{r}; c=q.caso; N=numel(q.x);
        idx=unique([1:128:N,N]); u=q.x(idx)/c.L; simulada=q.y(idx);
        teoria7=q.a0+zeros(size(u)); fft7=q.a0_fft+zeros(size(u));
        for n=1:7
            teoria7=teoria7+q.a(n)*cos(2*pi*n*u)+q.b(n)*sin(2*pi*n*u);
            fft7=fft7+q.a_fft(n)*cos(2*pi*n*u)+q.b_fft(n)*sin(2*pi*n*u);
        end
        csv=sprintf('onda_%02d.csv',r);
        writetable(table(u,simulada,teoria7,fft7),fullfile(pasta,csv));
        fprintf(fid,'\\begin{tikzpicture}\n\\begin{axis}[width=8cm,height=3.6cm,grid=major,');
        fprintf(fid,'xmin=0,xmax=1,xlabel={$x/L$},ylabel={$v$ (V)},title={%s},font=\\scriptsize]\n',nomeGrafico(r));
        fprintf(fid,'\\addplot[black,only marks,mark=*,mark size=0.35pt] table[x=u,y=simulada,col sep=comma]{resultados/%s};\n',csv);
        fprintf(fid,'\\addplot[blue,thin] table[x=u,y=teoria7,col sep=comma]{resultados/%s};\n',csv);
        fprintf(fid,'\\addplot[red,dashed,thin] table[x=u,y=fft7,col sep=comma]{resultados/%s};\n',csv);
        fprintf(fid,'\\end{axis}\n\\end{tikzpicture}\\\\\n');

        n=(1:7)'; teoriaC=hypot(q.a,q.b); fftC=hypot(q.a_fft,q.b_fft);
        csv=sprintf('espectro_caso_%02d.csv',r);
        writetable(table(n,teoriaC,fftC),fullfile(pasta,csv));
        fprintf(fid,'\\begin{tikzpicture}\n\\begin{axis}[width=8cm,height=3cm,grid=major,');
        fprintf(fid,'xtick={1,2,3,4,5,6,7},xlabel={Ordem $n$},ylabel={$C_n$ (V)},font=\\scriptsize]\n');
        fprintf(fid,'\\addplot[blue,ycomb,mark=*,mark size=1pt] table[x=n,y=teoriaC,col sep=comma]{resultados/%s};\n',csv);
        fprintf(fid,'\\addplot[red,only marks,mark=o,mark size=1.6pt] table[x=n,y=fftC,col sep=comma]{resultados/%s};\n',csv);
        fprintf(fid,'\\end{axis}\n\\end{tikzpicture}\\\\\n');
    end
    fimImagem(fid); fclose(fid); compilarImagem(ficheiro,imagens);
    incluirImagem(report,arquivos{g},[captions{g}, ...
        ' Para cada sinal, apresentam-se a onda e seu espectro. Preto: amostras; azul: teoria até a ordem 7; vermelho: FFT da simulação.'], ...
        sprintf('fig:comparacao%d',g));
end
fclose(report);
assert(isequal(sort(cobertos),1:24),'Algum sinal ficou sem onda ou espectro.');
end
function titulo=nomeGrafico(r)
basicos={'Rampa crescente uni.','Rampa decrescente uni.','Triangular uni.', ...
    'Rampa crescente bi.','Rampa decrescente bi.','Triangular bi.'};
if r<=18
    titulo=basicos{mod(r-1,6)+1};
    if r>6 && r<=12, titulo=[titulo,', $\delta=P$, $\phi=P/10$']; end
    if r>12, titulo=[titulo,', $\delta=0{,}8P$, $\phi=P/10$']; end
else
    titulo=['Retangular ',rotulo(r)];
end
end

function inicioImagem(fid)
fprintf(fid,['\\documentclass[border=2pt]{standalone}\n\\usepackage[utf8]{inputenc}\n' ...
    '\\usepackage[T1]{fontenc}\n\\usepackage{lmodern,amsmath,pgfplots}\n' ...
    '\\pgfplotsset{compat=1.18}\n\\begin{document}\n\\begin{tabular}{@{}c@{}}\n']);
end

function fimImagem(fid)
fprintf(fid,'\\end{tabular}\n\\end{document}\n');
end

function compilarImagem(ficheiro,pasta)
cmd=sprintf('xelatex -no-shell-escape -interaction=batchmode -halt-on-error -output-directory="%s" "%s"',pasta,ficheiro);
[status,out]=system(cmd); assert(status==0,'Falha ao compilar figura: %s',out);
[~,base]=fileparts(ficheiro);
cmd=sprintf('pdftoppm -singlefile -png -r 180 "%s" "%s"',fullfile(pasta,[base,'.pdf']),fullfile(pasta,base));
[status,out]=system(cmd); assert(status==0,'Falha ao exportar PNG: %s',out);
end

function incluirImagem(fid,base,caption,label)
fprintf(fid,'\\begin{figure}[H]\n\\centering\n');
fprintf(fid,'\\includegraphics[width=\\linewidth,height=0.86\\textheight,keepaspectratio]{resultados/figuras/%s.pdf}\n',base);
fprintf(fid,'\\caption{%s}\\label{%s}\n\\end{figure}\n',caption,label);
end
