# Gerar as imagens diretamente no MATLAB

Abra `T1_simulacoes_imagens.m` no MATLAB e clique em **Run**. Também pode executar, na pasta do projeto:

```matlab
addpath('scripts');
T1_simulacoes_imagens
```

O arquivo é independente dos outros scripts no modo padrão. Ele simula numericamente todos os 24 cenários, calcula a FFT e exporta uma figura por cenário, usando apenas as ferramentas gráficas do MATLAB. Não chama LaTeX nem Simulink nesse modo. Os parâmetros padrão são A=1 V, P=2π, omega=120π rad/s e N=65536.

Cada figura contém a onda simulada, as reconstruções teórica e por FFT até a ordem 7, o espectro das sete ordens e as energias teóricas e simuladas. As energias mostradas no gráfico estão em V²s; as tabelas CSV usam V²rad.

Os arquivos ficam em `resultados/imagens_matlab/`:

- **PNG**, com resolução padrão de 300 dpi;
- **PDF vetorial**, recomendado para inserir no relatório;
- **FIG**, para reabrir e editar no MATLAB;
- CSV de coeficientes e energias e MAT com os dados;
- `figuras_para_relatorio.tex`, com exemplos de inclusão das imagens.

Para usar as amostras das simulações do Simulink já realizadas:

```matlab
T1_simulacoes_imagens('Modo','salvos')
```

Para selecionar casos e manter suas figuras abertas:

```matlab
T1_simulacoes_imagens('Cenarios',[1 3 6 19 24],'MostrarFiguras',true)
```

Para alterar os parâmetros ou a pasta de saída:

```matlab
T1_simulacoes_imagens('A',2,'N',2^17,'PastaSaida',fullfile(pwd,'minhas_imagens'))
```

Os parâmetros A, P, Omega e N só se aplicam à geração nova no modo `matlab`. O modo `salvos` utiliza os parâmetros e as amostras contidos no arquivo MAT, sem alterá-los.

Para executar somente os cálculos:

```matlab
[coeficientes, energias] = T1_simulacoes_imagens('GerarImagens',false);
```

São contemplados os seis sinais básicos, os seis deslocados com delta=P, os seis com delta=0,8P e as retangulares unipolares e bipolares nos casos simétrico, assimétrico e alternado. Os números de caso, de 1 a 24, seguem a mesma ordem de `resultados/cobertura_sinais.csv`.

Para inserir uma imagem no LaTeX, copie o PDF desejado para `images/matlab/` e use:

```latex
\begin{figure}[htbp]
    \centering
    \includegraphics[width=\linewidth]{images/matlab/01_Fig1_sinal1.pdf}
    \caption{Rampa crescente unipolar: comparação entre teoria e simulação.}
    \label{fig:matlab_rampa_crescente}
\end{figure}
```

As figuras são geradas com `figure`, `tiledlayout`, `plot`, `stem`, `bar` e `exportgraphics`. Recomenda-se executar a exportação na sessão gráfica do MATLAB. O modo de cálculos pode ser executado sem abrir figuras.

Os CSV exportados incluem o erro absoluto de cada par de coeficientes,
o erro relativo da fundamental e da energia total e o THD teórico e simulado.
`THD7` considera as ordens 2 a 7; o THD total exclui a componente contínua
e usa a energia total. Na simulação, essa última medida é uma estimativa
limitada à aquisição finita. Todos os erros comparam teoria e simulação.
