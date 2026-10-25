# Reprodução do relatório

O documento principal é `Main.tex`; o PDF pronto é `Main.pdf`. A relação das alterações e suas justificativas está em `ALTERACOES_RELATORIO.md`.

## Produzir as imagens na interface do MATLAB

O arquivo independente `scripts/T1_simulacoes_imagens.m` contém as 24 simulações, a FFT, a teoria e a exportação gráfica nativa do MATLAB. Abra-o e clique em **Run**, ou execute:

```matlab
addpath('scripts');
T1_simulacoes_imagens
```

As imagens PNG, PDF e FIG serão gravadas em `resultados/imagens_matlab/`. Para utilizar as amostras já obtidas no Simulink, execute `T1_simulacoes_imagens('Modo','salvos')`. As instruções completas e as opções estão em `scripts/COMO_GERAR_IMAGENS_MATLAB.md`.

Os cálculos desse arquivo foram verificados nos 24 casos, tanto no modo MATLAB quanto no modo de dados salvos, conferindo 168 pares de coeficientes e as energias com os resultados existentes. A exportação por `plot`/`exportgraphics` deve ser executada na sessão gráfica do MATLAB; essa etapa não foi validada no ambiente remoto usado para conferir os cálculos.

## Gerar os resultados

Na pasta do projeto, execute no MATLAB:

```matlab
addpath('scripts');
T1_Fourier_MATLAB
```

O programa foi executado no MATLAB R2026a Update 4 e no Simulink. Ele usa A=1 V, P=2π, omega=120π rad/s e N=65536. Esses parâmetros são exemplos e podem ser alterados no início do script. As expressões analíticas do relatório continuam simbólicas.

São 24 cenários: seis ondas da Figura 1; seis da Figura 2 com delta=P; seis com delta=0,8P; e seis retangulares, cobrindo as versões unipolar e bipolar e os casos simétrico, assimétrico e alternado. O bloco MATLAB Function gera as ondas por suas equações; ele não usa os coeficientes teóricos para produzir o sinal a ser validado.

O modelo salvo em `resultados/T1_validacao_sinais.slx` representa o último cenário executado. Para percorrer todos os cenários e atualizar os parâmetros do bloco, execute o script.

Se o Simulink estiver indisponível em outra instalação, o programa executa somente o MATLAB e registra explicitamente essa condição. Isso não altera o registro entregue, que contém a execução real com Simulink disponível.

Para recriar somente as tabelas e os gráficos a partir dos dados já salvos:

```matlab
T1_Fourier_MATLAB('artefatos')
```

## Compilar o PDF

Na pasta do projeto, execute:

```powershell
xelatex -no-shell-escape -interaction=nonstopmode -halt-on-error Main.tex
biber Main
xelatex -no-shell-escape -interaction=nonstopmode -halt-on-error Main.tex
xelatex -no-shell-escape -interaction=nonstopmode -halt-on-error Main.tex
```

Também é possível usar uma receita de `latexmk` com XeLaTeX e Biber. Não é necessário executar novamente o MATLAB para compilar o relatório: os CSV e os arquivos TeX gerados acompanham o projeto.

Os gráficos usam as amostras reais das simulações e os coeficientes de FFT calculados no MATLAB. São exportados em PDF vetorial e PNG para `resultados/figuras/`, por meio de PGFPlots, XeLaTeX e `pdftoppm`; não dependem de janelas gráficas do MATLAB. O programa `scripts/T1_gerar_apresentacao.m` automatiza essa exportação e a criação das tabelas. Para regenerar as imagens, XeLaTeX e `pdftoppm` precisam estar disponíveis no PATH. Os caches de compilação do Simulink são encaminhados para a pasta temporária do sistema.

## Arquivos de resultados

| Arquivo | Conteúdo |
|---|---|
| `coeficientes_fft.csv` | 168 linhas: coeficientes teóricos, coeficientes da FFT das amostras simuladas e amplitudes |
| `energias.csv` | Energias e percentuais teóricos e simulados lado a lado, além da frequência, média e conversão temporal |
| `validacao_fft.csv` | Comparação direta da média, fundamental e valor quadrático médio: teoria versus simulação |
| `validacao_completa.mat` | Parâmetros, amostras simuladas, coeficientes e energias necessários à reprodução dos artefatos |
| `T1_validacao_sinais.slx` | Modelo Simulink com geração por equações, amostragem e exportação |
| `execucao_matlab.txt` | Registro das execuções; na revisão atual, os valores teóricos e simulados são impressos lado a lado |
| `ambiente.txt` | Versão do MATLAB, N, omega e informação sobre a execução do Simulink |
| `compilacao.txt` | Conferência final do PDF, dos resultados e hash SHA256 do documento entregue |
| `retangular_coef_1.tex` a `retangular_coef_3.tex` | Sete pares de coeficientes de cada exemplo retangular unipolar |
| `tabela_energias.tex`, `tabela_validacao.tex`, `resumo_validacao.tex` | Energias, percentuais, médias e fundamentais comparados diretamente com a teoria |
| `tabela_coeficientes_comparacao.tex` | 24 tabelas: comparação dos sete pares de coeficientes de todos os cenários |
| `onda_*.csv` | Amostras da simulação e reconstruções de ordem 7 com coeficientes teóricos e medidos pela FFT |
| `espectro_*.csv` | Amplitudes analíticas e de FFT para as figuras espectrais |
| `graficos_fourier.tex` | Inclusão das oito figuras comparativas, cada uma com três ondas e seus três espectros |
| `figuras/comparacao_*.pdf` e `figuras/comparacao_*.png` | Oito conjuntos atuais com onda e espectro dos 24 cenários, agrupados em três sinais por imagem |
| `cobertura_sinais.csv` | Conferência dos 24 cenários: tabela de coeficientes, energia, percentuais, onda e espectro de cada um |

Todos os resultados estão presentes no PDF, e não somente nos CSV: são 24 tabelas com os sete pares de coeficientes, quatro tabelas de energia absoluta, quatro de percentuais e 48 gráficos (24 ondas e 24 espectros). As quatro tabelas de média/fundamental e as quatro tabelas teóricas do desenvolvimento completam o total de 40 tabelas.

Os arquivos `onda_*.csv` preservam uma seleção das amostras reais para desenhar as curvas; a FFT e as energias usam todas as 65536 amostras. As tabelas identificam T (teoria) e S (simulação). A energia total simulada é `L*mean(y.^2)`; a energia da fundamental e das ordens 1–7 usa os coeficientes medidos pela FFT, e os percentuais simulados usam a energia total simulada no denominador. Não se compara o MATLAB com o Simulink: o Simulink gera as amostras, e o MATLAB as analisa e compara com a teoria.

`E_total_angular` e `E_total_sim` têm unidade V²rad. `E_total_temporal` e sua versão simulada têm unidade V²s; para joules sobre uma resistência R, divide-se a grandeza temporal por R.
