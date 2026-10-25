# Alterações e justificativas

## Cobertura integral dos sinais

A comparação foi ampliada para incluir todos os sinais no corpo do PDF. A seleção anterior de três exemplos foi removida.

- Os 24 cenários têm tabelas próprias com os sete pares de coeficientes teóricos e simulados: 168 pares ao todo.
- As energias total, da fundamental e das ordens 1–7 e os percentuais são apresentados para as seis ondas básicas, as seis deslocadas com delta=P, as seis com delta=0,8P e as seis retangulares.
- Todos os cenários têm a onda e o espectro de ordens 1–7, em oito figuras com seis gráficos cada: 48 gráficos ao todo.
- Foram conferidos 1098 valores nas 39 tabelas numéricas; a tabela simbólica completa o total de 40 tabelas.
- A planilha `resultados/cobertura_sinais.csv` identifica onde está a comparação de cada cenário.

Essas mudanças atendem à exigência de comparar todos os sinais e de apresentar suas energias, em vez de deixar parte da comparação apenas em arquivos externos. As simulações completas já salvas foram usadas para regenerar os resultados; nenhum caso foi preenchido por inferência a partir de outro sinal.

## Revisão atual: comparação direta entre teoria e simulação

Esta revisão substitui a apresentação anterior de erros por valores teóricos e simulados lado a lado. Os itens abaixo descrevem o estado atual; as demais seções registram a ampliação inicial com base no memorial.

| Alteração | Motivo |
|---|---|
| A FFT agora usa diretamente as amostras exportadas pelo Simulink | A comparação pedida é da teoria com a simulação |
| Removidas do relatório as diferenças entre gerações MATLAB e Simulink | As duas ferramentas são etapas da mesma análise, e não os dois resultados a comparar |
| Reexecutados os 24 casos e salvas as amostras no MAT | Garantir que coeficientes, energias e imagens simulados sejam medidos, não inferidos das fórmulas |
| Quatro tabelas com média e fundamental de todos os casos | Comparar a0 e C1 da teoria com os resultados simulados |
| 24 tabelas com todos os pares (an,bn) das ordens 1–7 | Mostrar valores e sinais de todos os cenários, incluindo zeros teóricos |
| Oito tabelas de energias e percentuais, com linhas T/S | Separar os valores analíticos dos estimados a partir das amostras e da FFT em todos os casos |
| Imagens em PDF e PNG com amostras, reconstruções e espectros | Facilitar a avaliação visual da concordância e do truncamento na ordem 7 |
| Reinserida a tabela dos ciclos alternados no desenvolvimento | A referência à Tabela 3 existia, mas sua inclusão havia sido removida |
| Corrigida a referência da equação do bipolar e a escrita de cos nas reconstruções | Eliminar referência indefinida e manter a notação matemática correta |
| Tabelas com vírgula decimal, unidades e identificação T/S; referências e largura revisadas | Melhorar legibilidade e conferir a representação no PDF |
| Ajustado o alinhamento vertical da seção de resultados | Evitar espaços verticais exagerados entre tabelas sem mudar a largura das colunas |

O código principal é `scripts/T1_Fourier_MATLAB.m`; a exportação das tabelas e imagens está em `scripts/T1_gerar_apresentacao.m`. As imagens são renderizadas a partir dos dados calculados no MATLAB, sem capturas de tela da interface.

O relatório foi completado a partir do arquivo `Memorial_Calculos_Detalhados.pdf` fornecido pelo usuário. As exigências acadêmicas descritas no início do memorial foram usadas como conteúdo da tarefa. Referências do memorial a scripts ou CSV de outro projeto não foram tratadas como arquivos disponíveis: a implementação e os resultados abaixo foram produzidos neste projeto.

## Atendimento ao enunciado

| Item do enunciado | O que foi entregue | Motivo |
|---|---|---|
| Oito ondas da Figura 1 | Mantidos os seis desenvolvimentos existentes; completada a retangular bipolar e a reconstrução das duas retangulares até n=7 | A seção do sinal 8 tinha apenas uma figura e faltavam as reconstruções retangulares |
| Amplitudes variáveis | Todas as expressões analíticas mantêm A como variável positiva | Os valores das figuras não fixam uma amplitude obrigatória |
| Casos phi=theta e phi diferente de theta | Coeficientes gerais e simétricos para as duas retangulares; exemplos com D=1/3; tabelas de todas as ordens 1–7 | Mostrar tanto a mudança de fase quanto a invariância das amplitudes para D fixo |
| Parâmetros primados nos pulsos | Modelo de dois ciclos que se repetem em 2P, com integrais dos dois pulsos e exemplo trabalhado | Um par de pulsos distintos não deve ser analisado automaticamente na base de período P |
| Seis ondas da Figura 2 | Equações periódicas, rotação dos coeficientes e seis séries truncadas | Essa parte do memorial ainda não estava no relatório |
| Casos delta=P e delta diferente de P | Período real, frequências e exemplos de 60 e 75 Hz | A mudança de período altera a base harmônica e a janela coerente da FFT |
| Energia da fundamental e das ordens 1–7 | Integração direta, Parseval, fórmulas e percentuais | O relatório anterior não calculava energia |
| Comparação com o sinal inteiro | Separação entre energia da média, das ordens 1–7 e da reconstrução com média | Evitar atribuir à soma dos harmônicos a energia da componente contínua |
| Comparação MATLAB/Simulink–FFT | 24 cenários executados, 168 pares de coeficientes e arquivos de erro | Produzir evidência reproduzível, em vez de apenas citar resultados do memorial |

## Mudanças por arquivo

| Arquivo | Alteração | Justificativa |
|---|---|---|
| `Main.tex` | Título definitivo; cabeçalho com departamento e Lista T1; preservados autores, logos e formato de duas colunas | Remover o título provisório e os dados de uma disciplina não relacionados ao material fornecido |
| `Main.tex` | Ativadas as seções de resultados, análise, conclusão e referências | Incorporar ao PDF o trabalho que antes ficava comentado |
| `Main.tex` | Corrigida a inclusão do arquivo `análise_resultados.tex`, respeitando o nome real | O caminho comentado anterior não correspondia ao arquivo existente |
| `Main.tex` | Usada a base bibliográfica `memorial.bib`; mantido `relatorio.bib` no disco | Referenciar a fonte realmente utilizada e não carregar uma base não citada, com aviso de sintaxe |
| `Capitulos/resumo.tex` | Preenchidos resumo e palavras-chave | O ambiente de resumo estava vazio |
| `Capitulos/introducao.tex` | Contextualização, objetivo e identificação do memorial | A seção continha somente o título |
| `Capitulos/definicao_problema.tex` | Requisitos, unidades e convenções de A, x, omega, P, L, d e delta | Explicitar as hipóteses necessárias antes das derivações |
| `Capitulos/solucao_problema.tex` | Inclusão dos fundamentos antes das ondas | Fixar a convenção de a0 como média e evitar fatores de 1/2 ambíguos |
| `Capitulos/solucao_problema.tex` | Triângulo unipolar escrito com `aligned`, separando a periodicidade | Manter o tamanho normal das fontes e a largura de uma coluna, sem `resizebox` |
| `Capitulos/solucao_problema.tex` | Explicitado o cancelamento das contribuições senoidais do triângulo unipolar | Detalhar a passagem das integrais para bn=0 |
| `Capitulos/solucao_problema.tex` | Numeração do sinal 3; redação corrigida de “sete termos” para ordem máxima 7; correção de “Desssa” | Incluir as ordens nulas na interpretação do truncamento e corrigir pequenos problemas editoriais |
| `Capitulos/solucao_problema.tex` | Largura retangular renomeada de delta para d, apenas nessa seção | Reservar delta ao período efetivo da Figura 2, como no memorial |
| `Capitulos/solucao_problema.tex` | Sinal 8 passou de subseção a subsubseção; inclusão do complemento após sua figura | Manter a mesma hierarquia dos demais sinais retangulares |
| `Capitulos/fundamentos_fourier.tex` | Série, integrais dos coeficientes, amplitude, fase e primitivas auxiliares | Dar suporte às simplificações usadas em todos os sinais |
| `Capitulos/complemento_calculos.tex` | Definição por trechos, média, integrais, limites e coeficientes do bipolar | Desenvolver os cálculos que faltavam, sem depender apenas de fórmulas finais |
| `Capitulos/complemento_calculos.tex` | Reconstruções retangulares gerais e exemplos simétrico/assimétrico | Explicitar a aplicação dos coeficientes até a ordem 7 |
| `Capitulos/complemento_calculos.tex` | Dois ciclos alternados e fundamental na base 2P | Considerar os intervalos primados e distinguir o período de repetição do intervalo de um pulso |
| `Capitulos/complemento_calculos.tex` | Seis ondas deslocadas, identidades de adição e tabela dos coeficientes | Demonstrar o efeito do deslocamento, em vez de apresentar apenas as expressões prontas |
| `Capitulos/complemento_calculos.tex` | Integrais de energia para rampas, triângulos e retangulares | Conferir independentemente o denominador dos percentuais de energia |
| `Capitulos/complemento_calculos.tex` | Parseval, somas finitas e percentuais com e sem CC | Atender à energia solicitada pelo enunciado sem misturar grandezas distintas |
| `Capitulos/complemento_calculos.tex` | Conversão de integral angular para V²s e para joules em uma resistência R | Não chamar uma integral de tensão ao quadrado de energia física sem definir a carga |
| `Capitulos/resultados_implementacao.tex` | Cenários, geração Simulink, derivação da normalização da FFT, erros, tabelas e gráficos | Descrever como os números foram realmente obtidos e como reproduzi-los |
| `Capitulos/análise_resultados.tex` | Comparações entre formas, média, deslocamento, período, alternância, Gibbs e discretização | Interpretar os resultados em relação ao enunciado; substitui um título genérico sem conteúdo |
| `Capitulos/conclusoes.tex` | Síntese do atendimento aos requisitos | A conclusão anterior estava vazia |
| `Capitulos/referencia.tex` e `memorial.bib` | Referência ao memorial por `biblatex`/Biber | Substituir referências de amplificadores operacionais que não eram usadas neste trabalho |
| `images/figura2_memorial-000.png` | Extração da figura original incorporada ao PDF | Documentar os seis sinais deslocados sem redesenhar ou alterar a figura fornecida |
| `scripts/T1_Fourier_MATLAB.m` | Programa completo de geração, cálculo analítico, FFT, simulação e exportação | Reproduzir os resultados do relatório e permitir novos valores dos parâmetros |
| `resultados/` | CSV, MAT, modelo Simulink, registro de execução, tabelas e dados dos gráficos | Entregar os dados por ordem e cenário e a evidência da execução |
| `README_VALIDACAO.md` | Instruções de reprodução e relação dos artefatos | Tornar o projeto utilizável sem depender desta conversa |
| `Main.pdf` | Versão compilada do relatório completo | Entregar o documento pronto para leitura |

Não foram adicionados ambientes `minipage` às derivações. As equações novas foram divididas com `aligned` e as tabelas e figuras mantêm a largura de uma coluna. Os `minipage` que já compunham o cabeçalho com os logos foram preservados.

## Conferências matemáticas e cuidados

1. A retangular bipolar é `v8=2*v7-A`; portanto, sua média é `A*(2D-1)`, e só é nula para D=1/2.
2. Os coeficientes harmônicos bipolares são o dobro dos unipolares. Não se multiplica a média unipolar por dois sem também subtrair A.
3. As fases seguem `v_n=C_n*cos(k_n*x-alpha_n)`, com `alpha_n=atan2(b_n,a_n)`. Harmônicos de amplitude nula não têm fase definida.
4. Os desenhos com intervalos primados foram tratados sob a hipótese explícita de repetição do par de ciclos. Se eles não se repetirem, dois pulsos desenhados não bastam para definir uma função periódica.
5. Foi explicitada a possível redução do período 2P para P quando os dois ciclos se tornam iguais.
6. Os percentuais “somente ordens 1–7” e “média mais ordens 1–7” são ambos apresentados. Nas três ondas unipolares em rampa/triângulo, a componente contínua corresponde a 75% da energia total.
7. A correção de fase `exp(-i*pi*n/N)` é aplicada porque os exemplos usam pontos médios. Ela deve ser omitida quando a aquisição começa exatamente em x=0.
8. Os valores de A, P, omega e dos intervalos dos exemplos foram identificados como escolhas de validação, não como dados impostos pelo enunciado.

## Validação executada

- MATLAB R2026a Update 4, com licença e instalação do Simulink verificadas.
- Execução completa dos 24 cenários no MATLAB e no Simulink, com N=65536 amostras por período real.
- Comparação de a0 e dos 168 pares `(an,bn)`, ordens 1 a 7.
- Média e fundamental teóricas e simuladas apresentadas para todos os 24 cenários.
- Energias totais medidas nas amostras e energias harmônicas extraídas dos coeficientes FFT, comparadas com as integrais analíticas.
- Removida da apresentação a comparação entre uma geração MATLAB e uma geração Simulink.
- Verificado que a energia da reconstrução até a ordem 7 não excede a energia total.
- As tabelas de energia coincidem com os exemplos do memorial, mantendo separada a parcela contínua.
- Os 12 exemplos da tabela final do memorial foram conferidos automaticamente, respeitando o arredondamento a quatro casas decimais.
- PDF final com 21 páginas e 40 tabelas, compilado com XeLaTeX e Biber, sem erros, referências indefinidas ou avisos de conteúdo excedendo a largura.
- Conferidos 1098 valores publicados nas 39 tabelas numéricas contra os CSV, respeitando o arredondamento. A tabela simbólica de deslocamentos e as páginas com as comparações foram revisadas visualmente.

Os valores publicados são resultados desta execução. As diferenças nas últimas casas entre teoria e simulação são compatíveis com a discretização das fronteiras dos sinais; não representam uma correção das fórmulas analíticas. Os registros e arquivos completos estão em `resultados/`.
