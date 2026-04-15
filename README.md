UERJ - Universidade do Estado do Rio de Janeiro, Rio de Janeiro (RJ), Brasil
PPG-CompMat - Programa de Pós-Graduação em Ciências Computacionais e Modelagem Matemática

Artigo: Uma Abordagem Computacional para Comparação Genética entre as variantes Delta e Omicron.
Autores: F. Farias, C. Fardilha, M. Gil, A. Sena.
Contato: farias.fabio@posgraduacao.uerj.br, camilafard@gmail.com, gilmanueljoao87@gmail.com, asena@ime.uerj.br

Esse projeto foi criado para atender atender a princípios valiosos da prática científica, são eles: transparência e reprodutibilidade.

A ordem para executar os procedimentos é a seguinte:

1) obter as sequências na base web do NCBI;
2) usar o arquivo 001_ncbi.R para ler as sequências e selecionar as amostras para as duas variantes e criar os arquivos fasta que serão usados no alinhamento das sequências;
3) adicionar no arquivo amostra.csv colunas com as posições iniciais e finais dos genes spike das sequências. Essas posições são encontradas ao digitar o código (accession) de cada sequência na barra de busca do NCBI;
4) usar o arquivo 002_spike.R para criar os arquivos fasta que serão usados no alinhamento do gene spike;
5) realizar os alinhamentos com ajuda da ferramenta MASA (linux), disponível em: https://github.com/edanssandes/MASACore, dois arquivos *.sh são usados: uma para executar o MASA para fazer os alinhamentos e outropara coletar os resultados;
6) usar o arquivo 003_analise.R para ler os resultados dos alinhamentos e faz toda a analise exploratória e inferência dos dados.
