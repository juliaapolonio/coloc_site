# Site Quarto — Atlas de Colocalização em AD

## Estrutura
```
site/
  _quarto.yml
  styles.css
  index.qmd
  manhattan.qmd
  bubble_explorer.qmd
  dev_adult.qmd
  cell_type.qmd
  locus_table.qmd
  data/                <- VOCÊ PRECISA CRIAR e colocar os arquivos aqui
    coloc_all_annotated.csv
    windows_com_nome_locus.tsv
    consensus_ADRD_withproxy.tsv
```

## Passo a passo

```bash
mkdir -p site/data
cp coloc_all_annotated.csv windows_com_nome_locus.tsv consensus_ADRD_withproxy.tsv site/data/

cd site
quarto render          # gera _site/ com o html estático
quarto preview          # pra conferir localmente antes de publicar
```

## Publicar (mais rápido: GitHub Pages ou Quarto Pub)

**Quarto Pub** (mais simples, um comando, hospedagem gratuita da Posit):
```bash
quarto publish quarto-pub
```

**GitHub Pages** (se já tem repositório):
```bash
quarto publish gh-pages
```

## Se o `consensus_ADRD_withproxy.tsv` for pesado demais para o navegador

O `manhattan.qmd` já filtra por p < 1e-4 antes de plotar, então o HTML final não carrega o arquivo bruto — só os pontos filtrados ficam embutidos. Mesmo assim, o `render` local vai processar o arquivo inteiro (pode levar alguns minutos). Se travar, reduza `p_display` para `1e-5` temporariamente só para testar o layout mais rápido.

## Prioridade se o tempo apertar

Para amanhã, a ordem de importância é: **index → bubble_explorer → manhattan**. As páginas de dev_adult/cell_type/locus_table são mais rápidas de renderizar (não leem o GWAS inteiro) — comece por elas se quiser garantir algo funcionando cedo.
