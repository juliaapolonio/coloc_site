#!/usr/bin/env Rscript
# =====================================================================
# TRIANGULAÇÃO: colocalização dev-only x cromatina de desenvolvimento
# =====================================================================
# Lógica: um par locus|gene dev-only ganha suporte adicional se, no MESMO
# contexto de desenvolvimento:
#   (a) a variante/janela do sinal de sQTL cai em cromatina aberta (ATAC-seq
#       fetal)
#   (b) esse elemento aberto faz um LOOP 3D (HiChIP/Hi-C fetal) até o
#       promotor do gene candidato
#
# Isso não prova causalidade, mas eleva "colocalização estatística" para
# "colocalização + elemento regulatório + conectividade física" -- é o
# tipo de evidência triangulada que costuma ser pedida por revisor.
#
# ENTRADAS QUE VOCÊ PRECISA APONTAR:
#   - peaks_files: 1+ BED de picos de ATAC-seq (idealmente já no formato
#     narrowPeak ou BED simples chr/start/end)
#   - loops_files: 1+ arquivo de loops HiChIP/Hi-C no formato .bedpe
#     (chr1 start1 end1 chr2 start2 end2 [...score])
#   - o splice_id/coordenada do evento dev-only (já preservado no
#     aggregate_coloc.R) -- é o que ancora a variante no espaço genômico
#
# SAÍDA: uma tabela por par locus|gene com colunas booleanas de suporte e
# um score agregado (0-3).
# =====================================================================

pkgs <- c("data.table", "stringr", "GenomicRanges", "org.Hs.eg.db", "AnnotationDbi")
to_install_cran <- setdiff(c("data.table","stringr"), rownames(installed.packages()))
if (length(to_install_cran) > 0) install.packages(to_install_cran, repos = "https://cloud.r-project.org")
to_install_bioc <- setdiff(c("GenomicRanges","org.Hs.eg.db","AnnotationDbi","GenomicFeatures",
                             "TxDb.Hsapiens.UCSC.hg38.knownGene"), rownames(installed.packages()))
if (length(to_install_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(to_install_bioc, update = FALSE, ask = FALSE)
}
library(data.table); library(stringr); library(GenomicRanges)
library(GenomicFeatures); library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db); library(AnnotationDbi)

## ---------------------------------------------------------------
## 1) Parâmetros                                          ### AJUSTAR ###
## ---------------------------------------------------------------
coloc_file <- "coloc_all_annotated.csv"

# um ou mais arquivos de picos de ATAC-seq de desenvolvimento (BED: chr,start,end,...)
# ex.: cortical plate, germinal zone, fetal microglia -- rotule por tecido/estágio
peaks_files <- c(
  cortical_plate = "atac_fetal_cortical_plate.bed",
  germinal_zone  = "atac_fetal_germinal_zone.bed"
)

# um ou mais arquivos de loop (.bedpe: chr1 start1 end1 chr2 start2 end2 [score])
loops_files <- c(
  hichip_fetal = "hichip_fetal_loops.bedpe"
)

promoter_flank <- 2000   # pb ao redor da TSS considerados "promotor"
loop_slop      <- 5000   # tolerância (pb) ao casar âncora de loop com peak/promotor
pph4_threshold <- 0.8

## ---------------------------------------------------------------
## 2) Carregar os pares dev-only e suas coordenadas
## ---------------------------------------------------------------
d <- fread(coloc_file)
d[is.na(gene_symbol) | gene_symbol == "", gene_symbol := gene_clean]
d[, gene_symbol := str_remove(as.character(gene_symbol), ";.*$")]
d[, pair := paste(locus, gene_symbol, sep = " | ")]
d[, stage := fifelse(str_detect(dataset, regex("^devbrain_", ignore_case = TRUE)), "Dev", "Adulto")]

sig <- d[PP.H4 > pph4_threshold]
w <- dcast(sig[, .(pp4 = max(PP.H4, na.rm = TRUE)), by = .(pair, stage)],
           pair ~ stage, value.var = "pp4", fill = 0)
dev_only_pairs <- w[Dev > pph4_threshold & Adulto <= pph4_threshold, pair]
message(length(dev_only_pairs), " pares dev-only para triangular")

# stopifnot(exists("splice_id")) -- ver aggregate_coloc.R; sem essa coluna,
# a coordenada do EVENTO não existe e a triangulação cai pra nível de janela
tem_splice_id <- "splice_id" %in% names(d)
if (!tem_splice_id) {
  warning("Coluna 'splice_id' ausente -- a triangulação vai usar a JANELA do ",
          "locus inteira em vez da coordenada exata do evento de splicing. ",
          "Reprocesse o aggregate_coloc.R com a correção do splice_id para ",
          "resolução fina (ver conversa anterior).")
}

evt <- sig[pair %in% dev_only_pairs & stage == "Dev"]
evt <- evt[order(-PP.H4), .SD[1], by = pair]  # 1 linha (melhor evento) por par

parse_coord <- function(x) {
  x <- as.character(x)
  lc <- str_match(x, "^([^:]+):([0-9]+):([0-9]+):(clu_[0-9]+)")
  bb <- str_match(x, "^(chr[0-9XYM]+)_([0-9]+)_([0-9]+)$")
  data.table(
    chr   = str_remove(fifelse(!is.na(lc[,2]), lc[,2], bb[,2]), "^chr"),
    start = as.integer(fifelse(!is.na(lc[,3]), lc[,3], bb[,3])),
    end   = as.integer(fifelse(!is.na(lc[,4]), lc[,4], bb[,4]))
  )
}

if (tem_splice_id) {
  coord <- parse_coord(evt$splice_id)
  falhou <- is.na(coord$chr)
  if (any(falhou)) {
    message(sum(falhou), " evento(s) sem coordenada parseável -- usando janela do locus")
    coord[falhou, `:=`(chr = evt$region_chr[falhou], start = evt$region_start[falhou],
                       end = evt$region_end[falhou])]
  }
} else {
  coord <- data.table(chr = evt$region_chr, start = evt$region_start, end = evt$region_end)
}
evt <- cbind(evt, coord)
evt_gr <- GRanges(paste0("chr", evt$chr), IRanges(evt$start, evt$end))

## ---------------------------------------------------------------
## 3) Promotor do gene candidato (TSS ± flanco)
## ---------------------------------------------------------------
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)
genes_gr$symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = genes_gr$gene_id,
                                         keytype = "ENTREZID", column = "SYMBOL",
                                         multiVals = "first")

tss <- ifelse(strand(genes_gr) == "+", start(genes_gr), end(genes_gr))
prom_gr <- GRanges(seqnames(genes_gr), IRanges(tss - promoter_flank, tss + promoter_flank),
                   strand = strand(genes_gr), symbol = genes_gr$symbol)

prom_evt <- prom_gr[match(evt$gene_symbol, prom_gr$symbol)]
sem_promotor <- is.na(start(prom_evt))
if (any(sem_promotor)) {
  message(sum(sem_promotor), " gene(s) sem promotor localizado no TxDb (symbol não bateu): ",
          paste(evt$gene_symbol[sem_promotor], collapse = ", "))
}

## ---------------------------------------------------------------
## 4) Suporte (a): evento cai em cromatina aberta?
## ---------------------------------------------------------------
ler_peaks <- function(f) {
  if (!file.exists(f)) { warning("Peaks não encontrado: ", f); return(NULL) }
  p <- fread(f, header = FALSE)
  setnames(p, 1:3, c("chr","start","end"))
  GRanges(ifelse(str_starts(p$chr, "chr"), p$chr, paste0("chr", p$chr)),
          IRanges(p$start, p$end))
}
peaks_list <- Filter(Negate(is.null), lapply(peaks_files, ler_peaks))

evt[, atac_suporte := FALSE]
evt[, atac_tecido := NA_character_]
for (nm in names(peaks_list)) {
  hit <- overlapsAny(evt_gr, peaks_list[[nm]])
  evt[hit == TRUE & atac_suporte == FALSE, atac_tecido := nm]
  evt[, atac_suporte := atac_suporte | hit]
}
message("Eventos com suporte de ATAC-seq fetal: ", sum(evt$atac_suporte), "/", nrow(evt))

## ---------------------------------------------------------------
## 5) Suporte (b): loop conectando o evento ao promotor do gene
## ---------------------------------------------------------------
ler_loops <- function(f) {
  if (!file.exists(f)) { warning("Loops não encontrado: ", f); return(NULL) }
  L <- fread(f, header = FALSE)
  setnames(L, 1:6, c("chr1","start1","end1","chr2","start2","end2"))
  anc1 <- GRanges(ifelse(str_starts(L$chr1,"chr"), L$chr1, paste0("chr", L$chr1)),
                  IRanges(L$start1, L$end1))
  anc2 <- GRanges(ifelse(str_starts(L$chr2,"chr"), L$chr2, paste0("chr", L$chr2)),
                  IRanges(L$start2, L$end2))
  list(anc1 = anc1, anc2 = anc2, id = seq_len(nrow(L)))
}
loops_list <- Filter(Negate(is.null), lapply(loops_files, ler_loops))

evt[, loop_suporte := FALSE]
evt[, loop_dataset := NA_character_]
for (nm in names(loops_list)) {
  Lp <- loops_list[[nm]]
  # evento numa âncora + promotor na outra (nas duas orientações), com folga
  ev_slop   <- evt_gr + loop_slop
  prom_slop <- prom_evt + loop_slop

  hit1 <- overlapsAny(Lp$anc1, ev_slop) & overlapsAny(Lp$anc2, prom_slop, maxgap = 0)
  hit2 <- overlapsAny(Lp$anc2, ev_slop) & overlapsAny(Lp$anc1, prom_slop, maxgap = 0)
  loop_ok <- hit1 | hit2

  # para cada evento, existe algum loop conectando-o ao SEU promotor específico?
  achou <- vapply(seq_len(nrow(evt)), function(i) {
    e_i <- ev_slop[i]; p_i <- prom_slop[i]
    if (is.na(start(p_i))) return(FALSE)
    a <- overlapsAny(Lp$anc1, e_i) & overlapsAny(Lp$anc2, p_i)
    b <- overlapsAny(Lp$anc2, e_i) & overlapsAny(Lp$anc1, p_i)
    any(a | b)
  }, logical(1))

  evt[achou & loop_suporte == FALSE, loop_dataset := nm]
  evt[, loop_suporte := loop_suporte | achou]
}
message("Eventos com suporte de loop 3D fetal: ", sum(evt$loop_suporte), "/", nrow(evt))

## ---------------------------------------------------------------
## 6) Score agregado e tabela final
## ---------------------------------------------------------------
evt[, score_triangulacao := 1 + atac_suporte + loop_suporte]  # 1 = só coloc; até 3
setorder(evt, -score_triangulacao, -PP.H4)

saida <- evt[, .(locus, gene_symbol, pair, PP.H4 = round(PP.H4, 3),
                 chr, start, end, atac_suporte, atac_tecido,
                 loop_suporte, loop_dataset, score_triangulacao)]
fwrite(saida, "triangulacao_dev_only.csv")
message("Tabela salva em: triangulacao_dev_only.csv")
print(saida[, .(pair, PP.H4, atac_suporte, loop_suporte, score_triangulacao)])

message("
Leitura do score:
  1 = só colocalização estatística (nenhum suporte epigenômico ainda)
  2 = colocalização + cromatina aberta (ATAC) OU + loop
  3 = colocalização + cromatina aberta + loop conectando ao promotor do gene
     (evidência mais forte -- variante regulatória num elemento ativo que
     fisicamente contata o promotor do gene candidato, no mesmo contexto
     de desenvolvimento em que o sinal genético foi detectado)
")
