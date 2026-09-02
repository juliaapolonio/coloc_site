#!/usr/bin/env Rscript
# =====================================================================
# FIGURAS PARA O PAPER
#   Fig 1 — Panorama: cobertura + taxa de detecção normalizada
#   Fig 2 — Mecanismo: eQTL vs sQTL no BigBrain (contraste pareado)
#   Fig 3 — Especificidade celular (heatmap estático, publicação)
#   Fig 4 — Desenvolvimento vs adulto (scatter com destaques)
# Fig 5 (nomenclatura / PP.H3) já está em fig_slide_pph3.R -- reusar de lá.
# =====================================================================

pkgs <- c("data.table", "stringr", "ggplot2", "ggrepel", "patchwork")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")
library(data.table); library(stringr); library(ggplot2); library(ggrepel); library(patchwork)

## ---------------------------------------------------------------
in_file <- "coloc_all_annotated.csv"     ### AJUSTAR ###
THR <- 0.8
out_dir <- "figuras_paper"; dir.create(out_dir, showWarnings = FALSE)

d <- fread(in_file)
d[is.na(gene_symbol) | gene_symbol == "", gene_symbol := gene_clean]
d[, gene_symbol := str_remove(as.character(gene_symbol), ";.*$")]
d[, pair := paste(locus, gene_symbol, sep = " | ")]
d[, family := fcase(
  str_detect(dataset, regex("^gtex_", ignore_case = TRUE)),        "GTEx",
  str_detect(dataset, regex("^bigbrain_", ignore_case = TRUE)),    "BigBrain",
  str_detect(dataset, regex("^devbrain_", ignore_case = TRUE)),    "DevBrain",
  str_detect(dataset, regex("^singlebrain_", ignore_case = TRUE)), "SingleBrain",
  default = "Outro")]
d[, mod   := fifelse(str_detect(dataset, regex("_sqtl$", ignore_case = TRUE)), "sQTL", "eQTL")]
d[, stage := fifelse(str_detect(dataset, regex("^devbrain_", ignore_case = TRUE)), "Dev", "Adulto")]
sig <- d[PP.H4 > THR]

tema <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

## =====================================================================
## FIG 1 — Panorama
## =====================================================================
rate <- merge(d[, .(n_testes = .N), by = .(dataset, family)],
              sig[, .(n_sig = .N), by = dataset], by = "dataset", all.x = TRUE)
rate[is.na(n_sig), n_sig := 0]
rate[, taxa := n_sig / n_testes * 1000]
setorder(rate, family, -taxa)
rate[, dataset := factor(dataset, levels = rev(dataset))]

p1a <- ggplot(rate, aes(x = taxa, y = dataset, fill = family)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = c(GTEx = "#4E79A7", BigBrain = "#F28E2B",
                               DevBrain = "#59A14F", SingleBrain = "#E15759"), name = NULL) +
  labs(x = "Colocalizações sig. por 1.000 testes", y = NULL,
       title = "A") +
  tema + theme(axis.text.y = element_text(size = 6), legend.position = "top")

cov <- data.table(
  cat = c("Loci testados", "Com colocalização", "Genes priorizados"),
  n = c(uniqueN(d$locus), uniqueN(sig$locus), uniqueN(sig$gene_symbol))
)
cov[, cat := factor(cat, levels = rev(cat))]
p1b <- ggplot(cov, aes(x = n, y = cat)) +
  geom_col(fill = "#1F6FB4", width = 0.6) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = NULL, title = "B") + tema

fig1 <- p1a + p1b + plot_layout(widths = c(2, 1))
ggsave(file.path(out_dir, "Fig1_panorama.pdf"), fig1, width = 11, height = 6)
ggsave(file.path(out_dir, "Fig1_panorama.png"), fig1, width = 11, height = 6, dpi = 300, bg = "white")

## =====================================================================
## FIG 2 — Mecanismo: eQTL vs sQTL (BigBrain)
## =====================================================================
bb <- d[family == "BigBrain"]
m <- dcast(bb[, .(pp4 = max(PP.H4, na.rm = TRUE)), by = .(pair, mod)],
           pair ~ mod, value.var = "pp4", fill = 0)
m <- m[eQTL > THR | sQTL > THR]
m[, classe := fcase(eQTL > THR & sQTL > THR, "Ambos",
                    sQTL > THR, "Só sQTL", default = "Só eQTL")]

lab2 <- m[classe == "Ambos"]

fig2 <- ggplot(m, aes(x = eQTL, y = sQTL, colour = classe)) +
  geom_hline(yintercept = THR, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = THR, linetype = "dashed", colour = "grey60") +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(data = lab2, aes(label = pair), size = 2.6, max.overlaps = 15,
                  show.legend = FALSE, seed = 1) +
  scale_colour_manual(values = c("Ambos" = "#1F6FB4", "Só sQTL" = "#C0392B",
                                 "Só eQTL" = "#95A5A6"), name = NULL) +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Melhor PP.H4 — eQTL", y = "Melhor PP.H4 — sQTL",
       title = paste0("BigBrain (mesma coorte): n = ", nrow(m), " pares")) +
  tema + theme(legend.position = "top")

ggsave(file.path(out_dir, "Fig2_mecanismo.pdf"), fig2, width = 7.5, height = 7)
ggsave(file.path(out_dir, "Fig2_mecanismo.png"), fig2, width = 7.5, height = 7, dpi = 300, bg = "white")

## =====================================================================
## FIG 3 — Especificidade celular (heatmap estático)
## =====================================================================
sc <- d[family == "SingleBrain"]
sc[, ct := str_remove(dataset, regex("^singlebrain_", ignore_case = TRUE))]
ctm <- sc[, .(pp4 = max(PP.H4, na.rm = TRUE)), by = .(pair, ct)]
ctm <- ctm[pair %in% unique(ctm[pp4 > THR, pair])]

ord3 <- ctm[, .(n_sig = sum(pp4 > THR), top = max(pp4)), by = pair][order(n_sig, -top)]
ctm[, pair := factor(pair, levels = ord3$pair)]

# limita a publicação a pares com maior sinal se a lista for muito longa
if (uniqueN(ctm$pair) > 60) {
  keep3 <- ord3[1:60, pair]
  ctm <- ctm[pair %in% keep3]
}

fig3 <- ggplot(ctm, aes(x = ct, y = pair, fill = pp4)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = ctm[pp4 > THR], shape = 21, size = 1, fill = "white", colour = "white") +
  scale_fill_gradientn(colours = c("grey95", "#c6dbef", "#4292c6", "#08306b"),
                       limits = c(0, 1), name = "PP.H4") +
  labs(x = "Tipo celular", y = NULL) +
  tema + theme(axis.text.y = element_text(size = 6), panel.grid = element_blank())

ggsave(file.path(out_dir, "Fig3_tipo_celular.pdf"), fig3, width = 7, height = 11)
ggsave(file.path(out_dir, "Fig3_tipo_celular.png"), fig3, width = 7, height = 11, dpi = 300, bg = "white")

## =====================================================================
## FIG 4 — Desenvolvimento vs adulto
## =====================================================================
w4 <- dcast(d[, .(pp4 = max(PP.H4, na.rm = TRUE)), by = .(pair, stage)],
            pair ~ stage, value.var = "pp4", fill = 0)
w4 <- w4[Adulto > THR | Dev > THR]
w4[, classe := fcase(Adulto > THR & Dev > THR, "Ambos",
                     Dev > THR, "Só desenvolvimento", default = "Só adulto")]

destaque <- c("CELF1/SPI1 | SLC39A13", "KAT8/BCKDK | BCKDK", "MYO15A/TOM1L2 | SHMT1",
             "ACE | ERN1", "APP | LINC00158")
lab4 <- w4[pair %in% destaque]

fig4 <- ggplot(w4, aes(x = Adulto, y = Dev)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey70") +
  geom_hline(yintercept = THR, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = THR, linetype = "dashed", colour = "grey60") +
  geom_point(aes(colour = classe), alpha = 0.7, size = 2) +
  geom_point(data = lab4, shape = 21, size = 3.5, colour = "black", fill = NA, stroke = 1) +
  geom_text_repel(data = lab4, aes(label = str_remove(pair, ".*\\| ")),
                  size = 3, fontface = "bold", seed = 1, max.overlaps = 20) +
  scale_colour_manual(values = c("Ambos" = "#7F8C8D", "Só adulto" = "#B7C9DB",
                                 "Só desenvolvimento" = "#C0392B"), name = NULL) +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Melhor PP.H4 — Adulto", y = "Melhor PP.H4 — Desenvolvimento",
       title = paste0("n = ", nrow(w4), " pares")) +
  tema + theme(legend.position = "top")

ggsave(file.path(out_dir, "Fig4_dev_adulto.pdf"), fig4, width = 7.5, height = 7)
ggsave(file.path(out_dir, "Fig4_dev_adulto.png"), fig4, width = 7.5, height = 7, dpi = 300, bg = "white")

message("Figuras salvas em: ", normalizePath(out_dir))
message("Fig 5 (nomenclatura/PP.H3): usar fig_slide_pph3.R (já existente).")
