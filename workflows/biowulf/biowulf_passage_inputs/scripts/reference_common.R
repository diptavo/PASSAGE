suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

canonical_markers <- list(
  breast = list(
    malignant_epithelial = c("EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "MUC1", "ERBB2", "ESR1", "PGR", "GATA3", "FOXA1", "KRT5", "KRT14", "KRT17", "MKI67", "TOP2A"),
    normal_epithelial = c("EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "KRT5", "KRT14", "TP63", "KRT17", "ACTA2", "MYH11"),
    t_nk = c("CD3D", "CD3E", "TRAC", "CD247", "NKG7", "GNLY", "KLRD1", "CCL5", "GZMB"),
    b_plasma = c("MS4A1", "CD79A", "CD74", "CD37", "MZB1", "JCHAIN", "IGHG1", "SDC1"),
    myeloid = c("LST1", "TYROBP", "AIF1", "FCER1G", "C1QA", "C1QB", "CTSS", "SPP1", "CD68"),
    fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "FAP", "CXCL12", "TAGLN"),
    endothelial = c("PECAM1", "VWF", "EMCN", "ENG", "KDR", "FLT1", "RAMP2", "PLVAP", "ACKR1"),
    perivascular = c("RGS5", "PDGFRB", "MCAM", "CSPG4", "NOTCH3", "ACTA2", "MYH11", "TAGLN", "DES")
  ),
  kidney = list(
    malignant_epithelial = c("CA9", "NDUFA4L2", "VEGFA", "EGLN3", "ANGPTL4", "SLC2A1", "ENO1", "PAX8", "EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "AMACR"),
    tubular_epithelial = c("LRP2", "CUBN", "SLC34A1", "SLC5A2", "ALDOB", "AQP1", "UMOD", "SLC12A1", "CLDN16", "KCNJ1", "SLC12A3", "CALB1", "TRPM6", "PVALB", "NPHS1", "NPHS2", "PODXL", "WT1", "SYNPO"),
    collecting_duct = c("AQP2", "AQP3", "KRT8", "KRT18", "ATP6V1B1", "SLC4A1", "FOXI1", "KRT19"),
    endothelial = c("PECAM1", "VWF", "EMCN", "ENG", "KDR", "FLT1", "RAMP2", "PLVAP"),
    pericyte = c("PDGFRB", "RGS5", "MCAM", "CSPG4", "NOTCH3", "ACTA2", "MYH11", "TAGLN"),
    t_nk = c("CD3D", "CD3E", "TRAC", "CD247", "NKG7", "GNLY", "KLRD1", "CCL5", "GZMB"),
    b_plasma = c("MS4A1", "CD79A", "CD74", "CD37", "CXCL13", "MZB1", "JCHAIN", "IGHG1", "SDC1"),
    myeloid = c("LST1", "TYROBP", "AIF1", "FCER1G", "C1QA", "C1QB", "CTSS", "SPP1", "CD68"),
    fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "FAP", "CXCL12", "TAGLN"),
    mast = c("TPSAB1", "TPSB2", "KIT", "CPA3", "MS4A2", "HDC")
  )
)

clean_symbol <- function(x) toupper(sub("[.][0-9]+$", "", trimws(as.character(x))))

map_breast_type <- function(x) {
  y <- tolower(as.character(x))
  out <- rep(NA_character_, length(y))
  out[grepl("cancer.*epithelial|malignant", y)] <- "malignant_epithelial"
  out[grepl("normal.*epithelial", y)] <- "normal_epithelial"
  out[grepl("(^|[^a-z])t[ _-]?cell|(^|[^a-z])nk|tcell", y)] <- "t_nk"
  out[grepl("b[ _-]?cell|plasma|plasmablast", y)] <- "b_plasma"
  out[grepl("myeloid|macro|monocyte|dendritic", y)] <- "myeloid"
  out[grepl("caf|fibro", y)] <- "fibroblast"
  out[grepl("endothelial", y)] <- "endothelial"
  out[grepl("pvl|perivascular|pericyte|smooth muscle", y)] <- "perivascular"
  out
}

map_kidney_type <- function(cell_association, state) {
  y <- tolower(as.character(cell_association))
  st <- tolower(as.character(state))
  out <- rep(NA_character_, length(y))
  out[grepl("malignant", y)] <- "malignant_epithelial"
  tubular <- grepl("tubular epithelial", y)
  out[tubular & !grepl("tumor", st)] <- "tubular_epithelial"
  out[tubular & grepl("tumor", st)] <- "malignant_epithelial"
  out[grepl("collecting duct", y)] <- "collecting_duct"
  out[grepl("endothelial", y)] <- "endothelial"
  out[grepl("pericyte", y)] <- "pericyte"
  out[grepl("(^|[^a-z])t cells|nk cells", y)] <- "t_nk"
  out[grepl("b/plasma|b cells|plasma", y)] <- "b_plasma"
  out[grepl("macrophage|myeloid|monocyte", y)] <- "myeloid"
  out[grepl("myofibroblast|fibroblast", y)] <- "fibroblast"
  out[grepl("mast", y)] <- "mast"
  out
}

read_table_any <- function(path, header = TRUE) {
  if (grepl("[.]gz$", path)) {
    fread(cmd = paste("gzip -dc", shQuote(path)), header = header, data.table = FALSE)
  } else {
    fread(path, header = header, data.table = FALSE)
  }
}

collapse_sparse_rows <- function(M, symbols) {
  symbols <- clean_symbol(symbols)
  ok <- !is.na(symbols) & nzchar(symbols)
  M <- M[ok, , drop = FALSE]
  symbols <- symbols[ok]
  if (anyDuplicated(symbols)) {
    levels_symbol <- unique(symbols)
    aggregate_rows <- sparseMatrix(
      i = seq_along(symbols),
      j = match(symbols, levels_symbol),
      x = 1,
      dims = c(length(symbols), length(levels_symbol))
    )
    M <- t(aggregate_rows) %*% M
    rownames(M) <- levels_symbol
  } else {
    rownames(M) <- symbols
  }
  as(M, "dgCMatrix")
}

group_pseudobulk_signature <- function(counts, donor, cell_type, top_n = 40L) {
  keep <- !is.na(donor) & nzchar(donor) & !is.na(cell_type) & nzchar(cell_type)
  counts <- counts[, keep, drop = FALSE]
  donor <- donor[keep]
  cell_type <- cell_type[keep]
  group <- paste(donor, cell_type, sep = "||")
  levels_group <- unique(group)
  indicator <- sparseMatrix(
    i = seq_along(group),
    j = match(group, levels_group),
    x = 1,
    dims = c(length(group), length(levels_group))
  )
  bulk <- counts %*% indicator
  lib <- pmax(Matrix::colSums(bulk), 1)
  bulk_cpm <- as.matrix(bulk %*% Diagonal(x = 1e6 / lib))
  group_type <- sub("^.*[|][|]", "", levels_group)
  types <- sort(unique(group_type))
  signature <- vapply(types, function(tt) rowMeans(bulk_cpm[, group_type == tt, drop = FALSE]), numeric(nrow(bulk_cpm)))
  rownames(signature) <- rownames(counts)
  colnames(signature) <- types

  eligible <- !grepl("^MT-|^RP[SL][0-9]|^HB[ABDG]", rownames(signature)) & rowMeans(signature) >= 1
  selected <- character(0)
  if (ncol(signature) > 1L) {
    for (jj in seq_len(ncol(signature))) {
      other <- rowMeans(signature[, -jj, drop = FALSE])
      specificity <- log1p(signature[, jj]) - log1p(other)
      idx <- which(eligible & is.finite(specificity))
      selected <- c(selected, rownames(signature)[idx[order(specificity[idx], decreasing = TRUE)[seq_len(min(top_n, length(idx)))]]])
    }
  }
  curated <- unique(unlist(canonical_markers$breast[intersect(types, names(canonical_markers$breast))], use.names = FALSE))
  selected <- unique(c(selected, intersect(curated, rownames(signature))))
  signature <- signature[selected, , drop = FALSE]
  list(
    signature = signature,
    group_table = data.frame(group = levels_group, donor = sub("[|][|].*$", "", levels_group), cell_type = group_type),
    cell_counts = as.data.frame(table(donor = donor, cell_type = cell_type), stringsAsFactors = FALSE)
  )
}

donor_balanced_dense_signature <- function(values, donor, cell_type) {
  keep <- !is.na(donor) & nzchar(donor) & !is.na(cell_type) & nzchar(cell_type)
  values <- values[, keep, drop = FALSE]
  donor <- donor[keep]
  cell_type <- cell_type[keep]
  group <- paste(donor, cell_type, sep = "||")
  group_means <- rowsum(t(values), group = group, reorder = FALSE) / as.numeric(table(group)[rownames(rowsum(t(values), group = group, reorder = FALSE))])
  group_means <- t(group_means)
  group_type <- sub("^.*[|][|]", "", colnames(group_means))
  types <- sort(unique(group_type))
  signature <- vapply(types, function(tt) rowMeans(group_means[, group_type == tt, drop = FALSE]), numeric(nrow(group_means)))
  rownames(signature) <- rownames(values)
  colnames(signature) <- types
  signature <- signature[apply(signature, 1L, function(z) all(is.finite(z)) && stats::sd(z) > 1e-8), , drop = FALSE]
  list(
    signature = signature,
    group_table = data.frame(group = colnames(group_means), donor = sub("[|][|].*$", "", colnames(group_means)), cell_type = group_type),
    cell_counts = as.data.frame(table(donor = donor, cell_type = cell_type), stringsAsFactors = FALSE)
  )
}
