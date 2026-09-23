#!/usr/bin/env Rscript

# Experiment-level normalization only. Never filter genes or combine studies here.
fail <- function(message) stop(message, call. = FALSE)

parse_args <- function(args) {
    required <- c("counts", "metadata", "qc", "annotation", "outdir")
    values <- list(smoke_test = FALSE)
    i <- 1L
    while (i <= length(args)) {
        key <- args[[i]]
        if (key == "--smoke-test") {
            if (values$smoke_test) fail("Duplicate --smoke-test")
            values$smoke_test <- TRUE
            i <- i + 1L
            next
        }
        if (!startsWith(key, "--")) fail(paste("Unexpected argument:", key))
        name <- substring(key, 3L)
        if (!(name %in% required) || !is.null(values[[name]])) {
            fail(paste("Unknown or duplicate option:", key))
        }
        if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
            fail(paste("Missing value for", key))
        }
        values[[name]] <- args[[i + 1L]]
        i <- i + 2L
    }
    missing <- required[vapply(required, function(x) is.null(values[[x]]), logical(1))]
    if (length(missing)) fail(paste("Missing options:", paste(missing, collapse = ", ")))
    for (name in required[required != "outdir"]) {
        if (!file.exists(values[[name]])) fail(paste("Missing input:", values[[name]]))
        values[[name]] <- normalizePath(values[[name]], mustWork = TRUE)
    }
    values$outdir <- normalizePath(values$outdir, mustWork = FALSE)
    values
}

sha256 <- function(path) {
    answer <- suppressWarnings(system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE))
    if (!length(answer) || !is.null(attr(answer, "status"))) {
        fail(paste("sha256sum failed for", path, paste(answer, collapse = " ")))
    }
    hash <- strsplit(answer[[1L]], "[[:space:]]+")[[1L]][[1L]]
    if (!grepl("^[0-9a-f]{64}$", hash)) fail(paste("Invalid SHA-256 for", path))
    hash
}

read_tsv <- function(path) {
    if (endsWith(path, ".gz")) {
        data.table::fread(cmd = paste("gzip -cd --", shQuote(path)),
                          sep = "\t", check.names = FALSE, integer64 = "double")
    } else {
        data.table::fread(path, sep = "\t", check.names = FALSE,
                          integer64 = "double")
    }
}

check_ids <- function(ids, label) {
    if (anyNA(ids) || any(!nzchar(ids))) fail(paste(label, "contains missing IDs"))
    if (anyDuplicated(ids)) fail(paste(label, "contains duplicate IDs"))
}

check_finite <- function(values, label, positive = FALSE) {
    if (anyNA(values) || any(!is.finite(values))) fail(paste(label, "contains NA/Inf"))
    if (positive && any(values <= 0)) fail(paste(label, "contains non-positive values"))
}

assert_nonempty_file <- function(path, label = basename(path)) {
    if (!file.exists(path)) fail(paste(label, "does not exist:", path))
    info <- file.info(path)
    size <- info$size[[1L]]
    if (is.na(size) || !is.finite(size) || size <= 0) {
        fail(paste(label, "is empty or has invalid size:", path))
    }
    invisible(size)
}

validate_gzip_file <- function(path, label = basename(path)) {
    assert_nonempty_file(path, label)
    answer <- suppressWarnings(system2(
        "gzip",
        c("-t", shQuote(path)),
        stdout = TRUE,
        stderr = TRUE
    ))
    status <- attr(answer, "status")
    if (!is.null(status) && status != 0L) {
        fail(paste(label, "failed gzip integrity validation:", paste(answer, collapse = " ")))
    }
    invisible(TRUE)
}

validate_json_file <- function(path, label = basename(path)) {
    assert_nonempty_file(path, label)
    value <- tryCatch(
        jsonlite::read_json(path, simplifyVector = FALSE),
        error = function(e) fail(paste(label, "is not valid JSON:", conditionMessage(e)))
    )
    if (!is.list(value)) fail(paste(label, "did not decode to a JSON object"))
    invisible(value)
}

validate_tsv_file <- function(path, expected_header, expected_rows, label) {
    assert_nonempty_file(path, label)
    con <- file(path, open = "rt")
    on.exit(close(con), add = TRUE)

    header <- readLines(con, n = 1L, warn = FALSE)
    if (length(header) != 1L || !identical(header[[1L]], expected_header)) {
        fail(paste(label, "has an unexpected header"))
    }

    rows <- 0L
    repeat {
        chunk <- readLines(con, n = 10000L, warn = FALSE)
        if (!length(chunk)) break
        rows <- rows + length(chunk)
    }

    if (rows != expected_rows) {
        fail(paste(label, "has", rows, "data rows; expected", expected_rows))
    }
    invisible(TRUE)
}

validate_matrix_gz_file <- function(path, genes, experiments, label) {
    validate_gzip_file(path, label)

    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)

    expected_header <- paste(c("Geneid", experiments), collapse = "\t")
    header <- readLines(con, n = 1L, warn = FALSE)

    if (length(header) != 1L || !identical(header[[1L]], expected_header)) {
        fail(paste(label, "has an unexpected matrix header"))
    }

    rows <- 0L
    repeat {
        chunk <- readLines(con, n = 5000L, warn = FALSE)
        if (!length(chunk)) break
        rows <- rows + length(chunk)
    }

    if (rows != length(genes)) {
        fail(paste(label, "has", rows, "gene rows; expected", length(genes)))
    }
    invisible(TRUE)
}

validate_checksum_manifest <- function(workdir) {
    manifest <- file.path(workdir, "SHA256SUMS.txt")
    assert_nonempty_file(manifest, "SHA256SUMS.txt")

    previous_wd <- getwd()
    on.exit(setwd(previous_wd), add = TRUE)
    setwd(workdir)

    answer <- suppressWarnings(system2(
        "sha256sum",
        c("--check", "SHA256SUMS.txt"),
        stdout = TRUE,
        stderr = TRUE
    ))
    status <- attr(answer, "status")

    if (!is.null(status) && status != 0L) {
        fail(paste("SHA256SUMS verification failed:", paste(answer, collapse = " ")))
    }

    invisible(TRUE)
}

write_matrix_gz <- function(values, genes, path) {
    con <- gzfile(path, open = "wb", compression = 6L)
    on.exit(close(con), add = TRUE)
    writeLines(paste(c("Geneid", colnames(values)), collapse = "\t"), con)
    # Bound the data-frame conversion to 128 rows rather than copying the matrix.
    for (start in seq.int(1L, nrow(values), by = 128L)) {
        rows <- start:min(start + 127L, nrow(values))
        block <- data.frame(Geneid = genes[rows], values[rows, , drop = FALSE],
                            check.names = FALSE)
        write.table(block, file = con, sep = "\t", row.names = FALSE,
                    col.names = FALSE, quote = FALSE, na = "NA")
    }
}

main <- function() {
    Sys.setenv(TZ = "UTC")
    args <- parse_args(commandArgs(trailingOnly = TRUE))
    for (package in c("edgeR", "data.table", "jsonlite")) {
        if (!requireNamespace(package, quietly = TRUE)) fail(paste("Missing R package:", package))
    }

    raw <- read_tsv(args$counts)
    if (ncol(raw) < 2L || names(raw)[[1L]] != "Geneid") fail("Matrix header must start with Geneid")
    genes <- raw[["Geneid"]]
    experiments <- names(raw)[-1L]
    check_ids(genes, "Matrix Geneid")
    check_ids(experiments, "Matrix experiments")
    n_genes <- length(genes)
    n_experiments <- length(experiments)
    if (!args$smoke_test && (n_genes != 27775L || n_experiments != 1473L)) {
        fail("Production matrix dimensions must be 27775 genes x 1473 experiments")
    }
    if (args$smoke_test && (n_genes < 2L || n_experiments < 2L ||
                            n_genes >= 27775L || n_experiments >= 1473L)) {
        fail("Smoke-test matrix must have at least two rows/columns and be a strict subset")
    }
    counts <- as.matrix(raw[, -1L, with = FALSE])
    rm(raw)
    if (!is.numeric(counts)) fail("Raw counts must be numeric")
    check_finite(counts, "Raw counts")
    if (any(counts < 0)) fail("Raw counts contain negative values")
    if (any(counts != trunc(counts))) fail("Raw counts are not integer-valued")
    if (any(counts >= 2^53)) fail("Raw counts reach the double-integer precision limit")
    dimnames(counts) <- list(genes, experiments)
    raw_library_size <- colSums(counts)
    check_finite(raw_library_size, "Raw library sizes", positive = TRUE)
    if (any(raw_library_size >= 2^53)) fail("Raw library size reaches the double-integer precision limit")
    n_nonzero <- rowSums(counts > 0)
    all_zero <- n_nonzero == 0L

    metadata <- read_tsv(args$metadata)
    if (!all(c("experiment_accession", "study_accession") %in% names(metadata))) {
        fail("Metadata lacks experiment_accession or study_accession")
    }
    check_ids(metadata$experiment_accession, "Metadata experiment_accession")
    if (nrow(metadata) != n_experiments ||
        !setequal(metadata$experiment_accession, experiments)) {
        fail("Matrix/metadata experiment sets or counts differ")
    }
    check_ids(unique(metadata$study_accession), "Metadata study_accession")
    if (!args$smoke_test && length(unique(metadata$study_accession)) != 92L) {
        fail("Production metadata must contain exactly 92 studies")
    }
    metadata <- metadata[match(experiments, metadata$experiment_accession)]

    annotation <- read_tsv(args$annotation)
    if (!all(c("Geneid", "Chr", "Start", "End", "Strand", "Length") %in% names(annotation))) {
        fail("Annotation lacks canonical featureCounts columns")
    }
    check_ids(annotation$Geneid, "Annotation Geneid")
    if (nrow(annotation) != n_genes || !setequal(annotation$Geneid, genes)) {
        fail("Annotation/matrix Geneid sets or counts differ")
    }
    if (!identical(annotation$Geneid, genes)) fail("Annotation/matrix Geneid order differs")
    lengths_bp <- annotation$Length
    if (!is.numeric(lengths_bp)) fail("Annotation Length must be numeric")
    check_finite(lengths_bp, "Annotation Length", positive = TRUE)

    qc <- read_tsv(args$qc)
    if (!all(c("experiment_accession", "total_assigned_gene_counts") %in% names(qc))) {
        fail("Raw-count QC lacks experiment_accession or total_assigned_gene_counts")
    }
    check_ids(qc$experiment_accession, "Raw-count QC experiment_accession")
    if (nrow(qc) != n_experiments || !setequal(qc$experiment_accession, experiments)) {
        fail("Matrix/raw-count QC experiment sets or counts differ")
    }
    qc_sizes <- qc$total_assigned_gene_counts[match(experiments, qc$experiment_accession)]
    if (!is.numeric(qc_sizes)) fail("Raw-count QC library sizes must be numeric")
    check_finite(qc_sizes, "Raw-count QC library sizes")
    if (!identical(as.numeric(raw_library_size), as.numeric(qc_sizes))) {
        fail("Matrix library sizes disagree with experiment_raw_count_qc.tsv")
    }
    rm(qc, annotation)

    y <- edgeR::DGEList(counts = counts, lib.size = raw_library_size)
    rm(counts)
    y <- edgeR::calcNormFactors(y, method = "TMM")
    factors <- y$samples$norm.factors
    effective_sizes <- raw_library_size * factors
    if (length(factors) != n_experiments) fail("Incorrect TMM factor count")
    check_finite(factors, "TMM factors", positive = TRUE)
    check_finite(effective_sizes, "Effective library sizes", positive = TRUE)

    logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2,
                         normalized.lib.sizes = TRUE)
    if (!identical(dim(logcpm), c(n_genes, n_experiments))) fail("Incorrect logCPM dimensions")
    check_finite(logcpm, "TMM logCPM")
    colnames(logcpm) <- experiments

    if (file.exists(args$outdir)) fail(paste("Output directory already exists:", args$outdir))
    parent <- dirname(args$outdir)
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    workdir <- tempfile(pattern = paste0(basename(args$outdir), ".tmp-"), tmpdir = parent)
    if (!dir.create(workdir)) fail(paste("Cannot create temporary output directory:", workdir))
    on.exit({
        # The source path disappears after a successful rename. Never touch outdir.
        if (dir.exists(workdir) && unlink(workdir, recursive = TRUE) != 0L) {
            warning(paste("Could not remove temporary output directory:", workdir))
        }
    }, add = TRUE)
    suffix <- paste0("_", n_experiments, ".tsv.gz")
    logcpm_name <- paste0("rnaseq_tmm_logcpm_by_experiment", suffix)
    tpm_name <- paste0("rnaseq_tpm_by_experiment", suffix)
    factor_name <- paste0("tmm_factors_", n_experiments, ".tsv")
    paths <- list(
        factors = file.path(workdir, factor_name),
        logcpm = file.path(workdir, logcpm_name),
        tpm = file.path(workdir, tpm_name),
        gene_qc = file.path(workdir, "gene_normalization_qc.tsv"),
        summary = file.path(workdir, "normalization_qc_summary.json"),
        provenance = file.path(workdir, "provenance.json")
    )
    data.table::fwrite(data.table::data.table(
        experiment_accession = experiments,
        study_accession = metadata$study_accession,
        raw_library_size = raw_library_size,
        tmm_norm_factor = factors,
        effective_library_size = effective_sizes
    ), paths$factors, sep = "\t")
    data.table::fwrite(data.table::data.table(
        Geneid = genes, Length = lengths_bp,
        all_zero_across_experiments = all_zero,
        n_nonzero_experiments = n_nonzero
    ), paths$gene_qc, sep = "\t")
    write_matrix_gz(logcpm, genes, paths$logcpm)
    rm(logcpm)
    gc(verbose = FALSE)

    # Use raw counts and canonical featureCounts Length; no TMM values enter TPM.
    tpm <- matrix(0, nrow = n_genes, ncol = n_experiments,
                  dimnames = list(genes, experiments))
    length_kb <- lengths_bp / 1000
    for (j in seq_len(n_experiments)) {
        rpk <- y$counts[, j] / length_kb
        denominator <- sum(rpk)
        if (!is.finite(denominator) || denominator <= 0) {
            fail(paste("Invalid TPM denominator for", experiments[[j]]))
        }
        tpm[, j] <- rpk / denominator * 1e6
    }
    if (!identical(dim(tpm), c(n_genes, n_experiments))) fail("Incorrect TPM dimensions")
    check_finite(tpm, "TPM")
    if (any(tpm < 0)) fail("TPM contains negative values")
    if (any(tpm[all_zero, , drop = FALSE] != 0)) fail("All-zero genes have nonzero TPM")
    tpm_sums <- colSums(tpm)
    tpm_tolerance <- 1e-4  # absolute TPM units; 1e-10 relative to 1,000,000
    if (any(abs(tpm_sums - 1e6) > tpm_tolerance)) {
        fail("TPM column sum outside absolute tolerance of 0.0001")
    }
    write_matrix_gz(tpm, genes, paths$tpm)
    rm(tpm)
    gc(verbose = FALSE)

    # Fail closed before any metadata/checksum publication.
    # A zero-byte or truncated output must never become a final result.
    validate_tsv_file(
        paths$factors,
        paste(c(
            "experiment_accession",
            "study_accession",
            "raw_library_size",
            "tmm_norm_factor",
            "effective_library_size"
        ), collapse = "\t"),
        n_experiments,
        "TMM factors"
    )
    validate_tsv_file(
        paths$gene_qc,
        paste(c(
            "Geneid",
            "Length",
            "all_zero_across_experiments",
            "n_nonzero_experiments"
        ), collapse = "\t"),
        n_genes,
        "gene normalization QC"
    )
    validate_matrix_gz_file(
        paths$logcpm,
        genes,
        experiments,
        "TMM logCPM matrix"
    )
    validate_matrix_gz_file(
        paths$tpm,
        genes,
        experiments,
        "TPM matrix"
    )

    range_summary <- function(x) list(min = min(x), median = median(x), max = max(x))
    dimensions <- list(genes = n_genes, experiments = n_experiments,
                       studies = length(unique(metadata$study_accession)))
    summary <- list(
        smoke_test = args$smoke_test,
        input_dimensions = dimensions,
        output_dimensions = list(tmm_logcpm = dimensions[1:2], tpm = dimensions[1:2]),
        globally_all_zero_genes = sum(all_zero),
        genes_nonzero_in_at_least_one_experiment = sum(!all_zero),
        tmm_factor = range_summary(factors),
        effective_library_size = range_summary(effective_sizes),
        tpm_column_sum = range_summary(tpm_sums),
        tpm_column_sum_target = 1e6,
        tpm_column_sum_absolute_tolerance = tpm_tolerance,
        all_zero_logcpm_note = "Finite logCPM from prior.count=2 does not establish expression. Use gene_normalization_qc.tsv."
    )
    jsonlite::write_json(summary, paths$summary, pretty = TRUE, auto_unbox = TRUE)
    validate_json_file(paths$summary, "normalization QC summary")

    input_paths <- list(matrix = args$counts, metadata = args$metadata,
                        raw_count_qc = args$qc, annotation = args$annotation)
    inputs <- lapply(input_paths, function(p) list(path = p, sha256 = sha256(p)))
    output_names <- c(factors = factor_name, logcpm = logcpm_name,
                      tpm = tpm_name, gene_qc = basename(paths$gene_qc),
                      summary = basename(paths$summary))
    output_hashes <- lapply(paths[names(output_names)], sha256)
    names(output_hashes) <- output_names
    git_commit <- suppressWarnings(system2("git", c("rev-parse", "HEAD"),
                                            stdout = TRUE, stderr = FALSE))
    if (!length(git_commit) || !is.null(attr(git_commit, "status"))) git_commit <- NA_character_
    script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    script_sha256 <- if (length(script_argument)) {
        sha256(normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE))
    } else {
        NA_character_
    }
    provenance <- list(
        smoke_test = args$smoke_test,
        inputs = inputs,
        output_sha256 = output_hashes,
        dimensions = dimensions,
        tmm_method = "TMM",
        logcpm_prior_count = 2,
        tpm_formula = "(raw_count / (featureCounts_Length_bp / 1000)) / column_sum_RPK * 1e6",
        tpm_length_source = args$annotation,
        gene_filtering = "none; all-zero genes retained",
        normalization_scope = "experiment-level only; no batch correction or study aggregation",
        r_version = R.version.string,
        package_versions = list(
            edgeR = as.character(utils::packageVersion("edgeR")),
            data.table = as.character(utils::packageVersion("data.table")),
            jsonlite = as.character(utils::packageVersion("jsonlite"))
        ),
        session_info = capture.output(utils::sessionInfo()),
        timestamp_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
        git_commit = git_commit[[1L]],
        normalization_script_sha256 = script_sha256
    )
    jsonlite::write_json(provenance, paths$provenance, pretty = TRUE, auto_unbox = TRUE)
    validate_json_file(paths$provenance, "provenance")

    checksum_names <- c(unname(output_names), basename(paths$provenance))
    checksum_paths <- file.path(workdir, checksum_names)
    checksum_lines <- vapply(seq_along(checksum_paths), function(i) {
        paste0(sha256(checksum_paths[[i]]), "  ", checksum_names[[i]])
    }, character(1))
    checksum_manifest <- file.path(workdir, "SHA256SUMS.txt")
    writeLines(checksum_lines, checksum_manifest)
    assert_nonempty_file(checksum_manifest, "SHA256SUMS.txt")
    validate_checksum_manifest(workdir)

    expected_final_names <- sort(c(checksum_names, "SHA256SUMS.txt"))
    actual_final_names <- sort(list.files(
        workdir,
        all.files = FALSE,
        full.names = FALSE,
        recursive = FALSE,
        no.. = TRUE
    ))
    if (!identical(actual_final_names, expected_final_names)) {
        fail(paste(
            "Unexpected final output set. Expected:",
            paste(expected_final_names, collapse = ","),
            "Observed:",
            paste(actual_final_names, collapse = ",")
        ))
    }

    if (!file.rename(workdir, args$outdir)) {
        fail(paste("Cannot publish output directory:", args$outdir))
    }
    cat("Validated and wrote", n_genes, "genes x", n_experiments,
        "experiments to", args$outdir, "\n")
}

tryCatch(main(), error = function(e) {
    message("ERROR: ", conditionMessage(e))
    quit(save = "no", status = 1L)
})
