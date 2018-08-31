#' Load SNP Calls
#'
#' This function loads the SNP Calls file from the working directory. If no file name is specified,
#' it assumes the file is named "snp_calls.csv".
#'
#' @param filename The name of your SNP Calls file. Defaults to NULL.
#' @export
loadCalls <- function(filename = NULL) {
    if (is.null(filename)) {
        read.table("snp_calls.csv", header = TRUE, sep = ",")
    } else {
        read.table(filename, header = TRUE, sep = ",")
    }
}



#' Count How DGRP Lines Have the SNP Variant
#'
#' This function shows in how many DGRP lines a given SNP variant is present, absent, or unknown.
#'
#' @param snp Your SNP Calls data.
#' @return A data frame displaying how many lines a given SNP variant is present, absent, or unknown.
#' @export
countSNP <- function(snp) {
    rownames_to_column <- tibble::rownames_to_column
    snp <- snp[,-1]
    present <- sapply(snp, function(y) length(which(y != 0 & !is.na(y))))
    absent <- sapply(snp, function(y) length(which(y == 0 & !is.na(y))))
    na <- sapply(snp, function(y) sum(is.na(y)))
    counts <- data.frame("Present" = present, "Absent" = absent, "N/A" = na)
    counts <- rownames_to_column(counts, var = "ID")
    return(counts)
}



#' Display Which DGRP Lines Have the SNP Variant
#'
#' This function lists all the DGRP lines where a given SNP variant is present, absent, or unknown.
#'
#' @param snp Your SNP Calls data.
#' @return A list of vectors showing the DGRP lines where a given SNP variant is present, absent, or unknown.
#' @export
snpLines <- function(snp) {
    `%>%` <- dplyr::`%>%`
    data <- snp
    countSnps <- function(x) {
        present <- which(x != 0 & !is.na(x)) %>% data$Line[.]
        absent <- which(x == 0 & !is.na(x)) %>% data$Line[.]
        na <- which(is.na(x)) %>% data$Line[.]
        lines <- list("Present" = present, "Absent" = absent, "N/A" = na)
        return(lines)
    }
    all <- lapply(snp[ , -1], countSnps)
    return(all)
}



#' List which DGRP Lines have a SNP and a Binary Phenotype.
#'
#' This function shows whether or not each DGRP line has a SNP listed in the SNP Calls document, as well as whether or not
#' each DGRP line has a specified binary phenotype.
#'
#' @param snp Your SNP Calls data.
#' @param phen A data table where column 1 lists all DGRP lines in ascending numerical order, and column 2 shows the presence
#' or absence of a binary phenotype coded as 1 or 0.
#' @return A data table showing whether each DGRP line has a given SNP ("Ynsp" or "Nsnp") and whether it has the binary phenotye
#' ("Yphen, Nphen".)
#' @export
snpPhen <- function(snp, phen){
    calls <- snp
    calls$Line <- as.numeric(gsub("line_", "", calls$Line))
    calls <- calls[order(calls$Line), ]
    comp <- cbind(Line = phen[ , 1], Phenotype = phen[ , 2], calls[ , -1])
    for (i in 3:ncol(comp)) {
        comp[ , i] <- ifelse(is.na(comp[ , i]), "NA", paste(
            ifelse(comp[ , i] > 0, "Ysnp", "Nsnp"),
            ifelse(comp[ , 2] == 1, "Yphen", "Nphen"), sep = ", "));
    }
    return(comp)
}



#' Generate 2x2 Contingency Tables for a Binary Phenotype
#'
#' This function creates 2x2 contingency tables for eacn SNP and a binary phenotype. With the default setting of table = FALSE,
#' the resutls are displayed vertically in a single data table. For table = TRUE, the function will instead generate a list of
#' 2x2 contingency tables.
#'
#' @param snp Your SNP Calls data.
#' @param phen A data table where column 1 lists all DGRP lines in ascending numerical order, and column 2 shows the presence
#' or absence of a binary phenotype coded as 1 or 0.
#' @param table Indicates the formatting for your results (single data frame or list of 2x2 contingency tables). Defaults to FALSE.
#' @return Contingency table(s) for all SNPs and a binary phenotype.
#' @export
phenFreq <- function(snp, phen, table = FALSE) {
    comp <- snpPhen(snp, phen)
    countmatrix <- function(x) {
        YY <- length(grep("Ysnp, Yphen", x))
        YN <- length(grep("Ysnp, Nphen", x))
        NY <- length(grep("Nsnp, Yphen", x))
        NN <- length(grep("Nsnp, Nphen", x))
        counts <- matrix(c(YY, YN, NY, NN), ncol = 2, byrow = TRUE)
        colnames(counts) <- c("Phen+", "Phen-")
        rownames(counts) <- c("SNP+", "SNP-")
        return(counts)
    }
    countdf <- function(x) {
        YY <- length(grep("Ysnp, Yphen", x))
        YN <- length(grep("Ysnp, Nphen", x))
        NY <- length(grep("Nsnp, Yphen", x))
        NN <- length(grep("Nsnp, Nphen", x))
        counts <- c(YY, YN, NY, NN)
        return(counts)
    }
    if (table == TRUE) {
        allcounts <- lapply(comp[ , -(1:2)], countmatrix)
        return(allcounts)
    } else {
        allcounts <- as.data.frame(sapply(comp[ , -(1:2)], countdf))
        colnames(allcounts) <- colnames(comp[ , -(1:2)])
        rownames(allcounts) <- c("SNP+, Phen+", "SNP+, Phen-", "SNP-, Phen+", "SNP-, Phen-")
        allcounts <- tibble::rownames_to_column(allcounts, var = "Result")
        return(allcounts)
    }
}



#' Compute Odds Ratios and p-Values Based on Contingency Tables
#'
#' This function computes the odds ratio for each SNP and a binary phenotype. It also shows the p-value based on a Fisher's
#' exact test. Note this p-value is just for approximation and does not correct for multiple testing.
#'
#' @param snp Your SNP Calls data.
#' @param phen A data table where column 1 lists all DGRP lines in ascending numerical order, and column 2 shows the presence
#' or absence of a binary phenotype coded as 1 or 0.
#' @return A data table showing the odds ratio and p-value for each SNP.
#' @export
phenOdds <- function(snp, phen) {
    fisher.test <- stats::fisher.test
    rownames_to_column <- tibble::rownames_to_column
    freq <- phenFreq(snp, phen, table = TRUE)
    fisher <- lapply(freq, fisher.test)
    results <- t(data.frame(sapply(fisher, function(x) c("Odds Ratio" = x[3], "p-value" = x[1]))))
    results <- rownames_to_column(as.data.frame(results))
    colnames(results) <- c("ID", "Odds Ratio", "p-value")
    return(results)
}



#' Summarize SNP Calls Data
#'
#' This function generates a complete summary of all information contained in the SNP Calls data, including "Calls" and "SNP
#' Counts". If a binary phenotype data is provided, the summary will also include information for "SNP vs. Phenotype Lines",
#' "SNP vs. Phenotype Total", and "Odds Ratio". If a file name is specified, this summary is exported to an Excel file.
#'
#' @param snp Your SNP Calls data.
#' @param phen A data table where column 1 lists all DGRP lines in ascending numerical order, and column 2 shows the presence
#' or absence of a binary phenotype coded as 1 or 0. Defaults to NULL.
#' @param filename The name of the exported Excel file. Defaults to NULL.
#' @return A list of data frames (or an Excel file with multiple sheets) containing a summaray of SNP Calls.
#' @export
summarizeCalls <- function(snp, phen = NULL, filename = NULL) {
    if (is.null(phen)) {
        combined <- list("Calls" = snp, "SNP Counts" = countSNP(snp))
    } else {
        combined <- list("Calls" = snp, "SNP Counts" = countSNP(snp), "SNP vs. Phenotype Lines" = snpPhen(snp, phen),
                     "SNP vs. Phenotype Total" = phenFreq(snp, phen), "Odds Ratio" = phenOdds(snp, phen))
    }
    if (is.null(filename)) {
        return(combined)
    } else {
        openxlsx::write.xlsx(combined, file = filename, colNames = TRUE, rowNames = FALSE, colWidths = "auto")
    }
}
