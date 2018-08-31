#' Load a Top Annotations File
#' 
#' This function loads the Top Annotations file from the working directory. If a file name is not specified, it assumes
#' that the file is named "gwas.top.annot.txt".
#' 
#' @param filename The name of your Top Annotations file in the working directory. Defaults to NULL.
#' @export
loadTA <- function(filename = NULL) {
    if (is.null(filename)) {
        read.table("gwas.top.annot.txt", header = TRUE, sep = " ")
    } else {
        read.table(filename, header = TRUE, sep = " ")
    }
}



#' Review the Info for Gene Hits
#' 
#' This function converts the messy "Gene Annotation" column into a more easily-readible summary. It generates a list
#' of data frames. The first data frame contains the FBgn ID, Gene Symbol, Site Class, and Distance to Gene for all
#' significant polymorphisms in or near a given gene. If a polymorphism is located in or near more than one gene,
#' the other possible hits are noted in the subsequent data frames.
#' 
#' @param TA Your Top Annotations data.
#' @return A list of data frames summarizing the info for your gene hits.
#' @export

geneInfo <- function(TA) {
    colsplit <- reshape2::colsplit
    str_count <- stringr::str_count
    gene.annot <- colsplit(TA$GeneAnnotation, pattern = ",", names = c("SiteClass", "TranscriptAnnot"))
    glen <- max(str_count(gene.annot$SiteClass, ";")) + 1
    genes <- colsplit(gene.annot$SiteClass, pattern = ";", names = paste("FBgn", c(1:glen), sep="_"))
    genes[genes == ""] <- NA
    genes$FBgn_1 <- gsub("SiteClass", "", genes$FBgn_1)
    genes <- apply(genes, 2, function(x) {
        gsub("\\[|\\]", "", x)
        })
    genes[genes == "|||"] <- NA
    split.genes <- lapply(split(t(genes), 1:glen), as.data.frame)
    gnames <- c("FBgn_ID", "Gene_Symbol", "Site_class", "Distance_to_Gene")
    split.genes <- lapply(split.genes, function(x) {
        colsplit(x[ ,1], pattern = "\\|", names = gnames)
        })
    split.genes[[1]][which(split.genes[[1]]$Site_class == ""), 3] <- "INTERGENIC"
    split.genes[[1]][which(split.genes[[1]]$FBgn_ID == ""), 1:2] <- NA
    return(split.genes)
}



#' List All the Unique Gene Hits
#' 
#' This function lists all the unqiue genes that had a significantly associated polymorphism. If a polymorphism was
#' in or near more than one gene, the alternate genes are also included in this list.
#' 
#' @param TA Your Top Annotations data.
#' @return A character vector containing all unique gene hits.
#' @export
uniqGenes <- function(TA) {
    y <- geneInfo(TA)
    uniq.genes <- lapply(y, function(z) {
        unique(z$Gene_Symbol)
        })
    uniq.genes <- unique(unlist(uniq.genes))
    uniq.genes <- uniq.genes[uniq.genes != "" & !is.na(uniq.genes)]
    return(uniq.genes)
}



#' Review the Info for Transcript Hits
#' 
#' This function converts the messy "Transcript Annotation" into a more easily readable summary.
#' 
#' @param TA Your Top Annotations Data.
#' @return A list of data frames containing a summary of transcript info.
#' @export
transcrInfo <- function(TA) {
    colsplit <- reshape2::colsplit
    str_count <- stringr::str_count
    gene.annot <- colsplit(TA$GeneAnnotation, pattern = ",", names = c("SiteClass", "TranscriptAnnot"))
    slen <- max(str_count(gene.annot$TranscriptAnnot, ";")) + 1
    trans <- colsplit(gene.annot$TranscriptAnnot, pattern = ";", names = paste("Transcript", c(1:slen), sep = "_"))
    trans[trans == ""] <- NA
    trans$Transcript_1 <- gsub("TranscriptAnnot", "", trans$Transcript_1)
    trans <- apply(trans, 2, function(x) {
        gsub("\\[|\\]|)", "", x)
        })
    split.trans <- lapply(split(t(trans), 1:ncol(trans)), as.data.frame)
    tnames <- c("Effect_Impact", "Functional_Class", "Codon_Change", "Amino_Acid_Change", "Amino_Acid_Length",
                "Gene_Name", "Gene_BioType", "Coding", "Transcript", "Exon", "Errors", "Warnings")
    split.trans <- lapply(split.trans, function(x) {
        colsplit(x[ , 1], pattern = "\\|", names=tnames)
        })
    return(split.trans)
}



#' Review the Regulation Info for Transcript Hits
#' 
#' This function converts the messy "Regulation Information" column into a more easily readable summary.
#' 
#' @param TA Your Top Annotation Data.
#' @return A list of data frames containing a summary of regulatory info.
#' @export
regInfo <- function(TA) {
    colsplit <- reshape2::colsplit
    str_count <- stringr::str_count
    rlen <- max(str_count(TA$RegulationAnnotation, ",")) + 1
    reg <- colsplit(TA$RegulationAnnotation, pattern = ",", names = paste("Reg_Feature", c(1:rlen), sep = "_"))
    reg$Reg_Feature_1 <- sub("-", "", reg$Reg_Feature_1)
    reg <- apply(reg, 2, function(x) {
        gsub("\\(|\\)", "", x)
        })
    split.reg <- lapply(split(t(reg), 1:ncol(reg)), as.data.frame)
    rnames <- c("Type", "Source", "Flybase_ID")
    split.reg <- lapply(split.reg, function(x) {
        colsplit(x[ , 1], pattern = "\\|", names = rnames)
        })
    return(split.reg)
}



#' Summarize Top Annotations
#' 
#' This function provides a full summary of all information contained in the Top Annotations File, including Gene Info,
#' Unique Genes, Transcript Info, andd Regulation Info. If a file name is specified, this summary is exported to an
#' Excel file.
#' 
#' @param TA Your Top Annotation Data
#' @param filename The name of the exported Excel file. Defaults to NULL.
#' @return A list of data frames (or an Excel file with multiple sheets) containing a summaray of Top Annotations.
#' @export
summarizeTA <- function(TA, filename=NULL) {
    setNames <- stats::setNames
    gnames <- c("FBgn_ID", "Gene_Symbol", "Site_class", "Distance_to_Gene")
    tnames <- c("Effect_Impact", "Functional_Class", "Codon_Change", "Amino_Acid_Change", "Amino_Acid_Length",
                "Gene_Name", "Gene_BioType", "Coding", "Transcript", "Exon", "Errors", "Warnings")
    rnames <- c("Type", "Source", "Flybase_ID")
    snp.info <- TA[ , 1:10]
    all.genes <- setNames(data.frame(matrix(unlist(geneInfo(TA)), nrow = nrow(TA), byrow = FALSE)),
                                 rep(gnames, length(geneInfo(TA))))
    uniq.genes <- data.frame("Unique_Genes" = uniqGenes(TA))
    all.trans <- setNames(data.frame(matrix(unlist(transcrInfo(TA)), nrow = nrow(TA), byrow = FALSE)),
                                 rep(tnames, length(transcrInfo(TA))))
    all.reg <- setNames(data.frame(matrix(unlist(regInfo(TA)), nrow = nrow(TA), byrow = FALSE)),
                               rep(rnames, length(regInfo(TA))))
    combined <- list("All Annotations" = TA, "SNP Info" = snp.info, "Genes Info" = all.genes,
                     "Unique Genes" = uniq.genes, "Transcripts Info" = all.trans, "Regulation Info" = all.reg)
    if (is.null(filename)) {
        return(combined)
    } else {
        (openxlsx::write.xlsx(combined, file = filename,  colNames = TRUE, rowNames = FALSE, keepNA = TRUE,
                              colWidths = "auto"))
    }
}