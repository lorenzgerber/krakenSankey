main <- function() {
  # Producing a Sanky Plot from a Kraken report file.
  #
  # Takes one argument:
  # The file name of the kraken report.
  #
  # Ex. usage:
  #   Rscript krakenSanky.R kraken.report
  #
  library(magrittr)
  library(networkD3)
  
  # get input args
  args <- commandArgs(trailingOnly = TRUE)
  filepath <- toString(args[1])
  
  # Load data
  report <- read_report(filepath)
  
  # Construct filename
  basename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(filepath))
  output_file <- paste(basename,'.html', sep = '')
  
  # Output to png file
  build_sankey_network(report) %>% saveNetwork(file = output_file)
}




#' Read Kraken-style and MetaPhlAn reports
#'
#' @param myfile Kraken-style or MetaPhlAn report file.
#' @param has_header If the kraken report has a header or not.
#' @param check_file If TRUE, only the first 5 lines of the file are loaded.
#'
#' @return report data.frame
#' @export
#'
read_report <- function(myfile, has_header=NULL, check_file = FALSE) {

  # TODO: Support for gzipped files ..
  #myfile <- file(myfile)
  #file_class <- summary(myfile)$class
  #if (file_class == "gzfile")
  #  myfile <- gzcon(myfile)

  first.line <- tryCatch( readLines(myfile,n=1, warn=FALSE),
                          error = function(e) { warning("Error reading ",myfile); return() })
  isASCII <-  function(txt) {
    if (length(txt) == 0)
      return(FALSE)
    raw <- charToRaw(txt)
    all(raw <= as.raw(127) && (raw >= as.raw(32) | raw == as.raw(9)))
  }
  if (length(first.line) == 0) {
    dmessage("Could not read ", myfile, ".")
    return(NULL)
  }
  tryCatch({
    if (nchar(first.line) == 0) {
      dmessage("First line of ", myfile, " is empty")
      return(NULL)
    }
  }, error = function(e) {
    dmessage(e)
    return(NULL)
  })

  if (!isTRUE(isASCII(first.line))) {
    dmessage(myfile," is not a ASCII file")
    return(NULL)
  }

  if (is.null(has_header)) {
    has_header <- grepl("^[a-zA-Z#%\"]",first.line)
  }

  is_krakenu_fmt <- grepl("^.?%\treads\ttaxReads\tkmers", first.line)
  is_kaiju_fmt <- grepl("^  *%\t  *reads", first.line)
  nrows <- ifelse(isTRUE(check_file), 5, -1)
  if (!is_krakenu_fmt && is_kaiju_fmt) {
    cont <- readLines(myfile)
    cont <- cont[!grepl("^-", cont)]
    cont <- sub(".*\t  *","", cont)
    cont <- sub("; ?$","", cont)
    report <- utils::read.delim(textConnection(cont), stringsAsFactors = FALSE)
    colnames(report) <- c("taxonReads", "taxLineage")
    report$cladeReads <- report$taxonReads

    report$taxLineage <- gsub("^","-_",report$taxLineage)
    report$taxLineage <- gsub("; ","|-_",report$taxLineage)
    report$taxLineage <- gsub("-_Viruses", "d_Viruses", report$taxLineage, fixed=T)
    report$taxLineage <- gsub("-_cellular organisms|-_Bacteria", "-_cellular organisms|d_Bacteria", report$taxLineage, fixed=T)
    report$taxLineage <- gsub("-_cellular organisms|-_Eukaryota", "-_cellular organisms|d_Eukaryota", report$taxLineage, fixed=T)
    report$taxLineage <- gsub("-_cellular organisms|-_Archaea", "-_cellular organisms|d_Archaea", report$taxLineage, fixed=T)
    report$taxLineage[1:(length(report$taxLineage)-1)] <- paste0("-_root|", report$taxLineage[1:(length(report$taxLineage)-1)])

    report$taxLineage[report$taxLineage=="-_unclassified"] <- "u_unclassified"

    new_counts <- integer(length = 0)
    for (j in seq_len(nrow(report))) {
      count <- report$cladeReads[j]
      tl <- report$taxLineage[j]
      tl2 <- sub("\\|[^|]*$","", tl)
      while (tl2 != tl) {
        if (tl2 %in% names(new_counts)) {
          new_counts[tl2] <- new_counts[tl2] + count
        } else {
          new_counts[tl2] <- count
        }
        tl <- tl2
        tl2 <- sub("\\|[^|]*$","", tl)
      }
    }
    report <- rbind(report,
                    data.frame(taxonReads=0,taxLineage=names(new_counts),cladeReads=as.integer(new_counts)))
    tl_order <- order(report$taxLineage)
    tl_order <- c(tl_order[length(tl_order)],tl_order[-length(tl_order)])
    report <- report[tl_order, c("taxLineage", "taxonReads", "cladeReads")]
  } else if (has_header) {
    report <- tryCatch({
      utils::read.table(myfile,sep="\t",header = T,
                        quote = "",stringsAsFactors=FALSE,
                        comment.char = "#", nrows = nrows,
                        check.names=FALSE)
    }, error = function(x) NULL, warning = function(x) NULL)
    if (is.null(report)) { return(NULL); }
    #colnames(report) <- c("percentage","cladeReads","taxonReads","taxRank","taxID","n_unique_kmers","n_kmers","perc_uniq_kmers","name")

    ## harmonize column names. TODO: Harmonize them in the scripts!
    colnames(report)[colnames(report) %in% c("#%","%","clade_perc","perc","percReadsClade")] <- "percentage"
    colnames(report)[colnames(report) %in% c("reads","numReadsClade","n_reads_clade","n.clade","n-clade")] <- "cladeReads"
    colnames(report)[colnames(report) %in% c("taxReads","numReadsTaxon","n_reads_taxo","n.stay","n-stay")] <- "taxonReads"
    colnames(report)[colnames(report) %in% c("rank","tax_taxRank","level")] <- "taxRank"
    colnames(report)[colnames(report) %in% c("tax","taxonid")] <- "taxID"
    colnames(report)[colnames(report) %in% c("indentedName","taxName")] <- "name"
    colnames(report)[colnames(report) %in% c("dup")] <- "kmerDuplicity"
    colnames(report)[colnames(report) %in% c("cov")] <- "kmerCoverage"
  } else {
    report <- tryCatch({
      utils::read.table(myfile,sep="\t",header = F,
                        col.names = c("percentage","cladeReads","taxonReads","taxRank","taxID","name"),
                        quote = "",stringsAsFactors=FALSE,
                        nrows = nrows)
    }, error=function(x) NULL, warning=function(x) NULL)
    if (is.null(report)) { return(NULL); }
  }

  if (ncol(report) < 2) {
    return(NULL)
  }
  if (colnames(report)[2] == "Metaphlan2_Analysis") {
    ## Metaphlan report
    colnames(report) <- c("taxLineage", "cladeReads")
    report <- report[order(report$taxLineage), ]
    report$taxLineage <- gsub("_"," ",report$taxLineage)
    report$taxLineage <- gsub("  ","_",report$taxLineage)
    report$taxLineage <- paste0("-_root|", report$taxLineage)

    report <- rbind(
      data.frame(taxLineage=c("u_unclassified","-_root"),"cladeReads"=c(0,100), stringsAsFactors = F),
      report)
  }

  if (all(c("name","taxRank") %in% colnames(report)) && !"taxLineage" %in% colnames(report)) {
    ## Kraken report
    report$depth <- nchar(gsub("\\S.*","",report$name))/2
    if (!all(report$depth == floor(report$depth))) {
      warning("Depth doesn't work out!")
      return(NULL)
    }
    report$name <- gsub("^ *","",report$name)

    ## 'fix' taxRank
    table(report$taxRank)
    allowed_taxRanks <- c("U", "S", "G", "F", "C", "D", "O", "K", "P")
    report$taxRank[report$taxRank=="class"] <- "C"
    report$taxRank[report$taxRank=="family"] <- "F"
    report$taxRank[report$taxRank=="genus"] <- "G"
    report$taxRank[report$taxRank=="superkingdom"] <- "D"
    report$taxRank[report$taxRank=="kingdom"] <- "K"
    report$taxRank[report$taxRank=="order"] <- "O"
    report$taxRank[report$taxRank=="phylum"] <- "P"
    report$taxRank[report$taxRank=="species"] <- "S"
    report$taxRank[report$name=="unclassified"] <- "U"
    report$taxRank[!report$taxRank %in% allowed_taxRanks] <- "-"

    report$name <- paste(tolower(report$taxRank),report$name,sep="_")

    rownames(report) <- NULL

    ## make taxLineage path
    report$taxLineage <- report$name
    n <- nrow(report)
    depths <- report$depth
    taxLineages <- report$name
    taxLineages_p <- as.list(seq_along(report$name))

    depth_row_tmp <- c(1:25)

    for (current_row in seq(from=1, to=nrow(report))) {
      dcr <- depths[current_row]
      depth_row_tmp[dcr+1] <- current_row
      if (dcr >= 1) {
        prev_pos <- depth_row_tmp[[dcr]]
        taxLineages_p[[current_row]] <- c(taxLineages_p[[prev_pos]], current_row)
      }
    }
    report$taxLineage <- sapply(taxLineages_p, function(x) paste0(taxLineages[x], collapse="|"))
    #report$taxLineage <- taxLineages

  } else if ("taxLineage" %in% colnames(report)) {
    taxLineages <- strsplit(report$taxLineage, "|", fixed=TRUE)

    if (!"name" %in% colnames(report))
      report$name <- sapply(taxLineages, function(x) x[length(x)])

    if (!"depth" %in% colnames(report)) {
      report$depth <- sapply(taxLineages, length) - 1
    }
    if (!"taxRank" %in% colnames(report))
      report$taxRank <- toupper(substr(report$name, 0, 1))
  }


  if (!all(c("name","taxRank") %in% colnames(report)) ||
      nrow(report) < 2 ||
      !((report[1,"name"] == "u_unclassified" && report[2,"name"] == "-_root") || report[1,"name"] == "-_root")) {
    message(paste("Warning: File",myfile,"does not have the required format"))
    print(utils::head(report))
    return(NULL)
  }


  if (!"taxonReads" %in% colnames(report)) {
    parent <- sub("^\\(.*\\)\\|.*$", "\\1", report$taxLineage)
    taxLineages <- strsplit(report$taxLineage, "|", fixed=TRUE)
    ## fix taxonReads
    report$taxonReads <- report$cladeReads - sapply(report$name, function(x) sum(report$cladeReads[parent == x]))
    #report$taxonReads[sapply(report$taxonReads, function(x) isTRUE(all.equal(x, 0)))] <- 0
    report$taxonReads[report$taxonReads <= 0.00001] <- 0  # fix for rounding in percentages by MetaPhlAn
  }

  report$percentage <- signif(report$cladeReads/sum(report$taxonReads),6) * 100
  if ('n_unique_kmers'  %in% colnames(report))
    report$kmerpercentage <- round(report$n_unique_kmers/sum(report$n_unique_kmers,na.rm=T),6) * 100
  #report$taxRankperc <- 100/taxRank(report$cladeReads)

  #report$depth <- NULL

  if ("taxID" %in% colnames(report)) {
    std_colnames <- c("percentage","cladeReads","taxonReads","taxRank", "taxID","name")
  } else {
    std_colnames <- c("percentage","cladeReads","taxonReads","taxRank","name")
  }
  stopifnot(all(std_colnames %in% colnames(report)))
  report[, c(std_colnames, setdiff(colnames(report), std_colnames))]
}



#' Read kraken or centrifuge-style report
#'
#' @param myfile kraken report file
#' @param collapse  should the results be collapsed to only those taxRanks specified in keep_taxRanks?
#' @param keep_taxRanks taxRanks to keep when collapse is TRUE
#' @param min.depth minimum depth
#' @param filter_taxon filter certain taxon names
#' @param has_header if the kraken report has a header or not
#' @param add_taxRank_columns if TRUE, for each taxRank columns are added
#'
#' @return report data.frame
#' @export
#'
read_report2 <- function(myfile,collapse=TRUE,keep_taxRanks=c("D","K","P","C","O","F","G","S"),min.depth=0,filter_taxon=NULL,
                         has_header=NULL,add_taxRank_columns=FALSE) {

  first.line <- readLines(myfile,n=1)
  isASCII <-  function(txt) all(charToRaw(txt) <= as.raw(127))
  if (!isASCII(first.line)) {
    dmessage(myfile," is no valid report - not all characters are ASCII")
    return(NULL)
  }
  if (is.null(has_header)) {
    has_header <- grepl("^[a-zA-Z]",first.line)
  }

  if (has_header) {
    report <- utils::read.table(myfile,sep="\t",header = T,
                                quote = "",stringsAsFactors=FALSE, comment.char="#")
    #colnames(report) <- c("percentage","cladeReads","taxonReads","taxRank","taxID","n_unique_kmers","n_kmers","perc_uniq_kmers","name")

    ## harmonize column names. TODO: Harmonize them in the scripts!
    colnames(report)[colnames(report)=="clade_perc"] <- "percentage"
    colnames(report)[colnames(report)=="perc"] <- "percentage"

    colnames(report)[colnames(report)=="n_reads_clade"] <- "cladeReads"
    colnames(report)[colnames(report)=="n.clade"] <- "cladeReads"

    colnames(report)[colnames(report)=="n_reads_taxo"] <- "taxonReads"
    colnames(report)[colnames(report)=="n.stay"] <- "taxonReads"

    colnames(report)[colnames(report)=="rank"] <- "taxRank"
    colnames(report)[colnames(report)=="tax_rank"] <- "taxRank"

    colnames(report)[colnames(report)=="taxonid"] <- "taxID"
    colnames(report)[colnames(report)=="tax"] <- "taxID"

  } else {
    report <- utils::read.table(myfile,sep="\t",header = F,
                                col.names = c("percentage","cladeReads","taxonReads","taxRank","taxID","name"),
                                quote = "",stringsAsFactors=FALSE, comment.char="#")
  }

  report$depth <- nchar(gsub("\\S.*","",report$name))/2
  report$name <- gsub("^ *","",report$name)
  report$name <- paste(tolower(report$taxRank),report$name,sep="_")


  ## Only stop at certain taxRanks
  ## filter taxon and further up the tree if 'filter_taxon' is defined
  kraken.tree <- build_kraken_tree(report)
  report <- collapse.taxRanks(kraken.tree,keep_taxRanks=keep_taxRanks,filter_taxon=filter_taxon)

  ## Add a metaphlan-style taxon string
  if (add_taxRank_columns) {
    report[,keep_taxRanks] <- NA
  }
  report$taxLineage = report$name
  rows_to_consider <- rep(FALSE,nrow(report))

  for (i in seq_len(nrow(report))) {
    ## depth > 2 correspond to taxRanks below 'D'
    if (i > 1 && report[i,"depth"] > min.depth) {
      ## find the maximal index of a row below the current depth
      idx <- report$depth < report[i,"depth"] & rows_to_consider
      if (!any(idx)) { next() }

      current.taxRank <- report[i,'taxRank']
      my_row <- max(which(idx))
      report[i,'taxLineage'] <- paste(report[my_row,'taxLineage'],report[i,'taxLineage'],sep="|")

      if (add_taxRank_columns) {
        if (report[my_row,'taxRank'] %in% keep_taxRanks) {
          taxRanks.cp <- keep_taxRanks[seq(from=1,to=which(keep_taxRanks == report[my_row,'taxRank']))]
          report[i,taxRanks.cp] <- report[my_row,taxRanks.cp]
        }

        report[i,report[i,'taxRank']] <- report[i,'name']
      }
    }
    
    rows_to_consider[i] <- TRUE
  }

  report <- report[report$depth >= min.depth,]

  report$percentage <- round(report$cladeReads/sum(report$taxonReads),6) * 100

  for (column in c("taxonReads", "cladeReads"))
    if (all(floor(report[[column]]) == report[[column]]))
      report[[column]] <- as.integer(report[[column]])

  if ('n_unique_kmers'  %in% colnames(report))
    report$kmerpercentage <- round(report$n_unique_kmers/sum(report$n_unique_kmers,na.rm=T),6) * 100
  #report$taxRankperc <- 100/taxRank(report$cladeReads)

  rownames(report) <- NULL

  report
}




#' Filter lines from a kraken report result based on the taxonomy name
#'
#' It updates the read_stay counts, and removes any children below the
#' entry, and any parent entries that have no "cladeReads" that stay
#'
#' @param report Report \code{data.frame}.
#' @param filter_taxon Name of entry to remove.
#' @param rm_clade If \code{TRUE}, remove all cladeReads at and below clade, otherwise just set the number of cladeReads that stay at taxon to zero.
#' @param do_message If \code{TRUE}, report how many rows and cladeReads were deleted.
#'
#' @return filtered report
#' @export
filter_taxon <- function(report, filter_taxon, rm_clade = TRUE, do_message=FALSE) {
  taxon_depth <- NULL
  taxonReads <- 0

  pos.taxons <- which(sub("._","",report$name) %in% filter_taxon)
  #pos.taxon <- which(report$name := filter_taxon)
  if (length(pos.taxons) == 0) {
    return(report)
  }

  row_seq <- seq_len(nrow(report))
  rows_to_delete <- rep(FALSE,nrow(report))

  taxon_depths <- report[pos.taxons,"depth"]
  if (isTRUE(rm_clade)) {
    taxonReads <- report[pos.taxons,"cladeReads"]
  } else {
    taxonReads <- report[pos.taxons,"taxonReads"]
    report[pos.taxons,"taxonReads"] <- 0
  }


  for (i in seq_along(pos.taxons)) {
    pos.taxon <- pos.taxons[i]
    if (pos.taxon == 1) {
      rows_to_delete[1] <- TRUE
      next
    }
    taxon_depth <- taxon_depths[i]
    taxonReads <- taxonReads[i]

    if (rm_clade) {
      tosum_below <-  row_seq >= pos.taxon & report$depth <= taxon_depth
      taxons_below <- cumsum(tosum_below) == 1
      rows_to_delete[taxons_below] <- TRUE
    }
    rows_to_update <- c(pos.taxon)

    taxons_above <- seq_len(nrow(report)) < pos.taxon & report$depth == taxon_depth

    any_stays <- FALSE
    prev_taxon_depth <- taxon_depth
    taxons_above <- c()
    for (i in seq(from=(pos.taxon-1),to=1)) {
      curr_taxon_depth <- report[i,"depth"]
      if (curr_taxon_depth < prev_taxon_depth) {
        if (!any_stays) {
          if (report[i,"cladeReads"] == taxonReads) {
            rows_to_delete[i] <- TRUE
            if (do_message)
              dmessage("Deleting ",report[i,"name"])
          } else {
            any_stays <- TRUE
          }
        }
        if (!rows_to_delete[i]) {
          rows_to_update <- c(rows_to_update, i)
          if (do_message)
            dmessage("Updating ",report[i,"name"])
        }
        prev_taxon_depth <- curr_taxon_depth
      } else {
        any_stays <- TRUE
      }
    }
    report[rows_to_update, "cladeReads"] <- report[rows_to_update, "cladeReads"] - taxonReads
  }

  #if (rm_clade)
  report[!rows_to_delete,]
  #else
  #  report
}




build_sankey_network <- function(my_report, taxRanks =  c("D","K","P","C","O","F","G","S"), maxn=10,
                                 zoom = F, title = NULL,
                                 ...) {
  stopifnot("taxRank" %in% colnames(my_report))
  if (!any(taxRanks %in% my_report$taxRank)) {
    warning("report does not contain any of the taxRanks - skipping it")
    return()
  }
  my_report <- subset(my_report, taxRank %in% taxRanks)
  my_report <- plyr::ddply(my_report, "taxRank", function(x) x[utils::tail(order(x$cladeReads,-x$depth), n=maxn), , drop = FALSE])

  my_report <- my_report[, c("name","taxLineage","taxonReads", "cladeReads","depth", "taxRank")]

  my_report <- my_report[!my_report$name %in% c('-_root'), ]
  #my_report$name <- sub("^-_root.", "", my_report$name)

  splits <- strsplit(my_report$taxLineage, "\\|")

  ## for the root nodes, we'll have to add an 'other' link to account for all cladeReads
  root_nodes <- sapply(splits[sapply(splits, length) ==2], function(x) x[2])

  sel <- sapply(splits, length) >= 3
  splits <- splits[sel]

  links <- data.frame(do.call(rbind,
                              lapply(splits, function(x) utils::tail(x[x %in% my_report$name], n=2))), stringsAsFactors = FALSE)
  colnames(links) <- c("source","target")
  links$value <- my_report[sel,"cladeReads"]

  my_taxRanks <- taxRanks[taxRanks %in% my_report$taxRank]
  taxRank_to_depth <- stats::setNames(seq_along(my_taxRanks)-1, my_taxRanks)


  nodes <- data.frame(name=my_report$name,
                      depth=taxRank_to_depth[my_report$taxRank],
                      value=my_report$cladeReads,
                      stringsAsFactors=FALSE)

  for (node_name in root_nodes) {
    diff_sum_vs_all <- my_report[my_report$name == node_name, "cladeReads"] - sum(links$value[links$source == node_name])
    if (diff_sum_vs_all > 0) {
      nname <- paste("other", sub("^._","",node_name))
      #nname <- node_name
      #links <- rbind(links, data.frame(source=node_name, target=nname, value=diff_sum_vs_all, stringsAsFactors = FALSE))
      #nodes <- rbind(nodes, nname)
    }
  }

  names_id = stats::setNames(seq_len(nrow(nodes)) - 1, nodes[,1])
  links$source <- names_id[links$source]
  links$target <- names_id[links$target]
  links <- links[links$source != links$target, ]

  nodes$name <- sub("^._","", nodes$name)
  links$source_name <- nodes$name[links$source + 1]

  if (!is.null(links))
    sankeyD3::sankeyNetwork(
      Links = links,
      Nodes = nodes,
      doubleclickTogglesChildren = TRUE,
      Source = "source",
      Target = "target",
      Value = "value",
      NodeID = "name",
      NodeGroup = "name",
      NodePosX = "depth",
      NodeValue = "value",
      #dragY = TRUE,
      xAxisDomain = my_taxRanks,
      numberFormat = "pavian",
      title = title,
      nodeWidth = 15,
      #linkGradient = TRUE,
      #nodeShadow = TRUE,
      nodeCornerRadius = 5,
      units = "cladeReads",
      fontSize = 12,
      iterations = maxn * 100,
      align = "none",
      highlightChildLinks = TRUE,
      orderByPath = TRUE,
      scaleNodeBreadthsByString = TRUE,
      zoom = zoom,
      ...
    )
}


main()
