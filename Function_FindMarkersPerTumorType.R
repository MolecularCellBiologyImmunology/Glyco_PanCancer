FindMarkersPerTumorType <- function (object, ident.1, ident.2 = NULL, grouping.var, assay = "RNA", 
                                 slot = "data", meta.method = metap::minimump, verbose = TRUE, 
                                 ...) {
  if (!require("plyr")) {
    print("FindMarkersPerTumorType require the package 'plyr'. Proceeding to install...")
    install.packages('plyr')
    library('plyr')
  }
  object.var <- FetchData(object = object, vars = grouping.var)
  object <- SetIdent(object = object, cells = colnames(x = object), 
                     value = paste(Idents(object = object), object.var[, 
                                                                       1], sep = "_"))
  levels.split <- names(x = sort(x = table(object.var[, 1])))
  num.groups <- length(levels.split)
  cells <- list()
  for (i in 1:num.groups) {
    cells[[i]] <- rownames(x = object.var[object.var[, 1] == 
                                            levels.split[i], , drop = FALSE])
  }
  marker.test <- list()
  ident.2.save <- ident.2
  for (i in 1:num.groups) {
    level.use <- levels.split[i]
    ident.use.1 <- paste(ident.1, level.use, sep = "_")
    ident.use.1.exists <- ident.use.1 %in% Idents(object = object)
    if (!all(ident.use.1.exists)) {
      bad.ids <- ident.1[!ident.use.1.exists]
      warning("Identity: ", paste(bad.ids, collapse = ", "), 
              " not present in group ", level.use, ". Skipping ", 
              level.use, call. = FALSE, immediate. = TRUE)
      next
    }
    ident.2 <- ident.2.save
    cells.1 <- WhichCells(object = object, idents = ident.use.1)
    if (is.null(x = ident.2)) {
      cells.2 <- setdiff(x = cells[[i]], y = cells.1)
      ident.use.2 <- names(x = which(x = table(Idents(object = object)[cells.2]) > 
                                       0))
      ident.2 <- gsub(pattern = paste0("_", level.use), 
                      replacement = "", x = ident.use.2)
      if (length(x = ident.use.2) == 0) {
        stop(paste("Only one identity class present:", 
                   ident.1))
      }
    }
    else {
      ident.use.2 <- paste(ident.2, level.use, sep = "_")
    }
    if (verbose) {
      message("Testing group ", level.use, ": (", paste(ident.1, 
                                                        collapse = ", "), ") vs (", paste(ident.2, collapse = ", "), 
              ")")
    }
    ident.use.2.exists <- ident.use.2 %in% Idents(object = object)
    if (!all(ident.use.2.exists)) {
      bad.ids <- ident.2[!ident.use.2.exists]
      warning("Identity: ", paste(bad.ids, collapse = ", "), 
              " not present in group ", level.use, ". Skipping ", 
              level.use, call. = FALSE, immediate. = TRUE)
      next
    }
    marker.test[[i]] <- FindMarkers(object = object, assay = assay, 
                                    slot = slot, ident.1 = ident.use.1, ident.2 = ident.use.2,
                                    verbose = verbose, ...)
    names(x = marker.test)[i] <- levels.split[i]
  }
  
  marker.test <- Filter(f = Negate(f = is.null), x = marker.test)
  genes.conserved <- Reduce(f = intersect, x = lapply(X = marker.test, 
                                                      FUN = function(x) {
                                                        return(rownames(x = x))
                                                      }))
  markers.conserved <- list()
  for (i in 1:length(x = marker.test)) {
    markers.conserved[[i]] <- as.data.frame(t(marker.test[[i]]))
    rownames(x = markers.conserved[[i]]) <- paste(names(x = marker.test)[i], 
                                                  rownames(x = markers.conserved[[i]]), sep = "_")
  }
  
  markers.combined <- Reduce(rbind.fill, markers.conserved)
  markers.combined <- t(markers.combined)
  colnames(markers.combined)<-  unlist(lapply(markers.conserved, rownames))
  
  markers.combined <- as.data.frame(markers.combined)
  return(markers.combined)
}
