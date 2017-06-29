### Pedigree.R from https://bioconductor.org/packages/release/bioc/html/GeneticsPed.html
###------------------------------------------------------------------------
### What: Pedigree class and methods
### $Id: pedigree.R 1165 2007-04-03 14:00:07Z ggorjan $
### Time-stamp: <2007-09-20 13:18:57 ggorjan>
###------------------------------------------------------------------------

### {{{ Pedigree

## setMethod(f="PedigreeS4", signature="data.frame",
##           def=function(x, subject="id", ascendant=c("father", "mother"),
##                        ascendantSex=c(1, 2), ascendantLevel=c(1, 1),
##                        unknown=NULL, generation=NULL, generationOrder="increasing",
##                        family=NULL, sex=NULL, sexName=c("male", "female"),
##                        dtBirth=NULL, other=NULL, gene=NULL, ...) {
##            x
##          })

## setMethod(f="PedigreeS4", signature="matrix",
##          def=function(x, ...) {
##            ## TODO
##            x
##          })

### }}}
### {{{ S3

###------------------------------------------------------------------------

Pedigree <- function(x, subject="id", ascendant=c("father", "mother"),
                     ascendantSex=c(1, 2), ascendantLevel=c(1, 1),
                     unknown=NA, sex=NA, dtBirth=NA, generation=NA,
                     family=NA, generationOrder="increasing", check=TRUE,
                     sort=FALSE, extend=FALSE, drop=TRUE, codes=FALSE) {
  
  ## --- Attributes ---
  
  class(x) <- c("Pedigree", "data.frame")
  
  ## FIXME: sex, dtBirth, generation, family, default to NULL and then
  ## if(is.null) NA
  
  attr(x, ".subject") <- subject
  attr(x, ".ascendant") <- ascendant
  attr(x, ".family") <- family
  attr(x, ".ascendantSex") <- ascendantSex
  attr(x, ".ascendantLevel") <- as.integer(ascendantLevel)
  attr(x, ".sex") <- sex
  attr(x, ".dtBirth") <- dtBirth
  attr(x, ".generation") <- generation
  attr(x, ".generationOrder") <- generationOrder
  attr(x, ".colClass") <- class(x[[subject]])[1]
  attr(x, ".checked") <- FALSE
  attr(x, ".sorted") <- FALSE
  attr(x, ".extended") <- FALSE
  attr(x, ".coded") <- FALSE
  attr(x, ".unknown") <- list(.id=NA, .family=NA, .sex=NA, .dtBirth=NA,
                              .generation=NA)
  
  ## --- Actions ---
  
  for(i in c(subject, ascendant, sex, dtBirth, generation)) {
    if(!is.na(i)) {
      if(is.null(x[[i]])) stop(sprintf("column %s does not exist", sQuote(i)))
    }
  }
  
  if(!(generationOrder %in% c("increasing", "decreasing"))) {
    stop(sprintf("'generationOrder' must %s or %s",
                 dQuote("increasing"), dQuote("decreasing")))
  }
  
  if(length(ascendant) != length(ascendantSex) |
     length(ascendant) != length(ascendantLevel) |
     length(ascendantSex) != length(ascendantLevel)) {
    stop("'ascendant', 'ascendantSex' and 'ascendantLevel' must have the same length")
  }
  
  #idClass(x)
  ## all id must have the same class: not needed when check
  ## action will be working FIXME
  
  if(any(!is.na(unlist(unknown)))) {
    ## FIXME - remove this when you will handle unknown stuff
    ##    x <- unknownToNA.Pedigree(x=x, unknown=unknown)
    for(i in ascendant) x[x[, i] %in% unknown, i] <- NA
  }
  
  if(!is.na(sex)) {
    if(!all(attr(x, ".ascendantSex") %in% x[[attr(x, ".sex")]]))
      stop("values of 'ascendantSex' must accord with values of 'sex' column")
  }
  
  ##  if(check) # Check
  ##    check(x)
  
  if(sort) # Sort
    x <- sort.Pedigree(x)
  
  if(extend) # Extend
    x <- extend(x)
  
  if(attr(x, ".colClass") == "factor" & drop) # drop unused levels for factors
    x <- dropLevels(x)
  
  if(codes) # Code Pedigree
    stop("not available")
  ## x <- as.integer.Pedigree(x)
  
  ## --- Return ---
  x
}

## --- print ---

if(FALSE) {
  print.Pedigree <- function(x, ...)
    summary(x)
}

### }}}
### {{{ is., as.

## --- is. ---

is.Pedigree <- function(x) inherits(x, "Pedigree")

## --- as. ---

as.Pedigree <- function (x, ...) UseMethod("as.Pedigree")

as.Pedigree.Pedigree <- function(x, ...) x

as.Pedigree.data.frame <- function (x, ...) Pedigree(x, ...)

as.Pedigree.matrix <- function (x, ...)
{
  if(is.null(colnames(x)))
    stop("matrix 'x' needs column names to proceed")
  Pedigree(as.data.frame(x), ...)
}

if(FALSE) {
  
  as.integer.Pedigree <- function(x, sort=FALSE, mapCha=FALSE, mapInt=FALSE, ...)
  {
    col <- NULL
    ## FIXME: what to do with other columns i.e. factors, sex, anything else?
    ## FIXME: get 1:n with sort?
    if(sort) x <- sort(x)
    subject <- attr(x, ".subject")
    if(is.null(col)) col <- c(subject, attr(x, ".ascendant"))
    ## FIXME: combine=TRUE will be OK only for id columns!
    mapLCha <- mapLevels(x[, col], codes=FALSE, combine=TRUE)
    ## Applying the same map so that integer codes will be consistent
    mapLevels(x[, col]) <- mapLCha
    mapLInt <- mapLevels(x[[subject]])
    x[, col] <- lapply(x[, col], as.integer)
    attr(x, ".colClass") <- "integer"
    if(mapCha | mapInt) {
      ret <- list()
      ret$pedigree <- x
      if(mapCha) ret$mapCha <- mapLCha
      if(mapInt) ret$mapInt <- mapLInt
      ret
    } else {
      x
    }
  }
  
  as.factor.Pedigree <- function(x, map, ...)
  {
    col <- NULL
    ## FIXME: what to do with other columns i.e. factors, sex, anything else?
    if(is.null(col)) col <- c(attr(x, ".subject"), attr(x, ".ascendant"))
    if(!missing(map)) {
      mapLevels(x[, col]) <- map
    } else {
      x[, col] <- lapply(x[, col], as.factor)
    }
    attr(x, ".colClass") <- "factor"
    x
  }
  
  as.character.Pedigree <- function(x, ...)
  {
    col <- NULL
    ## FIXME: what to do with other columns i.e. factors, sex, anything else?
    if(is.null(col)) col <- c(attr(x, ".subject"), attr(x, ".ascendant"))
    x[, col] <- lapply(x[, col], as.character)
    attr(x, ".colClass") <- "character"
    x
  }
  
}

### }}}
### {{{ [, [[

if(FALSE) {
  
  ## not easy as it is jumping from .data.frame to .Pedigree. Probably better
  ## to go with lists in general.
  
  "[.Pedigree" <- function(x, i, j, drop=TRUE)
  {
    ## --- Setup ---
    
    attrX <- attributes(x)                        ## save attributes
    classX <- class(x)
    class(x) <- classX[!(classX %in% "Pedigree")] ## remove Pedigree class
    
    nCol <- ncol(x)
    
    args <- vector(mode="list", length=3)
    names(args) <- c("i", "j", "drop")
    if(missing(drop)) args$drop <- NA else args$drop <- drop
    
    iMiss <- missing(i)
    jMiss <- missing(j)
    if(jMiss) args$j <- 1:nCol    else args$j <- j
    if(iMiss) args$i <- 1:nrow(x) else args$i <- i
    if(jMiss & iMiss) return(x[, , drop=drop])
    
    cat("Subsetting\n")
    
    ## --- Test if subseting changes extended status ---
    
    if(!iMiss) {
      n <- nrow(x[i, ])
      nAll <- nrow(extend(x[i, ]))
      if(n < nAll) attrX$.extended <- FALSE
      
      ## --- i can change sort ---
      
      rowName <- row.names(x)
      rowName <- rowName[rowName %in% i]
      if(!identical(order(i), order(rowName))) attrX$.sorted <- FALSE
    }
    
    ## --- j can remove some data ---
    
    ## which columns are selected
    colName <- names(x)
    if(is.character(args$j)) {
      colTest <- colName %in% args$j
    } else {
      colTest <- 1:nCol %in% args$j
    }
    
    if(!jMiss & !all(colTest)) {
      ## Which columns will be removed? - set their attributes to NA
      colName <- colName[!colTest]
      attrX[attrX %in% colName] <- NA
      ## Additional loss of information
      ascendant <- attr(x, ".ascendant")
      ascTest <- ascendant %in% colName
      attrX$.ascendant <- attrX$.ascendant[!ascTest]
      attrX$.ascendantLevel <- attrX$.ascendantLevel[!ascTest]
      attrX$.ascendantSex <- attrX$.ascendantSex[!ascTest]
      if(length(attrX$.ascendant) == 0) ## all ascendants were removed
        attrX$.ascendant <- attrX$.ascendantLevel <- attrX$.ascendantSex <- NA
      generation <- attr(x, ".generation")
      generTest <- generation %in% colName
      attrX$.generationOrder <- NA
    }
    
    ## --- Attributes and [.data.frame method ---
    
    args[is.na(args)] <- NULL ## remove unspecified args
    x <- do.call(what="[.data.frame", args=c(list(x=x), args))
    if(is.data.frame(x)) { ## set attributes back
      if(!iMiss) attrX$row.names <- row.names(x)
      if(!jMiss) attrX$names <- names(x)
      attributes(x) <- attrX
    }
    x
  }
  
  ## "[<-.Pedigree" <- function(x, i, j, value)
  ## {
  ##
  ## }
  ##
  ## "[[", "[[<-"m "$", "$<-"
}

### }}}
### {{{ Dear Emacs
## Local variables:
## folded-file: t
## End:
### }}}

###------------------------------------------------------------------------
### Pedigree.R ends here
