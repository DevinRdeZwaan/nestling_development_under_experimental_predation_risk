################################################################################################
##################### Small change to piecewiseSEM internal code ###############################

### Description:
### Latest version of piecewiseSEM will generate an error if your path model has both fixed
### effects that are categorical and a random effect (i.e., mixed model structure).
### This is because package 'car' was updated recently such that the F statistic column
### from an ANOVA output is called simply 'F', rather than the previous 'Chisq | F statistic'.
### Run the entire code below which makes a small adjustment in the column header's name 
### and changes the summary.psem function within your R environment.

### All code is in its original format, created by Dr. J. Lefchek except for the small by-pass.
citation('piecewiseSEM')



################################################################################################
### Run from here to end of code
################################################################################################


#' Remove random effects from all.vars
#' 
#' @keywords internal
#' 
all.vars.merMod <- function(formula.) {
  
  if(!any(class(formula.) %in% c("formula", "formula.cerror"))) formula. <- formula(formula.)
  
  if(class(formula.) == "formula.cerror")
    
    gsub(" " , "", unlist(strsplit(formula., "~~"))) else {
      
      n <- rownames(attr(terms(formula.), "factors"))
      
      if(any(grepl("\\|", n)))
        
        all.vars(lme4::nobars(formula.)) else
          
          all.vars(formula.)
      
    }
  
}

#' Get vector of untransformed variables
#' 
#' @keywords internal
#' 
all.vars_notrans <- function(formula.) {
  
  if(!all(class(formula.) %in% c("formula", "formula.cerror"))) formula. <- formula(formula.)
  
  if(class(formula.) == "formula") {
    
    if(any(grepl("\\|", formula.))) formula. <- lme4::nobars(formula.)
    
    formula. <- all.vars_trans(formula.)
    
    if(any(grepl(":", formula.))) {
      
      idx <- which(grepl(":", formula.))
      
      for(i in idx) formula.[i] <- paste(sapply(strsplit(formula.[i], ":"), stripTransformations), collapse = ":")
      
      for(j in (1:length(formula.))[-idx]) formula.[j] <- stripTransformations(formula.[j])
      
    } else {
      
      formula. <- sapply(formula., stripTransformations)
      
    }
    
  } else formula. <- unlist(strsplit(formula., " ~~ "))
  
  return(formula.)
  
}

#' Get vector of transformed variables
#' 
#' @keywords internal
#' 
all.vars_trans <- function(formula.) {
  
  if(!all(class(formula.) %in% c("formula", "formula.cerror"))) formula. <- formula(formula.)
  
  if(class(formula.) == "formula") {
    
    if(formula.[[3]] == 1) deparse(formula.[[2]]) else {
      
      if(any(grepl("\\|", formula.))) formula. <- lme4::nobars(formula.)
      
      c(rownames(attr(terms(formula.), "factors"))[1], labels(terms(formula.)))
      
    }
    
  } else unlist(strsplit(formula., " ~~ "))
  
}


#' Captures output table
#' 
#' @keywords internal
#' 
captureTable <- function(g, row.names = FALSE) {
  
  g1 <- capture.output(print(g, row.names = row.names))
  
  if(all(g1 == "data frame with 0 columns and 0 rows")) 
    
    g1 <- "No independence claims present. Tests of directed separation not possible."
  
  g1 <- paste0(g1, "\n")
  
  return(g1)
  
}

#' Bind data.frames of differing dimensions
#'
#' From: https://stackoverflow.com/a/31678079
#' 
#' @param ... data.frames to be bound, separated by commas
#' 
#'  @keywords internal
#'   
cbind_fill <- function(...) {
  
  nm <- list(...) 
  
  dfdetect <- grepl("data.frame|matrix", unlist(lapply(nm, function(cl) paste(class(cl), collapse = " ") )))
  
  vec <- data.frame(nm[!dfdetect])
  
  n <- max(sapply(nm[dfdetect], nrow)) 
  
  vec <- data.frame(lapply(vec, function(x) rep(x, n)))
  
  if (nrow(vec) > 0) nm <- c(nm[dfdetect], list(vec))
  
  nm <- lapply(nm, as.data.frame)
  
  do.call(cbind, lapply(nm, function (df1) 
    
    rbind(df1, as.data.frame(matrix(NA, ncol = ncol(df1), nrow = n-nrow(df1), dimnames = list(NULL, names(df1))))) )) 
  
}

#' Transform variables based on model formula and store in new data frame
#' 
#' @keywords internal
#' 
dataTrans <- function(formula., data) {
  
  notrans <- all.vars.merMod(formula.)
  
  if(class(formula.) == "formula.cerror") notrans <- gsub(".*\\((.*)\\)", "\\1", notrans)
  
  trans <- all.vars_trans(formula.)
  
  trans <- unlist(strsplit(trans, "\\:"))
  
  trans <- trans[!duplicated(trans)]
  
  if(any(grepl("scale\\(.*\\)", trans))) {
    
    trans[which(grepl("scale(.*)", trans))] <- notrans[which(grepl("scale(.*)", trans))]
    
    warning("`scale` applied directly to variable. Use argument `standardize = TRUE` instead.", call. = FALSE)
    
  }
  
  if(any(!notrans %in% trans)) {
    
    for(k in 1:length(notrans)) {
      
      if(is.factor(data[, notrans[k]])) next else {
        
        data[, notrans[k]] <-
          
          sapply(data[, notrans[k]], function(x) eval(parse(text = gsub(notrans[k], x, trans[k]))))
        
      }
      
    }
    
  }
  
  colnames(data) <- notrans
  
  return(data)
  
}

#' Get ANOVA results
#'
#' @keywords internal
#'
getAnova <- function(model, test.statistic = "F", test.type = "III") {
  
  ct <- summary(model)$coefficients
  
  krp <- as.data.frame(car::Anova(model, test.statistic = test.statistic, type = test.type))
  
  ret <- cbind.data.frame(
    ct[, 1:2], 
    DF = krp[, 3], 
    Crit.Value = krp[, 1], 
    P = krp[, ncol(krp)]
  )
  
  names(ret)[ncol(ret)] <- "Pr(>|t|)"
  
  return(ret)
  
}

#' Get random effects from lme
#' 
#' @keywords internal
#' 
findbars.lme <- function(model) {
  
  rand <- model$call$random
  
  sapply(rand, function(i) {
    
    i = gsub(".*\\|(.*)", "\\1", as.character(i)[2])
    
    strsplit(gsub(" " , "", i), "\\/")[[1]]
    
  } )
  
}

#' Get data from model list
#' 
#' @keywords internal
#' 
GetData <- function(modelList) {
  
  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)
  
  modelList <- removeData(modelList, formulas = 1)
  
  data.list <- lapply(modelList, GetSingleData)
  
  data.list <- data.list[!sapply(data.list, is.null)]
  
  data.list <- unname(data.list)
  
  if(all(sapply(data.list, class) == "comparative.data"))
    
    data <- data.list[[1]] else 
      
      data <- do.call(cbind_fill, data.list)
  
  data <- data[, !duplicated(colnames(data), fromLast = TRUE)]
  
  # colnames(data) <- gsub(".*\\((.*)\\).*", "\\1", colnames(data))
  
  data <- as.data.frame(data)
  
  rownames(data) <- 1:nrow(data)
  
  return(data)
  
}

#' Get data from one model
#' 
#' @keywords internal
#' 
GetSingleData <- function(model) {
  
  dat <- data.frame()
  
  switch(class(model)[1],
         "phylolm" = {
           stop("Please provide `data =` argument to `psem`.", call. = FALSE)
         },
         
         "phyloglm" = {
           stop("Please provide `data =` argument to `psem`.", call. = FALSE)
         },
         
         "lm" ={
           dat <- eval(getCall(model)$data, environment(formula(model)))
         },
         
         "negbin" = {
           dat <- eval(getCall(model)$data, environment(formula(model)))
         },
         
         "sarlm" = {
           dat <- eval(getCall(model)$data, environment(formula(model)))
         },
         
         "glm" = {
           dat <- model$data
         },
         
         "glmmPQL" = {
           dat <- model$data
         },
         
         "pgls" = {
           dat <- model$data
         },
         
         "lmerMod" = {
           dat <- model@frame
         },
         
         "glmerMod" = {
           dat <- model@frame
         },
         
         "lmerModLmerTest" = {
           dat <- model@frame
         },
         
         "merModLmerTest" = {
           dat <- model@frame
         },
         
         "gls" = {
           dat <-  nlme::getData(model)
         },
         
         "lme" = {
           dat <-  nlme::getData(model)
         }
         
  )
  
  dat
  
}

#' Obtain (observation-level) random effects from a generalized linear mixed model
#' 
#' RE = "all" all random effects are reported
#' RE = "RE" just group effects are reported
#' RE = "OLRE" just observation-level effects are reported
#' 
#' @keywords internal
#' 
GetOLRE <- function(sigma, model, X, data, RE = c("all", "RE", "OLRE")) {
  
  if(class(model) %in% c("lmerMod", "glmerMod")) {
    
    if(is.null(X)) X <- model.matrix(model)
    
    rand <- sapply(lme4::findbars(formula(model)), function(x) as.character(x)[3])
    
    rand <- rand[!duplicated(rand)] 
    
  } 
  
  # else if(class(model) %in% c("lme", "glmmPQL")) { }
  
  idx <- sapply(sapply(strsplit(rand, "\\:"), function(x) gsub("\\(|\\)", "", x)), function(x) {
    
    length(unique(data[, x])) == nrow(data)
    
  } )
  
  sigma.names <- unlist(names(sigma)) # unlist(strsplit(names(sigma), "\\."))
  
  idx. <- sapply(sigma.names, function(x) !any(x %in% rand[idx]))
  
  if(RE == "RE") 
    
    out <- sapply(sigma[idx.], function(i) {
      
      Z <- as.matrix(X[, rownames(i), drop = FALSE])
      
      sum(rowSums(Z %*% i) * Z) / nrow(X)
      
    } ) else if(RE == "OLRE") {
      
      if(all(idx == FALSE)) out <- 0 else {
        
        out <- sapply(sigma[idx], function(i) {
          
          Z <- as.matrix(X[, rownames(i), drop = FALSE])
          
          sum(rowSums(Z %*% i) * Z) / nrow(X)
          
        } ) } } else if(RE == "all")
          
          out <- sapply(sigma, function(i) {
            
            Z <- as.matrix(X[, rownames(i), drop = FALSE])
            
            sum(rowSums(Z %*% i) * Z) / nrow(X)
            
          } )
  
  if(length(out) == 0) out <- 0
  
  return(out)
  
}

#' Get random effects variance-covariance from lme
#' 
#' @keywords internal
#' 
GetVarCov <- function(model) {
  
  vc <- try(getVarCov(model), silent = TRUE)
  
  if(any(class(vc) == "try-error")) {
    
    vc <- nlme::VarCorr(model)
    
    v <- suppressWarnings(as.numeric(vc[, 1]))
    
    names(v) <- gsub(" =", "", rownames(vc))
    
    vm <- as.list(na.omit(v[-length(v)]))
    
    vl <- lapply(1:length(vm), function(i) matrix(vm[[i]], dimnames = list(names(vm)[i], names(vm)[i])))
    
    names(vl) <- names(which(is.na(v)))
    
    vl
    
  } else list(vc)
  
}

#' Assess significance
#' 
#' @keywords internal
#' 
isSig <- function(p) {
  
  ifelse(p > 0.01 & p < 0.05, "*",
         ifelse(p > 0.001 & p <= 0.01, "**",
                ifelse(p <= 0.001, "***", "")))
  
}

#' Recompute P-values using Kenward-Rogers approximation
#' 
#' @keywords internal
#' 
# KRp <- function(model, vars, data, intercepts = FALSE) {
# 
#   # if(any(grepl("\\*", all.vars_notrans(formula(model)))) & !all(grepl("\\*", vars))) {
# 
#     f <- all.vars_trans(formula(model))
# 
#     model <- update(model, as.formula(paste(f[1], " ~ ", paste(f[-1], collapse = " + "), " + ", paste(onlyBars(formula(model)), collapse = " + "))))
# 
#   # }
# 
#   out <- data.frame()
#   
#   for(x in vars) { #sapply(vars, function(x) {
# 
#     reduceModel <- update(model, as.formula(paste(". ~ . -", x)))
# 
#     if(nobs(model) != nobs(reduceModel)) stop("Different sample sizes for `KRmodcomp`. Remove all NAs and re-run")
#     
#     kr <- try(pbkrtest::KRmodcomp(model, reduceModel), silent = TRUE)
# 
#     if(class(kr) == "try-error") 
#       
#       stop("Cannot obtain P-values from `lmerMod` using `pbkrtest::KRmodcopm`. Consider fitting using `nlme::lme`") else {
#         
#         d <- round(kr$stats$ddf, 2)
#   
#         p <- kr$stats$p.valueU
#   
#         out <- rbind(out, data.frame(d, p))
#         
#       }
# 
#   } # )
# 
#   if(intercepts == TRUE) {
# 
#     reduceModelI <- update(model, as.formula(paste("~ . - 1")), data = data)
# 
#     krI <- try(pbkrtest::KRmodcomp(model, reduceModelI), silent = TRUE)
#     
#     if(class(krI) == "try-error") 
#       
#       stop("Cannot obtain P-values from `lmerMod` using `pbkrtest::KRmodcomp`. Consider re-fitting using `nlme::lme`")else {
#         
#         dI <- krI$stats$ddf
#     
#         pI <- krI$stats$p.valueU
#     
#         out <- rbind(data.frame(d = dI, p = pI), out)
#         
#       }
#       
#   }
#   
#   return(out)
# 
# }

#' Get list of formula from a `psem` object
#' 
#' @keywords internal
#' 


listFormula <- function(modelList, formulas = 0) {
  
  modelList <- removeData(modelList, formulas)
  
  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)
  
  fList <- lapply(modelList, function(i) if(any(class(i) %in% c("formula.cerror"))) i else formula(i) )
  
  fList <- lapply(fList, lme4::nobars)
  
  return(fList)
  
}

#' Get number of observations from a model
#' 
#' @keywords internal
#' 
nObs <- function(object, ...) if(any(class(object) %in% c("phylolm", "phyloglm", "sarlm"))) length(fitted(object)) else nobs(object, ...)

#' Get random effects from merMod
#' 
#' @keywords internal
#' 
onlyBars <- function(formula., slopes = TRUE) {
  
  f <- lme4::findbars(formula.)
  
  if(slopes == TRUE) paste(sapply(f, function(x) paste0("(", deparse(x), ")")), collapse = " + ") else {
    
    f <- f[sapply(f, function(x) grepl("1\\||1 \\|", deparse(x)))]
    
    paste(sapply(f, function(x) paste0("(", deparse(x), ")")), collapse = " + ")
    
  }
  
}

#' Do not print attributes with custom functions
#' 
#' @keywords internal
#' 
#' @method print attr
#' 
print.attr <- function(x, ...) {
  
  attributes(x) <- NULL
  
  noquote(x)
  
}

#' Remove data from the model list
#' 
#' formulas = 0, keep everything
#' formulas = 1, remove all formulas including correlated errors
#' formulas = 2, remove only formula but keep correlated errors
#' formulas = 3, remove correlated errors but keep formula
#' 
#' @keywords internal
#' 
removeData <- function(modelList, formulas = 0) {
  
  remove <- c("character", "matrix", "data.frame", "SpatialPointsDataFrame", "comparative.data")
  
  if(formulas == 1) remove <- c(remove, "formula", "formula.cerror")
  
  if(formulas == 2) remove <- c(remove, "formula")
  
  if(formulas == 3) remove <- c(remove, "formula.cerror")
  
  modelList[!sapply(modelList, function(x) any(class(x) %in% remove))]
  
}

#' Strip transformations
#' 
#' @keywords internal
#' 
stripTransformations <- function(x) {
  
  x <- gsub(".*\\((.*)\\).*", "\\1", x)
  
  gsub(" ", "", gsub("(.*)\\+.*", "\\1", x))
  
}

#' Get Response Name as a Character
#' 
#' @keywords internal
#' 
get_response <- function(mod){
  mod <- removeData(mod)
  f <- lapply(mod, formula)
  r <- lapply(f, function(x) x[[2]])
  
  return(as.character(r))
}









################## coefficient code



coefs <- function(modelList, standardize = "scale", standardize.type = "latent.linear", test.statistic = "F", test.type = "II", intercepts = FALSE) {
  
  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)
  
  if(class(modelList) == "psem") data <- modelList$data else data <- GetData(modelList)
  
  if(class(data) %in% c("SpatialPointsDataFrame")) data <- data@data
  
  if(class(data) %in% c("comparative.data")) data <- data$data
  
  if(all(standardize != "none")) { 
    
    ret <- stdCoefs(modelList, data, standardize, standardize.type, test.statistic, test.type, intercepts) 
    
  } else {
    
    ret <- unstdCoefs(modelList, data, test.statistic, test.type, intercepts)
    
  }
  
  ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 4)
  
  ret[is.na(ret)] <- "-"
  
  names(ret)[length(ret)] <- ""
  
  if(any(ret$Estimate == "-")) warning("Categorical variables detected. Please refer to documentation for interpretation of Estimates!", call. = FALSE)
  
  return(ret)
  
}

#' Get raw (undstandardized) coefficients from model
#' 
#' @keywords internal
#' 
#' @export
#' 
unstdCoefs <- function(modelList, data = NULL, test.statistic = "F", test.type = "II", intercepts = FALSE) {
  
  if(!all(class(modelList) %in% c("list", "psem"))) modelList <- list(modelList)
  
  if(is.null(data)) data <- GetData(modelList)
  
  modelList <- removeData(modelList, formulas = 2)
  
  ret <- do.call(rbind, lapply(modelList, function(i) {
    
    if(all(class(i) %in% c("formula.cerror"))) {
      
      ret <- cerror(i, modelList, data)
      
      names(ret) <- c("Response", "Predictor", "Estimate", "Std.Error", "DF", "Crit.Value", "P.Value", "")
      
    } else {
      
      ret <- getCoefficients(i, data, test.statistic, test.type)
      
      if(intercepts == FALSE) ret <- ret[ret$Predictor != "(Intercept)", ]
      
    }
    
    return(ret)
    
  } ) )
  
  rownames(ret) <- NULL
  
  return(ret)
  
}

#' Get coefficients from linear regression
#' 
#' @keywords internal
#' 
#' @export
#'
getCoefficients <- function(model, data = NULL, test.statistic = "F", test.type = "II") {
  
  if(is.null(data)) data <- GetData(model)
  
  if(all(class(model) %in% c("lm", "glm", "negbin", "lmerMod", "glmerMod", "lmerModLmerTest", "pgls", "phylolm", "phyloglm"))) {
    
    ret <- as.data.frame(summary(model)$coefficients)
    
    if(all(class(model) %in% c("lm", "glm", "negbin"))) ret <- cbind(ret[, 1:2], DF = summary(model)$df[2], ret[, 3:4])
    
    if(all(class(model) %in% c("glmerMod", "pgls"))) ret <- cbind(ret[, 1:2], DF = length(summary(model)$residuals), ret[, 3:4])
    
    if(all(class(model) %in% c("phylolm", "phyloglm"))) ret <- cbind(ret[, 1:2], DF = model$n, ret[, c(3, 6)])
    
    if(all(class(model) %in% c("lmerMod"))) {
      
      # krp <- KRp(model, vars[-1], data, intercepts = TRUE)
      
      ret <- getAnova(model)
      
      # ret <- cbind(ret[, 1:2], DF = nobs(model), ret[, 3], P.Value = NA)
      
    }
    
  }
  
  if(all(class(model) %in% c("sarlm"))) {
    
    ret <- as.data.frame(summary(model)$Coef)
    
    ret <- cbind(ret[, 1:2], DF = NA, ret[, 3:4])
    
  }
  
  if(all(class(model) %in% c("gls", "lme", "glmmPQL"))) {
    
    ret <- as.data.frame(summary(model)$tTable)
    
    if(ncol(ret) == 4 & all(class(model) %in% c("gls")))
      
      ret <- cbind(ret[, 1:2], DF = length(residuals(model)), ret[, 3:4])
    
  }
  
  ret <- cbind(ret, isSig(ret[, 5]))
  
  ret <- data.frame(
    Response = all.vars_trans(listFormula(list(model))[[1]])[1],
    Predictor = rownames(ret),
    ret
  )
  
  names(ret) <- c("Response", "Predictor", "Estimate", "Std.Error", "DF", "Crit.Value", "P.Value", "")
  
  # if(sum(grepl("\\:", ret$Predictor)) > 0) warning("Interactions present. Interpret with care.")
  
  ret <- handleCategoricalCoefs2(ret, model, data, test.statistic, test.type)
  
  rownames(ret) <- NULL
  
  ret$Response <- as.character(ret$Response)
  
  ret[is.na(ret$P.Value), "Response"] <- ""
  
  return(ret)  
  
}

#' Handles putting categorical variables into coefficient tables
#' for easy use in path analysis
#' 
#' @keywords internal
#' 
handleCategoricalCoefs2 <- function(ret, model, data, test.statistic = "F", test.type = "II") {
  
  vars <- names(data)
  
  f <- listFormula(list(model))[[1]]
  
  mf <- model.frame(f, data)
  
  catVars <- names(mf)[sapply(mf, class) %in% c("factor", "character")]
  
  if(length(catVars) == 0) return(ret) else {
    
    for(i in catVars) {
      
      meanFacts <- suppressMessages(lapply(i, function(v) emmeans::emmeans(model, specs = v, data = data)))
      
      meanFacts <- lapply(meanFacts, function(m) {
        
        m <- as.data.frame(m)
        
        rownames(m) <- paste0(i, m[, 1])
        
        m[, 1] <- paste(names(m)[1], "=", as.character(m[, 1]))
        
        m$Crit.Value <- with(m, emmean/SE)
        
        m$P.Value <- with(m, 2 * pt(abs(Crit.Value), df, lower.tail = FALSE))
        
        colnames(m)[1] <- "Predictor"
        
        return(m)
        
      } )
      
      meanFacts <- do.call(rbind, meanFacts)
      
      meanFacts <- cbind(data.frame(Response = ret$Response[1], meanFacts))
      
      names(meanFacts)[names(meanFacts) %in% c("emmean", "SE", "df")] <- c("Estimate", "Std.Error", "DF")
      
      meanFacts <- meanFacts[, -which(colnames(meanFacts) %in% c("lower.CL", "upper.CL", "asymp.LCL", "asymp.UCL"))]
      
      meanFacts <- cbind(meanFacts, isSig(meanFacts[, 7]))
      
      names(meanFacts)[8] <- "sig"
      
      atab <- as.data.frame(car::Anova(model, test.statistic = test.statistic, type = test.type))
      
      atab <- atab[match(i, rownames(atab)), ]
      
      atab <- data.frame(Response = ret$Response[1], Predictor = i, Estimate = NA, Std.Error = NA, DF = atab$Df, 
                         Crit.Value = atab[, min(which(grepl("F", colnames(atab))))], P.Value = atab[, ncol(atab)], 
                         isSig(atab[, ncol(atab)]))
      
      colnames(atab) <- colnames(meanFacts)
      
      ret2 <- rbind(atab, meanFacts)
      
      colnames(ret2)[ncol(ret2)] <- ""
      
      retsp <- split(ret, grepl(i, rownames(ret)))
      
      retsp[[2]] <- ret2
      
      ret <- rbind(
        
        retsp[[1]],
        
        retsp[[2]]
        
      )
      
    } 
    
    return(ret)
    
  }
  
  # removed categorical interactions for now sorry JEKB
  #  
  # #what are the factors and their interactions?
  # factTypes <- rownames(modanova)[grep(catVars, rownames(modanova))]
  # 
  # #figure out which rows contain factors AND interactions
  # hasFactInt <- grep("\\:", factTypes)
  # 
  # intVars <- factTypes[hasFactInt]
  # 
  # #if there are interactions, use emmeans to get either
  # if(length(intVars)>0){
  #   intFacts <- lapply(intVars, deparseInt, model = model, catVars = catVars, vars = vars)
  #   intFacts <- do.call(rbind, intFacts)
  #   intFacts <- cbind(data.frame(Response = ret$Response[1]), intFacts)
  #   names(intFacts)[8] <- ""
  #   ret <- rbind(ret, intFacts)
  # }
  
}

#' Calculate standardized regression coefficients
#' 
#' @keywords internal
#' 
#' @export
#' 
stdCoefs <- function(modelList, data = NULL, standardize = "scale", standardize.type = "latent.linear", 
                     test.statistic = "F", test.type = "II",intercepts = FALSE) {
  
  if(!all(class(modelList) %in% c("list", "psem"))) modelList <- list(modelList)
  
  if(is.null(data) & class(modelList) == "psem") data <- modelList$data 
  
  if(is.null(data)) data <- GetData(modelList)
  
  modelList <- removeData(modelList, formulas = 2)
  
  ret <- do.call(rbind, lapply(modelList, function(i) {
    
    if(all(class(i) %in% c("formula.cerror"))) {
      
      ret <- cerror(i, modelList, data)
      
      cbind.data.frame(ret[, 1:7], Std.Estimate = ret[, 3], sig = ret[, 8])
      
    } else {
      
      ret <- unstdCoefs(i, data, test.statistic, test.type, intercepts)
      
      vars <- all.vars.merMod(i)
      
      newdata <- data[, vars]
      
      if(any(class(newdata) %in% c("SpatialPointsDataFrame"))) newdata <- newdata@data
      
      newdata <- dataTrans(formula(i), newdata)
      
      numVars <- attr(terms(i), "term.labels")
      
      if(any(grepl("\\:", numVars))) {
        
        v <- sapply(numVars, function(x) {
          
          if(any(grepl("\\:", x))) {
            
            x. <- strsplit(x, ":")[[1]]
            
            x. <- gsub("(.*) \\+.*", "\\1", gsub(".*\\((.*)\\)", "\\1", x.))
            
            ifelse(!all(sapply(x., function(y) class(newdata[, y]) == "factor")), TRUE, FALSE)
            
          } else TRUE } )
        
        numVars <- numVars[v]
        
      }
      
      B <- ret[ret$Predictor %in% c("(Intercept)", numVars), "Estimate"]
      
      names(B) <- numVars
      
      sd.x <- GetSDx(i, modelList, newdata, standardize)
      
      sd.y <- GetSDy(i, newdata, standardize, standardize.type)
      
      if(intercepts == FALSE)
        
        Std.Estimate <- B * (sd.x / sd.y) else
          
          Std.Estimate <- c(0, B[-1] * (sd.x / sd.y))
      
      ret[which(ret$Predictor %in% c("(Intercept)", numVars)), "Std.Estimate"] <- Std.Estimate
      
      ret <- cbind(ret[, 1:7], ret[, 9, drop = FALSE], sig = ret[, 8])
      
      return(ret)
      
    }
    
  } ) )  
  
  rownames(ret) <- NULL
  
  return(ret)
  
}

#' Get standard deviation of predictor variables
#' 
#' @keywords internal
#' 
GetSDx <- function(model, modelList, data, standardize = "scale") {
  
  vars <- all.vars.merMod(model)
  
  numVars <- vars[which(sapply(data[, vars], class) != "factor")]
  
  if(all(standardize == "scale"))
    
    sd.x <- sapply(numVars[!grepl(":", numVars)][-1], function(x) sd(data[, x], na.rm = TRUE)) else
      
      if(all(standardize == "range"))
        
        sd.x <- sapply(numVars[!grepl(":", numVars)][-1], function(x) diff(range(data[, x], na.rm = TRUE))) else
          
          if(is.list(standardize)) {
            
            vars <- unlist(sapply(modelList, all.vars_notrans))
            
            vars <- vars[!grepl(":", vars)]
            
            if(!all(names(standardize) %in% vars))
              
              stop("Names in standardize list must match those in the model formula!")
            
            sd.x <- sapply(numVars[!grepl(":", numVars)][-1], function(x) {
              
              nm <- which(names(standardize) == x)
              
              if(sum(nm) == 0) {
                
                warning(paste0("Relevant range not specified for variable '", x, "'. Using observed range instead"), call. = FALSE)
                
                diff(range(data[, x], na.rm = TRUE)) 
                
              } else  diff(range(standardize[[nm]]))
              
            } )
            
          } else stop("`standardize` must be either 'scale' or 'range' (or a list of ranges).", call. = FALSE)
  
  if(any(grepl(":", all.vars_notrans(model)))) sd.x <- c(sd.x, scaleInt(model, data, standardize))
  
  if(length(sd.x) == 0) sd.x <- NA
  
  return(sd.x)
  
}

#' Properly scale standard deviations depending on the error distribution
#' 
#' @keywords internal
#' 
GetSDy <- function(model, data, standardize = "scale", standardize.type = "latent.linear") {
  
  vars <- all.vars.merMod(model)
  
  y <- vars[1]
  
  family. <- try(family(model), silent = TRUE)
  
  if(class(family.) == "try-error") family. <- try(model$family, silent = TRUE)
  
  if(class(family.) == "try-error" | is.null(family.) & all(class(model) %in% c("sarlm", "gls", "lme")))
    
    family. <- list(family = "gaussian", link = "identity")
  
  if(class(family.) == "try-error" | is.null(family.) | any(class(model) %in% c("glmerMod", "glmmPQL")))
    
    sd.y <- NA else {
      
      if(family.$family == "gaussian") {
        
        if(all(standardize == "scale")) sd.y <- sd(data[, y], na.rm = TRUE) else
          
          if(all(standardize == "range")) sd.y <- diff(range(data[, y], na.rm = TRUE)) else
            
            if(is.list(standardize)) {
              
              nm <- which(names(standardize) == y)
              
              if(sum(nm) == 0) {
                
                warning(paste0("Relevant range not specified for variable '", y, "'. Using observed range instead"), call. = FALSE)
                
                sd.y <- diff(range(data[, y], na.rm = TRUE)) 
                
              } else sd.y <- diff(range(standardize[[nm]]))
              
            }
          
      } else if(family.$family == "binomial")
        
        sd.y <- scaleGLM(model, standardize, standardize.type) else
          
          sd.y <- NA
        
    }
  
  return(sd.y)
  
}

#' Compute standard deviation or relevant range of response for GLMs
#' 
#' @keywords internal
#' 
scaleGLM <- function(model, standardize = "scale", standardize.type = "latent.linear") {
  
  preds <- predict(model, type = "link")
  
  if(standardize.type == "Menard.OE") {
    
    y <- all.vars_notrans(model)[1]
    
    data <- GetSingleData(model)
    
    R <- cor(data[, y], predict(model, type = "response"))
    
    sd.y <- sqrt(var(preds)) / R
    
  }
  
  if(standardize.type == "latent.linear") {
    
    link. <- family(model)$link
    
    if(link. == "logit") sigmaE <- pi^2/3 else
      
      if(link. == "probit") sigmaE <- 1
      
      sd.y <- sqrt(var(preds) + sigmaE)
      
  }
  
  if(all(standardize == "range") | is.list(standardize)) sd.y <- sd.y * 6
  
  return(sd.y)
  
}

#' Calculate standard deviation or relevant range for interaction terms
#' 
#' @keywords internal
#' 
scaleInt <- function(model, newdata, standardize) {
  
  v <- attr(terms(model), "term.labels")
  
  int <- v[grepl(":", v)]
  
  sapply(int, function(x) {
    
    x <- strsplit(x, ":")[[1]]
    
    x <- gsub("(.*) \\+.*", "\\1", gsub(".*\\((.*)\\)", "\\1", x))
    
    p <- apply(newdata[, x], 1, prod, na.rm = TRUE)
    
    if(standardize == "scale") sd(p) else if(standardize == "range")
      
      diff(range(p, na.rm = TRUE)) else if(is.list(standardize))
        
        stop("Relevant range standardization not applicable to models with interactions!")
    
  } )
  
}

################################# summary code

summary.psem <- function(object, ...,
                         basis.set = NULL, direction = NULL, conserve = FALSE, conditioning = FALSE,
                         add.claims = NULL,
                         standardize = "scale", standardize.type = "latent.linear", 
                         test.statistic = "F", test.type = "II",
                         intercepts = FALSE,
                         .progressBar = TRUE) {
  
  name <- deparse(match.call()$object)
  
  call <- paste(listFormula(object), collapse = "\n  ")
  
  dTable <- dSep(object, basis.set, direction, conserve, conditioning, .progressBar)
  
  Cstat <- fisherC(dTable, add.claims, direction, conserve, conditioning, .progressBar)
  
  IC <- infCrit(object, Cstat, add.claims, direction, conserve, conditioning, .progressBar)
  
  coefficients <- coefs(object, standardize, standardize.type, test.statistic, test.type, intercepts)
  
  R2 <- rsquared(object)
  
  R2[, which(sapply(R2, is.numeric))] <- round(R2[, which(sapply(R2, is.numeric))], 2)
  
  if(length(dTable) > 0)
    
    dTable[, which(sapply(dTable, is.numeric))] <- round(dTable[, which(sapply(dTable, is.numeric))], 4)
  
  l <- list(name = name, call = call, dTable = dTable, Cstat = Cstat, IC = IC, coefficients = coefficients, R2 = R2)
  
  class(l) <- "summary.psem"
  
  l
  
} 



### This line is critical. It makes the change to the internal function in the R workspace.
assignInNamespace("handleCategoricalCoefs", handleCategoricalCoefs2, ns="piecewiseSEM")
