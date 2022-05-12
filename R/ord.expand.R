#' Expand data for model fitting
#'
#' This function expand each response value into indicator for model fitting
#'
#' @param formula an R formula
#' @param times a vector of times which occur within subject clusters
#' @param data a data frame including response data and covariates
#' @param subjects a character string specifying the subject variable
#' @param categories a numeric variable indicating the number of ordinal levels
#' @return an expanded data frame
ord.expand <- function(formula, times, data, subjects, categories){

  # remove missing data
  reorder <- match(c(sapply(attr(terms(formula), "variables"), deparse,
             width.cutoff = 500)[-1L], subjects), names(data))
  if (any(is.na(reorder))){stop("data: model frame and formula mismatch in model.matrix()")}
  data <- data[, sort(reorder)]
  split.data <- split(data, data[subjects])
  new.subjects  <- length(unique(data[,subjects]))

  # model formula
  orig.formula <- as.formula(formula)
  form.vars <- attr(terms.formula(as.formula(formula)), "variables")
  resp.var <- attr(terms.formula(as.formula(formula)), "response")
  term.labels <- attr(terms.formula(as.formula(formula)), "term.labels")
  if (resp.var == 0){stop("No response variable in formula")}
  resp.var <- resp.var + 1
  resp.label <- as.character(form.vars[[resp.var]])
  if (is.element("cuts", term.labels)){stop("data: term name cuts is not permitted")}
  formula <- as.formula(paste(resp.label, "~",
                              paste(c("cuts - 1", paste(term.labels, collapse = "+")), collapse = "+")))

  # new data frame
  scores <- resp.label
  namvars <- names(data)
  nvars <- length(names(data))
  categs1 <- categories - 1
  nlen <- categs1 * dim(data)[1]
  ndata <- data[1:nlen, ]
  rownames(ndata) <- 1:nlen
  for (j in 1:nvars){
  if (namvars[j] == scores){
   boscores <- matrix(0, ncol = categs1, nrow = dim(data)[1])
   for (i in 1:categs1){
    boscores[,i] <- as.numeric(data[, j] <= i)
    }
   ndata[,j] <- as.vector(t(boscores))
   cuts <- factor(rep(1:categs1, dim(data)[1]))
   }
  else {
   ndata[,j] <- rep(data[,j], each = categs1)
   }
   if (namvars[j] == subjects) {names(ndata)[j] <- "subjects"}
  }
  ndata <- data.frame(cuts = cuts,ndata)
  levels(ndata$cuts) <- paste(as.character(1:categories)[-categories],
                              as.character(1:categories)[-1], sep = "|")


  # output
  list(data = ndata, formula = formula)
}
