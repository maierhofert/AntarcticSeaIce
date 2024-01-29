
my_posterior_smooths.btl = function (object, fit, smooth, newdata = NULL, nsamples = NULL, 
          subset = NULL, ...) {
  # smooth <- (as_one_characte(smooth))
  smef <- brms:::tidy_smef(object, fit$data)
  smterms <- unique(smef$term)
  if (!smooth %in% smterms) {
    stop2("Term '", smooth, "' cannot be found. Available ", 
          "smooth terms are: ", collapse_comma(smterms))
  }
  sub_smef <- brms:::subset2(smef, term = smooth)
  covars <- all_vars(sub_smef$covars[[1]])
  byvars <- all_vars(sub_smef$byvars[[1]])
  req_vars <- c(covars, byvars)
  sdata <- standata(fit, newdata, re_formula = NA, internal = TRUE, 
                    check_response = FALSE, req_vars = req_vars)
  subset <- brms:::subset_samples(fit, subset, nsamples)
  samples <- as.matrix(fit, subset = subset)
  prep_args <- nlist(x = object, samples, sdata, data = fit$data)
  prep <- do_call(prepare_predictions, prep_args)
  i <- which(smterms %in% smooth)[1]
  J <- which(smef$termnum == i)
  scs <- unlist(attr(prep$sm$fe$Xs, "smcols")[J])
  prep$sm$fe$Xs <- prep$sm$fe$Xs[, scs, drop = FALSE]
  prep$sm$fe$bs <- prep$sm$fe$bs[, scs, drop = FALSE]
  prep$sm$re <- prep$sm$re[J]
  prep$family <- brmsfamily("gaussian")
  brms:::predictor(prep, i = NULL)
}
assignInNamespace("posterior_smooths.btl", my_posterior_smooths.btl, ns = "brms")
