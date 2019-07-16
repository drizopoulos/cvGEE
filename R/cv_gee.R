cv_gee <- function (object, rule = c("all", "quadratic", "logarithmic", "spherical"), 
                    max_count = 500, K = 5L, M = 10L, seed = 1L, return_data = FALSE) {
    if (!inherits(object, "geeglm")) {
        stop("function 'cv_gee()' currently only works for 'geeglm' objects.")
    }
    rule <- match.arg(rule)
    predict_geeglm <- function (object, newdata, orig_data) {
        form <- formula(object)
        namesVars <- all.vars(form)
        orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
        Terms <- terms(form)
        mfX <- model.frame(Terms, data = orig_data)
        mfX_new <- model.frame(Terms, newdata, xlev = .getXlevels(Terms, mfX))
        X <- model.matrix(Terms, mfX_new)
        betas <- coef(object)
        eta <- c(X %*% betas)
        resp <- model.response(mfX_new)
        n <- length(resp)
        max_count <- rep(max_count, length.out = n)
        if (object$family$family == "gaussian") {
            warning("the family of 'object' is 'gaussian'; argument 'rule'", 
                    " is set to 'quadratic'.")
            rule <- "quadratic"
            (resp - object$family$linkinv(eta))^2
        } else if (object$family$family == "binomial") {
            if (NCOL(resp) == 2) {
                N <- max_count <- resp[, 1] + resp[, 2]
                resp <- resp[, 1]
            } else {
                N <- max_count <- rep(1, n)
            }
            max_count_seq <- lapply(max_count, seq, from = 0)
            probs <- object$family$linkinv(eta)

            log_p_y <- dbinom(resp, size = N, prob = probs, log = TRUE)
            quad_fun_binomial <- function (c1, c2, p) {
                sum(exp(2 * dbinom(x = c1, size = c2, prob = p, log = TRUE)))
            }
            quadrat_p <- mapply(quad_fun_binomial, c1 = max_count_seq, c2 = N, p = probs)
            switch(rule,
                   "all" = c(log_p_y, 2 * exp(log_p_y) - quadrat_p, 
                             exp(log_p_y - 0.5 * log(quadrat_p))),
                   "logarithmic" = log_p_y,
                   "quadratic" = 2 * exp(log_p_y) - quadrat_p, 
                   "spherical" = exp(log_p_y - 0.5 * log(quadrat_p)))
        } else if (object$family$family == "poisson") {
            counts <- object$family$linkinv(eta)
            log_p_y <- dpois(resp, lambda = counts, log = TRUE)
            max_count_seq <- lapply(max_count, seq, from = 0)
            quad_fun_poisson <- function (c1, c2) {
                sum(exp(2 * dpois(c1, lambda = c2, log = TRUE)))
            }
            quadrat_p <- mapply(quad_fun_poisson, c1 = max_count_seq, c2 = counts)
            switch(rule,
                   "all" = c(log_p_y, 2 * exp(log_p_y) - quadrat_p, 
                             exp(log_p_y - 0.5 * log(quadrat_p))),
                   "logarithmic" = log_p_y,
                   "quadratic" = 2 * exp(log_p_y) - quadrat_p, 
                   "spherical" = exp(log_p_y - 0.5 * log(quadrat_p)))
        }
    }
    ##################################################
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    data <- object$data
    id <- object$id
    unq_id <- unique(id)
    n <- length(unq_id)
    mfX_orig <- model.frame(terms(formula(object)), data = data)
    if (!is.null(nas <- attr(mfX_orig, "na.action"))) {
      data <- data[-nas, ]
    }
    Q <- if (rule == "all") nrow(data) * 3 else nrow(data)
    out <- matrix(0.0, Q, M)
    for (m in seq_len(M)) {
        splits <- split(seq_len(n), sample(rep(seq_len(K), length.out = n)))
        for (i in seq_along(splits)) {
            id.i <- unq_id[splits[[i]]]
            train_data <- data[!id %in% id.i, ]
            test_data <- data[id %in% id.i, ]
            fit_i <- update(object, data = train_data)
            out[id %in% id.i, m] <- predict_geeglm(fit_i, test_data, train_data)
        }
    }
    scores <- rowMeans(out, na.rm = TRUE)
    if (return_data) {
        if (rule == "all") {
            data <- do.call('rbind', rep(list(data), 3))
            labs <- c("logarithmic", "quadratic", "spherical")
            data[[".rule"]] <- gl(3, nrow(data)/3, labels = labs)
        }
        data[[".score"]] <- scores
        data
    } else {
        if (rule == "all") {
            rules <- rep(c("logarithmic", "quadratic", "spherical"), each = nrow(data))
            split(scores, rules)
        } else {
            scores
        }
    }
}