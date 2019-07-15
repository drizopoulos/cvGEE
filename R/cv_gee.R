cv_gee <- function (object, rule = c("all", "quadratic", "logarithmic", "spherical"), 
                    max_count = 500, conf_int = FALSE, level = 0.95, 
                    M = 30, fold = 5, seed = 1L, 
                    return_data = FALSE) {
    if (!inherits(object, "geeglm")) {
        stop("function 'cv_gee()' currently only works for 'geeglm' objects.")
    }
    rule <- match.arg(rule)
    if (conf_int && M < 100) {
        warning("for a more stable estimation of the confidence interval for ", 
                "the scoring rule use a value of argument 'M' of 100 or greater.")
    }
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
            p_y <- dbinom(resp, size = N, prob = probs)
            log_p_y <- dbinom(resp, size = N, prob = probs, log = TRUE)
            quadrat_p <- numeric(n)
            for (i in seq_len(n)) {
                quadrat_p[i] <- sum(dbinom(max_count_seq[[i]], size = N[i], prob = probs[i])^2)
            }
            switch(rule, 
                   "all" = c(log_p_y, 2 * p_y + quadrat_p, p_y / sqrt(quadrat_p)),
                   "logarithmic" = log_p_y,
                   "quadratic" = 2 * p_y + quadrat_p, 
                   "spherical" = exp(log_p_y - 0.5 * log(quadrat_p)))
        } else if (object$family$family == "poisson") {
            counts <- object$family$linkinv(eta)
            log_p_y <- dpois(resp, lambda = counts, log = TRUE)
            max_count_seq <- lapply(max_count, seq, from = 0)
            quad_fun <- function (c1, c2) {
                sum(exp(2 * dpois(c1, lambda = c2, log = TRUE)))
            }
            quadrat_p <- mapply(quad_fun, c1 = max_count_seq, c2 = counts)
            switch(rule,
                   "all" = c(log_p_y, 2 * exp(log_p_y) + quadrat_p, 
                             exp(log_p_y - 0.5 * log(quadrat_p))),
                   "logarithmic" = log_p_y,
                   "quadratic" = 2 * exp(log_p_y) + quadrat_p, 
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
    #mfX_orig <- model.frame(terms(formula(object)), data = data)
    #if (!is.null(nas <- attr(mfX_orig, "na.action"))) {
    #  data <- data[-nas, ]
    #}
    K <- if (rule == "all") nrow(data) * 3 else nrow(data)
    out <- matrix(0.0, K, M)
    for (m in seq_len(M)) {
        splits <- split(seq_len(n), sample(rep(seq_len(fold), length.out = n)))
        for (i in seq_along(splits)) {
            id.i <- unq_id[splits[[i]]]
            train_data <- data[!id %in% id.i, ]
            test_data <- data[id %in% id.i, ]
            fit_i <- update(object, data = train_data)
            out[id %in% id.i, m] <- predict_geeglm(fit_i, test_data, data)
        }
    }
    if (return_data) {
        if (rule == "all") {
            data <- do.call('rbind', rep(list(data), 3))
            data$rule <- gl(3, nrow(data)/3, 
                            labels = c("logarithmic", "quadratic", "spherical"))
        }
        data$score <- rowMeans(out, na.rm = TRUE)
        if (conf_int) {
            data$low_score <- apply(out, 1, quantile, probs = (1 - level) / 2, na.rm = TRUE)
            data$upp_score <- apply(out, 1, quantile, probs = (1 + level) / 2, na.rm = TRUE)
        }
        data
    } else {
        rules <- rep(c("logarithmic", "quadratic", "spherical"), each = nrow(data))
        score <- split(rowMeans(out, na.rm = TRUE), rules)
        if (conf_int) {
            low_score <- split(apply(out, 1, quantile, probs = (1 - level) / 2, na.rm = TRUE), rules)
            upp_score <- split(apply(out, 1, quantile, probs = (1 - level) / 2), na.rm = TRUE, rules)
            score <- list(score = score, low_score = low_score, upp_score = upp_score)
        }
        score
    }
}