project_number <- function(value_2011, value_2023, years = 2025:2050) {
  g <- (value_2023 / value_2011)^(1 / 12) - 1
  g <- pmin(0.02, pmax(-0.95, g))
  out <- value_2023 * (1 + g)^(years - 2023)
  out[!is.finite(out)] <- 0
  out
}

morbidity_proxies <- function(d) {
  med <- function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
  D <- ifelse(d$deaths > 0 & d$yll > 0, d$yll / d$deaths, NA_real_)
  for (group in list(interaction(d$country, d$year, d$sex, drop = TRUE),
                     interaction(d$year, d$sex, drop = TRUE), rep(1, nrow(d)))) {
    fill <- ave(D, group, FUN = med)
    D[!is.finite(D)] <- fill[!is.finite(D)]
  }
  d$D <- pmax(0.1, D)
  d$m <- d$yld / (d$population * d$D)
  d$p <- pmin(0.999999, pmax(0, 1 - d$incidence / d$prevalence))
  d$m[!is.finite(d$m)] <- 0
  d$p[!is.finite(d$p)] <- 0
  d
}

counterfactual_cell <- function(population, lfpr, sigma_history,
                                m_history, p_history) {
  stopifnot(all(is.finite(sigma_history)), all(sigma_history >= 0),
            all(sigma_history < 1), length(m_history) == length(p_history))
  m_history[!is.finite(m_history)] <- 0
  p_history[!is.finite(p_history)] <- 0
  persistence <- pmin(0.999999, pmax(0, p_history))
  retained <- 1 - m_history * persistence^(seq_along(m_history) - 1)
  stopifnot(all(retained > 0), all(retained <= 1))
  c(population_cf = population / prod(1 - sigma_history),
    lfpr_cf = min(0.999999, pmin(0.999, pmax(0.001, lfpr)) / prod(retained)))
}

effective_human_capital <- function(cells) {
  x <- cells$age_midpoint - cells$schooling - 5
  h <- exp(0.091 * cells$schooling + 0.1301 * x - 0.0023 * x^2)
  aggregate(cbind(H = h * cells$lfpr * cells$population,
                  H_cf = h * cells$lfpr_cf * cells$population_cf),
            list(year = cells$year), sum)
}

forecast_cost <- function(year, cost, log_scale = FALSE, level = 95) {
  stopifnot(requireNamespace("forecast", quietly = TRUE), length(cost) >= 10,
            length(year) == length(cost), all(diff(year) == 1), max(year) < 2024,
            all(is.finite(cost)), all(cost >= 0), !log_scale || all(cost > 0))
  predict_model <- function(y, h, model) {
    z <- ts(if (log_scale) log(y) else y, frequency = 1)
    fit <- switch(model,
      holt = forecast::holt(z, h = h, damped = TRUE, level = level),
      ets = forecast::forecast(forecast::ets(z), h = h, level = level),
      arima = forecast::forecast(forecast::auto.arima(z), h = h, level = level))
    out <- cbind(estimate = as.numeric(fit$mean),
                 lower = as.numeric(fit$lower[, 1]), upper = as.numeric(fit$upper[, 1]))
    if (log_scale) out <- exp(out)
    pmax(out, 0)
  }
  candidates <- c("holt", "ets", "arima")
  train <- head(cost, -5)
  test <- tail(cost, 5)
  errors <- t(vapply(candidates, function(model) {
    pred <- tryCatch(predict_model(train, 5, model)[, "estimate"],
                     error = function(e) rep(NA_real_, 5))
    if (any(!is.finite(pred))) return(c(smape = Inf, mae = Inf))
    denominator <- abs(test) + abs(pred)
    c(smape = 100 * mean(ifelse(denominator == 0, 0,
                                2 * abs(test - pred) / denominator)),
      mae = mean(abs(test - pred)))
  }, c(smape = 0, mae = 0)))
  stopifnot(any(is.finite(errors[, "smape"])))
  best <- candidates[order(errors[, "smape"], errors[, "mae"])[1]]
  horizon <- 2024 - max(year)
  list(model = best, validation = errors,
       baseline_2024 = predict_model(cost, horizon, best)[horizon, ])
}

simulate_gdp <- function(year, gdp, H, H_cf, K0, s, beta, cost_2024, delta = 0.05) {
  n <- length(year)
  stopifnot(n >= 2, all(diff(year) == 1), length(gdp) == n,
            length(H) == n, length(H_cf) == n,
            all(is.finite(c(gdp, H, H_cf, K0, s, beta, cost_2024))),
            all(gdp > 0), all(H > 0), all(H_cf > 0), length(K0) == 1, K0 > 0,
            length(s) == 1, length(beta) == 1, beta > 0, beta < 1,
            length(cost_2024) == 1, cost_2024 >= 0, year[1] == 2025,
            length(delta) == 1, is.finite(delta), delta >= 0, delta <= 1)
  K <- K_cf <- Y_cf <- A <- numeric(n)
  K[1] <- K_cf[1] <- K0
  TC <- cost_2024 * gdp / gdp[1]
  for (t in seq_len(n)) {
    if (K[t] <= 0 || K_cf[t] <= 0) stop("Non-positive capital; check inputs.")
    A[t] <- gdp[t] / (K[t]^(1 - beta) * H[t]^beta)
    Y_cf[t] <- A[t] * K_cf[t]^(1 - beta) * H_cf[t]^beta
    if (t < n) {
      K[t + 1] <- (1 - delta) * K[t] + s * gdp[t]
      K_cf[t + 1] <- (1 - delta) * K_cf[t] + s * Y_cf[t] + s * TC[t]
    }
  }
  data.frame(year, gdp, gdp_cf = Y_cf, loss = Y_cf - gdp, K, K_cf, A, TC)
}

summarise_loss <- function(annual, population_2025, r = 0.03) {
  stopifnot(length(r) == 1, is.finite(r), r > -1, population_2025 > 0)
  weight <- (1 + r)^-(annual$year - 2025)
  loss <- sum(annual$loss * weight)
  discounted_gdp <- sum(annual$gdp * weight)
  c(loss = loss, discounted_gdp = discounted_gdp,
    share_percent = 100 * loss / discounted_gdp, per_capita = loss / population_2025)
}

impute_countries <- function(d) {
  med <- function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
  observed <- is.finite(d$estimate)
  missing <- !observed
  stopifnot(all(d$lower[observed] <= d$estimate[observed]),
            all(d$estimate[observed] <= d$upper[observed]))
  pc <- d$discounted_gdp / d$population_2025
  for (i in which(!is.finite(d$discounted_gdp))) {
    d$discounted_gdp[i] <- med(pc[d$region == d$region[i]]) * d$population_2025[i]
  }
  stopifnot(all(is.finite(d$discounted_gdp)), all(d$discounted_gdp > 0))
  d$share <- d$estimate / d$discounted_gdp
  train <- d[observed & complete.cases(d[c("share", "incidence", "deaths")]), ]
  fit <- if (nrow(train) > 3) lm(share ~ incidence + deaths, data = train) else NULL
  predicted <- rep(NA_real_, nrow(d))
  if (!is.null(fit) && fit$rank == 3) predicted <- as.numeric(predict(fit, newdata = d))
  for (i in which(missing & !is.finite(predicted))) {
    predicted[i] <- med(d$share[observed & d$region == d$region[i]])
    if (!is.finite(predicted[i])) predicted[i] <- med(d$share[observed])
  }
  spread <- (d$upper - d$lower) / (2 * qnorm(0.975) * d$discounted_gdp)
  sd_mean <- mean(spread[observed & is.finite(spread)], na.rm = TRUE)
  stopifnot(all(is.finite(predicted[missing])), !any(missing) || is.finite(sd_mean))
  centre <- pmax(0, predicted[missing])
  d$estimate[missing] <- centre * d$discounted_gdp[missing]
  d$lower[missing] <- pmax(0, predicted[missing] - qnorm(0.975) * sd_mean) *
    d$discounted_gdp[missing]
  d$upper[missing] <- pmax(centre, predicted[missing] + qnorm(0.975) * sd_mean) *
    d$discounted_gdp[missing]
  d$share <- d$estimate / d$discounted_gdp
  d$imputed <- missing
  d
}

ramsey_rate <- function(gdp_pc_2025, gdp_pc_2050) {
  stopifnot(all(is.finite(c(gdp_pc_2025, gdp_pc_2050))),
            all(gdp_pc_2025 > 0), all(gdp_pc_2050 > 0))
  pmax(0, 2 * ((gdp_pc_2050 / gdp_pc_2025)^(1 / 25) - 1))
}

ramsey_loss <- function(loss, r, rates = c(0, 0.02, 0.03, 0.05)) {
  stopifnot(length(loss) == length(rates), length(rates) >= 2,
            all(diff(rates) > 0), is.finite(r), r >= 0, all(is.finite(loss)))
  if (all(loss == 0)) return(0)
  if (any(loss <= 0)) stop("Log interpolation requires positive losses; use annual results.")
  j <- max(1, min(length(rates) - 1, findInterval(r, rates)))
  exp(log(loss[j]) + (r - rates[j]) / (rates[j + 1] - rates[j]) *
        (log(loss[j + 1]) - log(loss[j])))
}

result_psa <- function(results, rho_country = 0.5, rho_disease = c(0, 0.5, 1),
                       n_country = 1000, n_global = 100000, unit = 1, seed = 2025) {
  stopifnot(!anyDuplicated(results[c("country", "disease")]),
            all(is.finite(as.matrix(results[c("estimate", "lower", "upper")]))),
            all(results$lower >= 0), all(results$lower <= results$estimate),
            all(results$estimate <= results$upper), unit > 0,
            rho_country >= 0, rho_country <= 1, all(rho_disease >= 0),
            all(rho_disease <= 1), n_country >= 2, n_global >= 2)
  set.seed(seed)
  disease_draws <- lapply(split(results, results$disease), function(d) {
    n <- nrow(d)
    z <- sqrt(rho_country) * matrix(rnorm(n_country), n_country, n) +
      sqrt(1 - rho_country) * matrix(rnorm(n_country * n), n_country, n)
    mu <- log1p(d$estimate / unit)
    sd_lo <- (mu - log1p(d$lower / unit)) / qnorm(0.975)
    sd_hi <- (log1p(d$upper / unit) - mu) / qnorm(0.975)
    deviation <- ifelse(z < 0, sweep(z, 2, sd_lo, `*`), sweep(z, 2, sd_hi, `*`))
    rowSums(pmax(expm1(sweep(deviation, 2, mu, `+`)), 0) * unit)
  })
  n_disease <- length(disease_draws)
  stopifnot(n_disease > 0)
  totals <- lapply(rho_disease, function(rho) {
    z <- sqrt(rho) * matrix(rnorm(n_global), n_global, n_disease) +
      sqrt(1 - rho) * matrix(rnorm(n_global * n_disease), n_global, n_disease)
    joint <- vapply(seq_len(n_disease), function(j)
      as.numeric(quantile(disease_draws[[j]], pnorm(z[, j]), names = FALSE)),
      numeric(n_global))
    q <- quantile(rowSums(joint), c(0.025, 0.5, 0.975), names = FALSE)
    data.frame(rho_disease = rho, lower = q[1], median = q[2], upper = q[3])
  })
  do.call(rbind, totals)
}
