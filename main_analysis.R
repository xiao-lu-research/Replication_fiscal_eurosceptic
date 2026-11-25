##################################################################################################################
# Replication file for the main results in 
# "Rules Under Pressure: Euroscepticism and the Enforcement of EU Fiscal Discipline"
# to be published by Journal of European Public Policy
# edited by Xiao Lu
##################################################################################################################
# Please set to the current working directory.
rm(list=ls())
options(stringsAsFactors = FALSE)

library(ggplot2)
library(lmtest)
library(sandwich)
library(survival)

load("replication_data.RData")

output_dir <- "outputs"
dir.create(output_dir, showWarnings = FALSE)

dat <- replication_data[
  replication_data$year >= 1999L & replication_data$year <= 2020L,
]

model_vars <- c(
  "fiscal_fundamentals_index_0_1",
  "eurobarometer_euroskepticism_fraction",
  "nonoverlap_surveillance_0_1",
  "country_power_size_index_z",
  "dist_comcou",
  "years_to_next_election",
  "governing_party_count",
  "government_left_right_position",
  "government_eu_anti_pro_position",
  "country_real_gdp_growth_pct",
  "eu_growth_control_pct",
  "euro_area_member_binary",
  "post_crisis_2008_2012_binary"
)

sample_vars <- c(
  "country", "year", "outcome_edp_initiation",
  "outcome_upgrade_sanction", model_vars
)
d <- dat[complete.cases(dat[, sample_vars]), ]

# Logit model
fit_logit <- function(outcome, data) {
  fit <- glm(reformulate(model_vars, response = outcome),
             data = data, family = binomial())
  vc <- sandwich::vcovCL(fit, cluster = data$country, type = "HC1")
  ct <- lmtest::coeftest(fit, vcov. = vc)

  list(
    fit = fit,
    vcov = vc,
    table = ct,
    terms = data.frame(
      term = rownames(ct),
      estimate = unname(ct[, 1]),
      std_error = unname(ct[, 2]),
      p_value = unname(ct[, 4]),
      stringsAsFactors = FALSE
    ),
    n = nrow(data),
    events = sum(data[[outcome]]),
    countries = length(unique(data$country)),
    informative_strata = NA_integer_,
    aic = AIC(fit)
  )
}

# Conditional logit model
fit_conditional_logit <- function(outcome, data) {
  fit <- survival::clogit(
    as.formula(paste(
      outcome, "~",
      paste(c(model_vars, "strata(country)", "cluster(country)"),
            collapse = " + ")
    )),
    data = data,
    method = "efron",
    model = TRUE,
    x = TRUE,
    control = survival::coxph.control(iter.max = 100)
  )
  fit_summary <- summary(fit)$coefficients
  informative <- tapply(
    data[[outcome]],
    data$country,
    function(x) length(unique(x)) == 2L
  )
  informative_countries <- names(informative)[informative]

  list(
    fit = fit,
    vcov = vcov(fit),
    terms = data.frame(
      term = rownames(fit_summary),
      estimate = unname(fit_summary[, "coef"]),
      std_error = unname(fit_summary[, "robust se"]),
      p_value = unname(fit_summary[, "Pr(>|z|)"]),
      stringsAsFactors = FALSE
    ),
    n = nrow(data),
    events = sum(data[[outcome]]),
    countries = length(unique(data$country)),
    informative_strata = sum(informative),
    informative_n = sum(data$country %in% informative_countries),
    aic = AIC(fit)
  )
}

opening_annual <- fit_logit("outcome_edp_initiation", d)
stepping_annual <- fit_logit("outcome_upgrade_sanction", d)

opening_fe_data <- d[d$opening_eligible == 1L, ]
opening_fe <- fit_conditional_logit(
  "outcome_edp_initiation",
  opening_fe_data
)

stepping_fe_data <- d[d$open_edp_risk_set == 1L, ]
stepping_fe <- fit_conditional_logit(
  "outcome_upgrade_sanction",
  stepping_fe_data
)

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, digits = digits, format = "f"))
}

star <- function(p) {
  ifelse(
    p < 0.001, "***",
    ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ifelse(p < 0.10, "+", "")))
  )
}

reporting_p_value <- function(term, estimate, std_error, two_sided_p) {
  if (term == "fiscal_fundamentals_index_0_1") {
    return(pnorm(estimate / std_error, lower.tail = FALSE))
  }
  if (term == "eurobarometer_euroskepticism_fraction") {
    return(pnorm(estimate / std_error, lower.tail = TRUE))
  }
  two_sided_p
}

var_labels <- c(
  outcome_edp_initiation = "EDP Opening",
  outcome_upgrade_sanction = "EDP Stepping-up",
  fiscal_fundamentals_index_0_1 = "Fiscal vulnerability",
  eurobarometer_euroskepticism_fraction = "Domestic Euroscepticism",
  nonoverlap_surveillance_0_1 = "Surveillance intensity",
  country_power_size_index_z = "Member-state power",
  dist_comcou = "Commission-Council conflict",
  years_to_next_election = "Years to next election",
  governing_party_count = "Governing-party count",
  government_left_right_position = "Government left-right position",
  government_eu_anti_pro_position = "Government anti-/pro-EU position",
  country_real_gdp_growth_pct = "Domestic GDP growth",
  eu_growth_control_pct = "EU GDP growth",
  euro_area_member_binary = "Euro-area member",
  post_crisis_2008_2012_binary = "Crisis period, 2008--2012"
)

# Figure 1. Timing of EDP opening, stepping-up, and strict escalation.
country_alpha <- sort(unique(as.character(dat$country)), method = "radix")
country_plot_levels <- rev(country_alpha)
dat$country_plot <- factor(dat$country, levels = country_plot_levels)

country_rows <- split(dat, dat$country)
coverage <- data.frame(
  country = country_alpha,
  start_year = vapply(
    country_rows[country_alpha],
    function(x) as.integer(min(x$year)),
    integer(1)
  ),
  end_year = vapply(
    country_rows[country_alpha],
    function(x) as.integer(max(x$year)),
    integer(1)
  ),
  stringsAsFactors = FALSE
)
coverage$country_plot <- factor(coverage$country, levels = country_plot_levels)

opening_only <- dat[
  dat$outcome_edp_initiation == 1L &
    dat$outcome_upgrade_sanction == 0L,
]
stepping_only <- dat[
  dat$outcome_edp_initiation == 0L &
    dat$outcome_upgrade_sanction == 1L,
]
same_year_both <- dat[
  dat$outcome_edp_initiation == 1L &
    dat$outcome_upgrade_sanction == 1L,
]
stepping_only_non_strict <- stepping_only[
  stepping_only$outcome_upgrade_strict_escalation == 0L,
]
stepping_only_strict <- stepping_only[
  stepping_only$outcome_upgrade_strict_escalation == 1L,
]
same_year_both_non_strict <- same_year_both[
  same_year_both$outcome_upgrade_strict_escalation == 0L,
]
same_year_both_strict <- same_year_both[
  same_year_both$outcome_upgrade_strict_escalation == 1L,
]

p_timeline <- ggplot() +
  geom_segment(
    data = coverage,
    aes(
      x = start_year, xend = end_year,
      y = country_plot, yend = country_plot,
      color = "Member-state years in sample"
    ),
    linewidth = 1.15,
    lineend = "round"
  ) +
  geom_point(
    data = opening_only,
    aes(x = year, y = country_plot, shape = "EDP opening"),
    size = 3.25, stroke = 0.75, color = "black"
  ) +
  geom_point(
    data = stepping_only_non_strict,
    aes(x = year, y = country_plot, shape = "EDP stepping-up"),
    size = 2.75, stroke = 0.75, color = "black"
  ) +
  geom_point(
    data = stepping_only_strict,
    aes(x = year, y = country_plot, shape = "Strict escalation"),
    size = 2.75, color = "black"
  ) +
  geom_point(
    data = same_year_both,
    aes(x = year, y = country_plot, shape = "EDP opening"),
    size = 4.0, stroke = 0.75, color = "black"
  ) +
  geom_point(
    data = same_year_both_non_strict,
    aes(x = year, y = country_plot, shape = "EDP stepping-up"),
    size = 2.35, stroke = 0.75, color = "black"
  ) +
  geom_point(
    data = same_year_both_strict,
    aes(x = year, y = country_plot, shape = "Strict escalation"),
    size = 2.35, color = "black"
  ) +
  scale_color_manual(
    values = c("Member-state years in sample" = "grey72")
  ) +
  scale_shape_manual(
    breaks = c("EDP opening", "EDP stepping-up", "Strict escalation"),
    values = c(
      "EDP opening" = 1,
      "EDP stepping-up" = 2,
      "Strict escalation" = 17
    )
  ) +
  scale_x_continuous(
    breaks = c(1999, 2002, 2005, 2008, 2011, 2014, 2017, 2020),
    minor_breaks = 1999:2020,
    limits = c(1998.5, 2020.5),
    expand = c(0, 0),
    sec.axis = dup_axis(name = NULL)
  ) +
  labs(x = "Year", y = NULL, color = NULL, shape = NULL) +
  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(linewidth = 1.4)
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(size = 3.0, stroke = 0.75)
    )
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.28),
    panel.grid.minor.x = element_line(color = "grey94", linewidth = 0.20),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 8.7, color = "black"),
    axis.text.x = element_text(size = 8.5, color = "black"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = grid::unit(1.2, "lines"),
    legend.text = element_text(size = 8.8),
    plot.margin = margin(7, 10, 7, 7)
  )

ggsave(
  file.path(output_dir, "Fig_1.pdf"),
  p_timeline,
  width = 7.4,
  height = 8.2,
  units = "in"
)

# Figure 2. Distributions of the two main explanatory variables.
stage_event_values <- rbind(
  data.frame(
    Stage = "EDP Opening",
    Variable = "Fiscal Vulnerability",
    Value = d$fiscal_fundamentals_index_0_1[d$outcome_edp_initiation == 1L]
  ),
  data.frame(
    Stage = "EDP Stepping-up",
    Variable = "Fiscal Vulnerability",
    Value = d$fiscal_fundamentals_index_0_1[d$outcome_upgrade_sanction == 1L]
  ),
  data.frame(
    Stage = "EDP Opening",
    Variable = "Domestic Euroscepticism",
    Value = d$eurobarometer_euroskepticism_fraction[
      d$outcome_edp_initiation == 1L
    ]
  ),
  data.frame(
    Stage = "EDP Stepping-up",
    Variable = "Domestic Euroscepticism",
    Value = d$eurobarometer_euroskepticism_fraction[
      d$outcome_upgrade_sanction == 1L
    ]
  )
)
stage_event_values$Variable <- factor(
  stage_event_values$Variable,
  levels = c("Fiscal Vulnerability", "Domestic Euroscepticism")
)
stage_event_values$Stage <- factor(
  stage_event_values$Stage,
  levels = c("EDP Opening", "EDP Stepping-up")
)

p_explanatory_desc <- ggplot(
  stage_event_values,
  aes(x = Stage, y = Value)
) +
  geom_violin(
    fill = "grey86", color = "black", linewidth = 0.35,
    trim = FALSE, width = 0.82
  ) +
  geom_boxplot(
    width = 0.16, fill = "white", color = "black",
    linewidth = 0.35, outlier.shape = NA
  ) +
  facet_wrap(~ Variable, scales = "free_y", nrow = 1) +
  scale_x_discrete(labels = c(
    "EDP Opening" = "EDP\nOpening",
    "EDP Stepping-up" = "EDP\nStepping-up"
  )) +
  labs(x = NULL, y = "Observed Value") +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(
      fill = "grey92", color = "black", linewidth = 0.3
    ),
    axis.text.x = element_text(size = 9.0),
    panel.spacing.x = grid::unit(0.8, "lines")
  )

ggsave(
  file.path(output_dir, "Fig_2.pdf"),
  p_explanatory_desc,
  width = 7.4,
  height = 4.0
)

# Table 1. Descriptive statistics.
desc_vars <- c("outcome_edp_initiation", "outcome_upgrade_sanction", model_vars)
desc <- data.frame(
  Variable = unname(var_labels[desc_vars]),
  N = sapply(d[desc_vars], function(x) sum(!is.na(x))),
  Mean = sapply(d[desc_vars], mean),
  SD = sapply(d[desc_vars], sd),
  Min = sapply(d[desc_vars], min),
  Max = sapply(d[desc_vars], max),
  row.names = NULL
)


# Table 2. Main annual and conditional logit models.
coef_labels <- c(
  fiscal_fundamentals_index_0_1 = "Fiscal vulnerability",
  eurobarometer_euroskepticism_fraction = "Domestic Euroscepticism",
  nonoverlap_surveillance_0_1 = "Surveillance intensity",
  country_power_size_index_z = "Member-state power",
  dist_comcou = "Commission-Council conflict",
  years_to_next_election = "Years to next election",
  governing_party_count = "Governing-party count",
  government_left_right_position = "Government left-right position",
  government_eu_anti_pro_position = "Government anti-/pro-EU position",
  country_real_gdp_growth_pct = "Domestic GDP growth",
  eu_growth_control_pct = "EU GDP growth",
  euro_area_member_binary = "Euro-area member",
  post_crisis_2008_2012_binary = "Crisis period, 2008--2012",
  `(Intercept)` = "Constant"
)
table_models <- list(
  opening_annual,
  opening_fe,
  stepping_annual,
  stepping_fe
)

get_term <- function(model, term) {
  out <- model$terms[model$terms$term == term, , drop = FALSE]
  if (!nrow(out)) return(NULL)
  out[1L, , drop = FALSE]
}

coef_cell <- function(model, term) {
  out <- get_term(model, term)
  if (is.null(out)) return("--")
  p <- reporting_p_value(
    term,
    out$estimate,
    out$std_error,
    out$p_value
  )
  paste0(fmt(out$estimate), star(p))
}

se_cell <- function(model, term) {
  out <- get_term(model, term)
  if (is.null(out)) return("")
  paste0("(", fmt(out$std_error), ")")
}

table_rows <- character()
for (term in names(coef_labels)) {
  table_rows <- c(
    table_rows,
    paste0(
      coef_labels[[term]], " & ",
      paste(vapply(table_models, coef_cell, character(1), term = term),
            collapse = " & "),
      " \\\\"
    ),
    paste0(
      " & ",
      paste(vapply(table_models, se_cell, character(1), term = term),
            collapse = " & "),
      " \\\\"
    )
  )
}

table2 <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{\\hl{Main Logit Models of EDP Opening and Stepping-up}}",
  "\\label{tab:main-models}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  " & \\multicolumn{2}{c}{EDP Opening} & \\multicolumn{2}{c}{EDP Stepping-up} \\\\ ",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  " & Annual Logit & Conditional Logit & Annual Logit & Conditional Logit \\\\ ",
  "\\midrule",
  table_rows,
  "\\midrule",
  paste0(
    "Observations & ",
    paste(vapply(table_models, function(x) x$n, integer(1)), collapse = " & "),
    " \\\\"
  ),
  paste0(
    "Events & ",
    paste(vapply(table_models, function(x) x$events, integer(1)), collapse = " & "),
    " \\\\"
  ),
  paste0(
    "Countries & ",
    paste(vapply(table_models, function(x) x$countries, integer(1)), collapse = " & "),
    " \\\\"
  ),
  paste0(
    "AIC & ",
    paste(vapply(table_models, function(x) fmt(x$aic), character(1)), collapse = " & "),
    " \\\\"
  ),
  "\\bottomrule",
  "\\end{tabular}%",
  "}",
  "\\end{table}"
)
writeLines(table2, file.path(output_dir, "Tab_2.tex"))


profiled_probability_data <- function(model, data, outcome, focal, values) {
  varies <- tapply(
    data[[outcome]], data$country,
    function(x) length(unique(x)) == 2L
  )
  prediction_sample <- data[
    data$country %in% names(varies)[varies],
    , drop = FALSE
  ]
  
  beta <- coef(model$fit)
  x <- as.matrix(prediction_sample[, names(beta), drop = FALSE])
  eta <- as.vector(x %*% beta)
  groups <- split(
    seq_len(nrow(prediction_sample)),
    as.character(prediction_sample$country)
  )
  
  country_intercepts <- vapply(groups, function(index) {
    event_count <- sum(prediction_sample[[outcome]][index])
    uniroot(
      function(intercept) {
        sum(plogis(intercept + eta[index])) - event_count
      },
      c(-max(eta[index]) - 40, -min(eta[index]) + 40)
    )$root
  }, numeric(1L))
  intercept <- unname(
    country_intercepts[as.character(prediction_sample$country)]
  )
  
  intercept_gradient <- matrix(
    0, nrow = nrow(prediction_sample), ncol = length(beta),
    dimnames = list(NULL, names(beta))
  )
  for (index in groups) {
    probability <- plogis(intercept[index[1L]] + eta[index])
    weight <- probability * (1 - probability)
    gradient <- -colSums(x[index, , drop = FALSE] * weight) / sum(weight)
    intercept_gradient[index, ] <- matrix(
      gradient, nrow = length(index), ncol = length(beta), byrow = TRUE
    )
  }
  
  do.call(rbind, lapply(values, function(value) {
    profile_matrix <- x
    profile_matrix[, focal] <- value
    probability <- plogis(intercept + as.vector(profile_matrix %*% beta))
    weight <- probability * (1 - probability)
    gradient <- colMeans(
      sweep(profile_matrix + intercept_gradient, 1L, weight, `*`)
    )
    standard_error <- sqrt(as.numeric(
      t(gradient) %*% model$vcov %*% gradient
    ))
    mean_probability <- mean(probability)
    data.frame(
      value = value,
      probability = mean_probability,
      lower = max(0, mean_probability - 1.96 * standard_error),
      upper = min(1, mean_probability + 1.96 * standard_error)
    )
  }))
}

conditional_prediction_data <- function(focal, values) {
  specifications <- list(
    "EDP Opening" = list(
      model = opening_fe,
      data = opening_fe_data,
      outcome = "outcome_edp_initiation"
    ),
    "EDP Stepping-up" = list(
      model = stepping_fe,
      data = stepping_fe_data,
      outcome = "outcome_upgrade_sanction"
    )
  )
  do.call(rbind, lapply(names(specifications), function(stage) {
    specification <- specifications[[stage]]
    result <- profiled_probability_data(
      specification$model,
      specification$data,
      specification$outcome,
      focal,
      values
    )
    result$Stage <- stage
    result
  }))
}
# Figure 3. conditional-logit probabilities by fiscal vulnerability.
pred_fiscal <- conditional_prediction_data(
  "fiscal_fundamentals_index_0_1",
  seq(0, 1, length.out = 101)
)
pred_fiscal$Stage <- factor(
  pred_fiscal$Stage,
  levels = c("EDP Opening", "EDP Stepping-up")
)

p_fs <- ggplot(pred_fiscal, aes(x = value, y = probability)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = Stage),
    alpha = 0.20, color = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.8, color = "black") +
  facet_wrap(~ Stage, nrow = 1) +
  scale_fill_manual(values = c(
    "EDP Opening" = "grey75",
    "EDP Stepping-up" = "grey45"
  )) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Fiscal Vulnerability", y = "Predicted Probability") +
  guides(fill = "none") +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(
      fill = "grey92", color = "black", linewidth = 0.3
    )
  )

ggsave(
  file.path(output_dir, "Fig_3.pdf"),
  p_fs,
  width = 7.4,
  height = 4.2
)

# Figure 4. conditional-logit probabilities by Euroscepticism.
euroscepticism_range <- quantile(
  d$eurobarometer_euroskepticism_fraction,
  c(0.05, 0.95),
  na.rm = TRUE
)
pred_euro <- conditional_prediction_data(
  "eurobarometer_euroskepticism_fraction",
  seq(
    euroscepticism_range[1],
    euroscepticism_range[2],
    length.out = 101
  )
)
pred_euro$Stage <- factor(
  pred_euro$Stage,
  levels = c("EDP Opening", "EDP Stepping-up")
)

p_eu <- ggplot(
  pred_euro,
  aes(x = value, y = probability, linetype = Stage, group = Stage)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = Stage),
    alpha = 0.20, color = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.8, color = "black") +
  scale_linetype_manual(values = c(
    "EDP Opening" = "solid",
    "EDP Stepping-up" = "dashed"
  )) +
  scale_fill_manual(values = c(
    "EDP Opening" = "grey75",
    "EDP Stepping-up" = "grey45"
  )) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x = "Domestic Euroscepticism",
    y = "Predicted Probability",
    linetype = NULL,
    fill = NULL
  ) +
  guides(
    fill = "none",
    linetype = guide_legend(override.aes = list(linewidth = 1.1))
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.key.width = grid::unit(1.8, "cm"),
    legend.spacing.x = grid::unit(0.35, "cm")
  )

ggsave(
  file.path(output_dir, "Fig_4.pdf"),
  p_eu,
  width = 7.0,
  height = 4.2
)
