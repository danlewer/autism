library(haven) # for reading dta
library(data.table)
library(RColorBrewer)
add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)
form2 <- function (x, n = 2) format(round(x, n), digits = n, nsmall = n)

# read data

rd <- function (file, pd = 'any', id = F) {
  d <- read_dta(file, col_select = c('exposure', 'sex', 'age', 'year', '_D', '_Y'))
  names(d) <- c('autism', 'sex', 'age', 'year', 'diag', 'py')
  d$pd <- pd
  d$id <- id
  as.data.table(d)
}

d <- rbind(
  rd("Any_personality_0_irr_rts_yr_No_ID.dta", pd = 'any', id = F),
  rd("Any_personality_1_irr_rts_yr_ID.dta", pd = 'any', id = T),
  rd("Borderline_0_irr_rts_yr_No_ID.dta", pd = 'bpd', id = F),
  rd("Borderline_1_irr_rts_yr_ID.dta", pd = 'bpd', id = T)
)

# recode variables

d[, sex := factor(sex, 1:2, c('m', 'f'))]
d[, id := factor(id, c(F, T), c('noID', 'ID'))]
d[, autism := factor(autism, 0:1, c('control', 'autism'))]
d[, autism_id := paste0(autism, '_', id)]
d[, autism_id := factor(autism_id, c('control_noID', 'control_ID', 'autism_noID', 'autism_ID'))]


# standardised rates

pf <- function (model, nd) {
  f <- model$family$linkinv
  p <- predict(model, newdata = nd, se.fit = T)
  r <- cbind(nd, 
             pred = f(p$fit), 
             lower = f(p$fit - qnorm(0.975) * p$se.fit),
             upper = f(p$fit + qnorm(0.975) * p$se.fit))
  r$form <- paste0(form2(r$pred), ' (', form2(r$lower), '-', form2(r$upper), ')')
  r$form <- gsub('\\(', ' (', gsub(' ', '', r$form))
  r
}

# -----------------------------
# time trends by sex and autism
# -----------------------------

inputs <- expand.grid(sex = c('m', 'f'), autism = c('autism', 'control'), pd = c('any', 'bpd'))
f <- function (SEX, AUTISM, PD) glm(diag ~ poly(year, 2) + poly(age, 2) + offset(log(py)), data = d[sex == SEX & autism == AUTISM & pd == PD], family = 'poisson')
models <- mapply(f, SEX = inputs$sex, AUTISM = inputs$autism, PD = inputs$pd, SIMPLIFY = F)
nd <- data.frame(py = 1, year = 2000:2018, age = 20)
predicted_rates <- lapply(models, pf, nd = nd)
labs <- lapply(split(inputs, 1:nrow(inputs)), function (x) x[rep(1, nrow(nd)),])
predicted_rates <- mapply(cbind, labs, predicted_rates, SIMPLIFY = F)
predicted_rates <- rbindlist(predicted_rates)
predicted_rates$py <- NULL
obs <- rbind(d[, .(pd = 'any', diag = sum(diag), py = sum(py)), c('autism', 'sex', 'year')],
             d[pd == 'bpd', .(pd = 'bpd', diag = sum(diag), py = sum(py)), c('autism', 'sex', 'year')])
predicted_rates <- obs[predicted_rates, on = c('autism', 'sex', 'year', 'pd')]
predicted_rates[, obs := diag / py]
predicted_rates[, pch := ifelse(autism == 'control', 1, 4)]

# predicted_rates <- fread("https://raw.githubusercontent.com/danlewer/autism/main/pd/time_trends_sex_autism_21aug2024.csv")

lf <- function(SEX, PD, AUTISM, col = 'red', ymax = 80) {
  with(predicted_rates[sex == SEX & pd == PD & autism == AUTISM], {
    lines(year, pred, col = col, lwd = 1.5)
    polygon(x = c(year, rev(year)), y = pmin(c(lower, rev(upper)), ymax), col = add.alpha(col, 0.3), border = NA)
    points(year, obs, pch = pch, col = col)
  })
}
plot2 <- function (.) {
  plot(1, type = 'n', xlim = c(2000, 2018), ylim = c(0, 80), axes = F, xlab = NA, ylab = NA)
  rect(2000, 0, 2018, 80)
  axis(1, seq(2000, 2015, 5), pos = 0)
  axis(2, 0:8 * 10, las = 2, pos = 2000)
}

cols <- brewer.pal(3, 'Set1')
par(mfrow = c(2, 2))

plot2()
lf('m', 'any', 'autism', col = cols[1])
lf('m', 'any', 'control', col = cols[2])
text(2009, 70, 'Male\nAny personality disorder\ndiagnosis')

plot2()
lf('m', 'bpd', 'autism', col = cols[1])
lf('m', 'bpd', 'control', col = cols[2])
text(2009, 70, 'Male\nBorderline or emotionally unstable\npersonality disorder')

plot2()
lf('f', 'any', 'autism', col = cols[1])
lf('f', 'any', 'control', col = cols[2])
text(2009, 70, 'Female\nAny personality disorder\ndiagnosis')

plot2() 
lf('f', 'bpd', 'autism', col = cols[1])
lf('f', 'bpd', 'control', col = cols[2])
text(2009, 70, 'Female\nBorderline or emotionally unstable\npersonality disorder')
