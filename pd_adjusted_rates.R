library(haven) # for reading dta
library(data.table)
library(RColorBrewer)
add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)
form2 <- function (x, n = 2) format(round(x, n), digits = n, nsmall = n)

# setwd("C:/Users/rmhidle/OneDrive - University College London/liz/pd/adjusted_rates")

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

m1 <- glm(diag ~ autism_id + sex + poly(age, 2) + poly(year, 2) + offset(log(py)), data = d[pd == 'any'], family = 'poisson')
m2 <- glm(diag ~ autism_id + sex + poly(age, 2) + poly(year, 2) + offset(log(py)), data = d[pd == 'bpd'], family = 'poisson')

nd <- CJ(autism_id = c('control_noID', 'control_ID', 'autism_noID', 'autism_ID'), sex = c('m', 'f'), py = 1, year = 2009, age = 20)
stand_any <- pf(m1, nd = nd)
stand_bpd <- pf(m2, nd = nd)

stand_any$pd <- 'any'
stand_bpd$pd <- 'bpd'
stand <- rbind(stand_any, stand_bpd)
stand[, autism_id := factor(autism_id, c('autism_noID', 'autism_ID', 'control_noID', 'control_ID'))]
stand[, sex := factor(sex, c('m', 'f'))]
stand <- stand[order(autism_id, sex, pd)]
write.csv(stand, 'standardised_rates.csv')

# ------------------------------------
# time trends in PD autism vs. control
# ------------------------------------ 

m1 <- glm(diag ~ poly(year, 2) + poly(age, 2) + offset(log(py)), data = d[autism == 'autism' & pd == 'any'], family = 'poisson')
m2 <- glm(diag ~ poly(year, 2) + poly(age, 2) + offset(log(py)), data = d[autism == 'control' & pd == 'any'], family = 'poisson')
nd <- CJ(py = 1, year = 2000:2018, age = 20)
time_autism_any <- pf(m1, nd = nd)
time_control_any <- pf(m2, nd = nd)

# sanity check

d[autism == 'autism' & pd == 'any', .(cases = sum(diag), py = sum(py), rate = sum(diag) / sum(py)), year]
d[autism == 'control' & pd == 'any', .(cases = sum(diag), py = sum(py), rate = sum(diag) / sum(py)), year]

# plot

cols <- brewer.pal(3, 'Set1')[1:2]

png('time_trends_any.png', height = 6, width = 8, units = 'in', res = 300) 

par(mar = c(5, 5, 1, 10), xpd = F)
plot(1, type = 'n', xlim = c(2000, 2018), ylim = c(0, 60), xlab = 'Year', ylab = 'Incident PD diagnoses per 10,000 person-years')
with(time_autism_any, {
  polygon(x = c(year, rev(year)), y = c(lower, rev(upper)), col = add.alpha(cols[1], 0.3), border = NA)
  lines(year, pred, col = cols[1], lwd = 2)
})
with(time_control_any, {
  polygon(x = c(year, rev(year)), y = c(lower, rev(upper)), col = add.alpha(cols[2], 0.3), border = NA)
  lines(year, pred, col = cols[2], lwd = 2)
})
ys <- seq(20, 30, length.out = 6)
par(xpd = NA)
rect(2020, ys[c(1, 4)], 2022, ys[c(3, 6)], col = add.alpha(cols[2:1], alpha = 0.3), border = NA)
segments(2020, ys[c(2, 5)], x1 = 2022, col = cols[2:1], lwd = 2)
text(2022.5, ys[c(2, 5)], c('Comparison\ngroup', 'Autistic\ngroup'), adj = 0)

dev.off()

# BPD

m1 <- glm(diag ~ poly(year, 2) + poly(age, 2) + offset(log(py)), data = d[autism == 'autism' & pd == 'bpd'], family = 'poisson')
m2 <- glm(diag ~ poly(year, 2) + poly(age, 2) + offset(log(py)), data = d[autism == 'control' & pd == 'bpd'], family = 'poisson')
time_autism_bpd <- pf(m1, nd = nd)
time_control_bpd <- pf(m2, nd = nd)

png('time_trends_bpd.png', height = 6, width = 8, units = 'in', res = 300) 

par(mar = c(5, 5, 1, 10), xpd = F)
plot(1, type = 'n', xlim = c(2000, 2018), ylim = c(0, 60), xlab = 'Year', ylab = 'Incident PD diagnoses per 10,000 person-years')
with(time_autism_bpd, {
  polygon(x = c(year, rev(year)), y = c(lower, rev(upper)), col = add.alpha(cols[1], 0.3), border = NA)
  lines(year, pred, col = cols[1], lwd = 2)
})
with(time_control_bpd, {
  polygon(x = c(year, rev(year)), y = c(lower, rev(upper)), col = add.alpha(cols[2], 0.3), border = NA)
  lines(year, pred, col = cols[2], lwd = 2)
})
ys <- seq(20, 30, length.out = 6)
par(xpd = NA)
rect(2020, ys[c(1, 4)], 2022, ys[c(3, 6)], col = add.alpha(cols[2:1], alpha = 0.3), border = NA)
segments(2020, ys[c(2, 5)], x1 = 2022, col = cols[2:1], lwd = 2)
text(2022.5, ys[c(2, 5)], c('Comparison\ngroup', 'Autistic\ngroup'), adj = 0)

dev.off()

# panel chart

png('time_trends_panel.png', height = 6, width = 12, units = 'in', res = 300)

par(mfrow = c(1, 2), mar = c(5, 5, 1, 0), oma = c(0, 0, 0, 10), xpd = F)

plot(1, type = 'n', xlim = c(2000, 2018), ylim = c(0, 60), xlab = 'Year', ylab = 'Incident diagnoses per 10,000 person-years')
with(time_autism_any, {
  polygon(x = c(year, rev(year)), y = c(lower, rev(upper)), col = add.alpha(cols[1], 0.3), border = NA)
  lines(year, pred, col = cols[1], lwd = 2)
})
with(time_control_any, {
  polygon(x = c(year, rev(year)), y = c(lower, rev(upper)), col = add.alpha(cols[2], 0.3), border = NA)
  lines(year, pred, col = cols[2], lwd = 2)
})
text(2009, 58, 'Any personality\ndisorder diagnosis')

plot(1, type = 'n', xlim = c(2000, 2018), ylim = c(0, 60), xlab = 'Year', ylab = NA)
with(time_autism_bpd, {
  polygon(x = c(year, rev(year)), y = c(lower, rev(upper)), col = add.alpha(cols[1], 0.3), border = NA)
  lines(year, pred, col = cols[1], lwd = 2)
})
with(time_control_bpd, {
  polygon(x = c(year, rev(year)), y = c(lower, rev(upper)), col = add.alpha(cols[2], 0.3), border = NA)
  lines(year, pred, col = cols[2], lwd = 2)
})
text(2009, 58, 'Borderline or emotionally\nunstable personality disorder')

ys <- seq(40, 60, length.out = 6)
par(xpd = NA)
rect(2020, ys[c(1, 4)], 2022, ys[c(3, 6)], col = add.alpha(cols[2:1], alpha = 0.3), border = NA)
segments(2020, ys[c(2, 5)], x1 = 2022, col = cols[2:1], lwd = 2)
text(2022.5, ys[c(2, 5)], c('Comparison\ngroup', 'Autistic\ngroup'), adj = 0)

dev.off()

# table for supplementary

time_autism_any[, pd_type := 'any']
time_autism_bpd[, pd_type := 'bpd']
time_control_any[, pd_type := 'any']
time_control_bpd[, pd_type := 'bpd']

time_autism_any[, group := 'autistic']
time_autism_bpd[, group := 'autistic']
time_control_any[, group := 'control']
time_control_bpd[, group := 'control']

sup <- rbind(time_autism_any, time_autism_bpd, time_control_any, time_control_bpd)
sup <- sup[, c('group', 'pd_type', 'year', 'form')]
fwrite(sup, 'supplementary_table_time_trends.csv')
