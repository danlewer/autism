library(RColorBrewer)
d <- read.csv('https://raw.githubusercontent.com/danlewer/autism/main/pd/time_trends_sex_6aug2024.csv')

time_male <- d[d$sex == 'm',]
time_female <- d[d$sex == 'f',]

png('time_trends_sex.png', height = 6, width = 6, units = 'in', res = 300)

par(xpd = NA, mar = c(5, 5, 1, 7))
plot(1, type = 'n', xlim = c(2000, 2018), ylim = c(0, 30), xlab = 'Year', ylab = 'Incident diagnoses per 10,000 person-years')
with(time_male, {
  polygon(x = c(year, rev(year)), y = c(lower, rev(upper)), col = add.alpha(cols[1], 0.3), border = NA)
  lines(year, pred, col = cols[1], lwd = 2)
})
with(time_female, {
  polygon(x = c(year, rev(year)), y = c(lower, rev(upper)), col = add.alpha(cols[2], 0.3), border = NA)
  lines(year, pred, col = cols[2], lwd = 2)
})
ys <- seq(20, 30, length.out = 6)
par(xpd = NA)
rect(2020, ys[c(1, 4)], 2022, ys[c(3, 6)], col = add.alpha(cols[1:2], alpha = 0.3), border = NA)
segments(2020, ys[c(2, 5)], x1 = 2022, col = cols[1:2], lwd = 2)
text(2022.5, ys[c(2, 5)], c('Male', 'Female'), adj = 0)

dev.off()
