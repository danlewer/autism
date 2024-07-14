library(RColorBrewer)
add.alpha <- function(cols, alpha) rgb(t(col2rgb(cols)/255), alpha = alpha)

# read data

d <- read.csv('https://raw.githubusercontent.com/danlewer/autism/main/pd/time_trends.csv')
cols <- brewer.pal(3, 'Set1')[1:2]

# panel chart

png('time_trends_panel.png', height = 6, width = 12, units = 'in', res = 300)

par(mfrow = c(1, 2), mar = c(5, 5, 1, 0), oma = c(0, 0, 0, 10), xpd = F)

plot(1, type = 'n', xlim = c(2000, 2018), ylim = c(0, 60), xlab = 'Year', ylab = 'Incident diagnoses per 10,000 person-years')
with(d[d$pd_type == 'any' & d$group == 'autistic',], {
  polygon(x = c(year, rev(year)), y = c(lower, rev(upper)), col = add.alpha(cols[1], 0.3), border = NA)
  lines(year, pred, col = cols[1], lwd = 2)
})
with(d[d$pd_type == 'any' & d$group == 'control',], {
  polygon(x = c(year, rev(year)), y = c(lower, rev(upper)), col = add.alpha(cols[2], 0.3), border = NA)
  lines(year, pred, col = cols[2], lwd = 2)
})
text(2009, 58, 'Any personality\ndisorder diagnosis')

plot(1, type = 'n', xlim = c(2000, 2018), ylim = c(0, 60), xlab = 'Year', ylab = NA)
with(d[d$pd_type == 'bpd' & d$group == 'autistic',], {
  polygon(x = c(year, rev(year)), y = c(lower, rev(upper)), col = add.alpha(cols[1], 0.3), border = NA)
  lines(year, pred, col = cols[1], lwd = 2)
})
with(d[d$pd_type == 'bpd' & d$group == 'control',], {
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
