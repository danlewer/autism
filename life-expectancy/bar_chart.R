library(RColorBrewer)

d <- data.frame(
  group = rep(c('noID', 'ID'), each = 2),
  sex = rep(c('female', 'male'), 2)
)
d$le <- c(76.84, 74.57, 69.61, 71.66)
d$leC <- c(83.29, 80.71, 84.20, 78.94)
d$leD <- c(6.45, 6.14, 14.59, 7.28)

d$le_lower <- c(72.23, 71.94, 66.04, 68.97)
d$le_upper <- c(81.49, 77.60, 74.27, 75.00)

d$leC_lower <- c(81.32, 79.46, 81.85, 77.82)
d$leC_upper <- c(85.35, 82.11, 86.88, 80.16)

d$leD_lower <- c(1.37, 2.84, 9.45, 3.78)
d$leD_upper <- c(11.58, 9.07, 19.02, 10.27)

cols <- brewer.pal(3, 'Set3')
xl <- c(matrix(1:12, ncol = 3, byrow = T) + 0:3)

png('autism_life_expectancy_bar_chart.png', height = 5, width = 8, units = 'in', res = 300)

par(xpd = NA, mar = c(3, 5, 1, 12))
plot(1, type = 'n', xlim = c(0, 17), ylim = c(0, 100), xlab = NA, ylab = NA, axes = F)
rect(xl, 0, xl+1, c(d$le, d$leC, d$leD), col = rep(cols, each = 4))
arrows(xl + 0.5, c(d$le_lower, d$leC_lower, d$leD_lower), y1 = c(d$le_upper, d$leC_upper, d$leD_upper), code = 3, angle = 90, length = 0.05)
axis(2, 0:10 * 10, pos = 0, las = 2)
rect(0, 0, 17, 100)
segments(8.5, 0, y1 = 100)
text(c(4.25, 12.75), 92, c('Without intellectual\ndisability', 'With intellectual\ndisability'))
text(apply(matrix(xl, ncol = 3), 1, mean) + 0.5, -5, rep(c('Female', 'Male'), 2))
title(ylab = 'Life expectancy at age 18', line = 2)
ys <- c(40, 60, 70)
rect(18, ys, 19, ys+7, col = rev(cols))
text(19.5, ys+5, rev(c('Autistic participants', 'Comparison group\nwithout autism or\nintellectual disability', 'Deficit')), adj = c(0, 1))

dev.off()
