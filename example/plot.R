library(ggplot2)
library(relater)

coal <- read.coal("example.coal")
p <- ggplot(coal) + geom_step(aes(x = epoch.start, y = 0.5/haploid.coalescence.rate)) + scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10")

ggsave(p, file = "plot.pdf", width = 10, height = 10)
