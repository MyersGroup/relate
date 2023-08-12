library(ggplot2)
library(relater)
library(cowplot)

coal <- read.coal("example_bypop.pairwise.coal")
p1 <- ggplot(coal) + geom_step(aes(x = epoch.start, y = 0.5/haploid.coalescence.rate)) + scale_x_continuous(trans = "log10") + scale_y_continuous(limit = c(0,5e4)) + theme_bw()

coal <- read.coal("example_bypop_mrate.pairwise.coal")
p2 <- ggplot(coal) + geom_step(aes(x = epoch.start, y = 0.5/haploid.coalescence.rate)) + scale_x_continuous(trans = "log10") + scale_y_continuous(limit = c(0,5e4)) + theme_bw()

p <- plot_grid(p1,p2,ncol = 1)
ggsave(p, file = "plot.pdf", width = 10, height = 10)
