# plot subclass groups in Workplaces
library(tidyverse)
library(data.table)
library(gridExtra)
library(grid)

load(file = "~/ClimBehav_bj/results/data/processed/pClass.RData")

p <- grid.arrange(pClassSub[[7]], pClassSub[[8]], pClassSub[[9]], pClassSub[[10]], nrow = 1)

#save as eps
cairo_ps(file = "~/ClimBehav_bj/results/plots/figureS2.eps",width = 16,height = 6)
pushViewport(viewport(layout = grid.layout(1,1)))
print(ggdraw(p), vp = viewport(layout.pos.row = 1,layout.pos.col = 1))
dev.off()