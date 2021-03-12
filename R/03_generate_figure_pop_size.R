require(ggplot2)
require(ggpubr)
library(lubridate)
setwd("")


Pop_size <- read.csv("Population_size.csv",header=T)
Pop_size[, "Month_Year"] <-as_date(dmy(Pop_size$Month_Year))
Pop_size$Study_location = factor(Pop_size$Study_location, levels = c("Pescara", "Lviv"), ordered = TRUE) 

facet_wrap_fig = ggplot(Pop_size, aes(x=Month_Year, y=Mean)) +
  geom_point() +
  scale_x_date(breaks = as.Date(c("2018-04-01", "2018-07-01", "2018-10-01", "2019-01-01", "2019-04-01", "2019-07-01")),
               labels = c("Apr-18", "Jul-18", "Oct-18", "*Jan-19", "Apr-19", "Jul-19")) +
  geom_errorbar(aes(ymin=Lower_95CI, ymax=Upper_95CI), linetype=3, width=8) +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  theme_bw() +
  geom_line() +
  xlab("Primary period") +
  ylab("Estimated population size") + 
  facet_grid(Study_location ~ Study_site, scales = "free_y")


pdf("Facet_wrap.pdf", width = 14, height =8)
facet_wrap_fig
dev.off()
