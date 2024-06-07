library(ggpubr)
data1 <- read.csv('CD14.csv')
data2 <- read.csv('CD16.csv')
data3 <- read.csv('NK proliferating.csv')

p_value <- compare_means(Proportion ~ Time, data = data1[1:21,],method = "wilcox.test")

# CD14和CD16占比变化箱线图
p1 <- ggpaired(data1, x = "Time", y = "Proportion",id="donor",point.size = 1.4,xlab = "Time", ylab = "CD14 Mono Proportion",
              color = "Time", palette = c("#1663A9", "#D2352C","#369F2D"),
              line.color = "gray", line.size = 0.4,font.label = list(size = 9, color = "black"),
              facet.by = "model", short.panel.labs = FALSE)
p2 <- ggpaired(data2, x = "Time", y = "Proportion",id="donor",point.size = 1.4,xlab = "Time", ylab = "CD16 Mono Proportion",
               color = "Time", palette = c("#1663A9", "#D2352C","#369F2D"),
               line.color = "gray", line.size = 0.4,font.label = list(size = 9, color = "black"),
               facet.by = "model", short.panel.labs = FALSE)
p3 <- ggpaired(data3, x = "Time", y = "Proportion",id="donor",point.size = 1.4,xlab = "Time", ylab = "NK proliferating Proportion",
               color = "Time", palette = c("#1663A9", "#D2352C","#369F2D"),
               line.color = "gray", line.size = 0.4,font.label = list(size = 9, color = "black"),
               facet.by = "model", short.panel.labs = FALSE)

p1 + stat_compare_means(label = "p.format", paired = TRUE)


ggsave('CD14_proportion.png',plot = p1,dpi = 600,width = 180,height = 150, units = "mm")
ggsave('CD16_proportion.png',plot = p2,dpi = 600,width = 180,height = 150, units = "mm")
ggsave('NK_proliferating_proportion.png',plot = p3,dpi = 600,width = 180,height = 150, units = "mm")
