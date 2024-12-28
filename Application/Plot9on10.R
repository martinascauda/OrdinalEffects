load("Application/DataInt910_For_Plots.RData")
library(ggplot2)
library(ggrain)
#Effects of changing variable 9 to the four levels of variable 10 
ggplot(data, aes(x = Level, y = OCE, fill = 	Change)) +
  geom_rain(alpha = .5, rain.side = 'f', cov = "Change", 
            boxplot.args = list(outlier.shape = NA, alpha = .8),
            violin.args = list(alpha = .5, color = NA),
            point.args=list(size=0.2),
            point.args.pos = list(position = position_jitterdodge(
              jitter.width = 0.18,
              jitter.height = 0,
              dodge.width = 0.5,
              seed = 42
            )),
            boxplot.args.pos = list(width = .1,
                                    position = ggpp::position_dodgenudge(width = 0.5,
                                                                         x = 0.045)),
            violin.args.pos = list(width = 0.9,
                                   position = position_nudge(x = c(rep(rep(-.3, 256*2),12), rep(rep(.3, 256*2),12)
                                   ))
            )) +
  theme_classic() +
  scale_y_continuous(
    breaks = seq(-0.6,0.6, by = 0.1), labels = scales::comma # Adjust break interval
  ) +
  stat_summary(fun = mean, geom = "line", aes(group = Change, color = Change),
               position = position_dodge(width = 0.5), linewidth = 0.3, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", aes(group = Change, color = Change),
               position = position_dodge(width = 0.5), size = 1.5, pch = 23, alpha = 1, color = "black") +
  scale_fill_manual(name="Intervention \non Anhedonia ",
                    labels=c(expression(1 %->% 2),expression(1 %->% 3),expression(1 %->% 4), expression(2 %->% 3 ), expression(2 %->% 4), expression(3 %->% 4)), values=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02" ,"#A6761D","#666666"))  +
  scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02" ,"#A6761D","#666666"))+
  guides(
    fill = guide_legend(position = "inside"), color = 'none')+
  theme(legend.direction="horizontal", legend.position.inside = c(0.25, 0.9), #legend.key.size = unit(1.5, 'cm'), legend.title = element_text(size = 15),
        #legend.text = element_text(size = 10)
        )+
  labs(x="Ordinal Levels of Fatigue", y = "Ordinal Causal Effect (OCE)")

ggsave(file="Application/Effect910.pdf", width = 210, height = 149.5, units = "mm")
