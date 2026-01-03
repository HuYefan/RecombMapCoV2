library(readxl)
library(writexl)
library(tidyverse)
library(cowplot)

Sys.setlocale('LC_TIME', 'C')
DateLabel<- format(Sys.Date(), "%y%m%d")

# 加载recombination数据框
recombination<-read_excel("Data/Recombinant_Variants.xlsx",col_names =TRUE)

# 数据预处理
recombination_final <- recombination %>%
  separate_rows(Breakpoint, sep = ",") %>%
  filter(Breakpoint != "") %>%
  separate(Breakpoint, into = c("BP_start", "BP_end"), sep = "-", convert = TRUE) %>%
  mutate(across(c(BP_start, BP_end), as.numeric)) 

recombination_final$Lineage <- factor(recombination_final$Lineage,
                                      level=rev(unique(recombination$Lineage)))
# 图表1: 按月统计Lineage的count
# 数据预处理 - 提取月份
monthly_counts <- recombination %>%
  mutate(Date = as.Date(Date),  # 确保Date是日期类型
         YearMonth = floor_date(Date, "month")) %>%  # 提取年月(每月第一天)
  group_by(YearMonth) %>%
  summarise(Lineage_count = n_distinct(Lineage))  # 计算每月独特Lineage数量

# 绘制月度统计图
monthlyplot<-ggplot(monthly_counts, aes(x = YearMonth, y = Lineage_count)) +
  geom_col(fill = "#4DBBD5", alpha = 0.8) +  # 使用柱状图
  geom_text(aes(label = Lineage_count), vjust = -0.5, size = 6, size.unit = "pt") +  # 添加数值标签
  labs(title = "",
       x = "",
       y = "Number of Designation") +
  theme_bw() +
  scale_y_continuous(limits = c(0, 13),
                     breaks = seq(1:7)*2-2)+
  scale_x_date(date_labels = "%Y-%m", 
               limits = as.Date(c("2019-12-26","2025-12-31")),
               date_breaks = "3 month") +  # 按月显示
  theme(text = element_text(size=7),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank())


# 图表2: Breakpoint范围图
BPcovplot<-ggplot(recombination_final, aes(x = BP_start, xend = BP_end, 
                               y = Lineage, yend = Lineage)) +
  geom_segment(color = "#4DBBD5", linewidth = 1) +
  labs(title = "",
       x = "",
       y = "Lineage") +
  theme_bw() +
  scale_x_continuous(limits = c(1, 29903),
                     breaks = c(1, 5000, 10000, 15000, 20000, 25000, 29903))+
  theme(text = element_text(size=7),
        panel.grid.major.x = element_line(color = "grey80"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank()
        )


# 创建包含所有基因组位置的数据框
genome_positions <- data.frame(Position = 1:29903)  # 假设基因组长度为29903bp

# 统计每个位置的Breakpoint覆盖次数
position_counts <- genome_positions %>%
  mutate(
    BP_count = sapply(Position, function(pos) {
      sum(recombination_final$BP_start <= pos & recombination_final$BP_end >= pos)
    })
  )

# 图表3: 基础折线图
BPstatplot<-ggplot(position_counts, aes(x = Position, y = BP_count)) +
  geom_line(linewidth = 0.5, color = "#4DBBD5") +
  #scale_color_gradient(low = "#91D1C299", high = "#00A087", name = "Count") +
  labs(title = "",
       x = "",
       y = "Breakpoint\n Coverage Count") +
  theme_bw() +
  scale_x_continuous(limits = c(1, 29903),
                     breaks = c(1, 5000, 10000, 15000, 20000, 25000, 29903))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # 避免y轴从0开始
  theme(text = element_text(size=7),
        legend.position = "none",
        panel.grid.minor = element_blank())

# 定义新冠病毒ORF区域 (参考SARS-CoV-2基因组注释)
orf_regions <- read_xlsx("Data/SARS-CoV-2_ORF.xlsx")

# 计算每个ORF区域的平均重组次数
RegionStats <- function(Regions){
  Regions %>%
    rowwise() %>%
    mutate(
      # 获取该ORF区域内所有位置的BP_count
      BP_data = list(position_counts %>% 
                       filter(Position >= Start & Position <= End) %>% 
                       pull(BP_count)),
      # 计算统计量
      Mean_BP = mean(BP_data),
      SD_BP = sd(BP_data),
      SE_BP = SD_BP / sqrt(length(BP_data)),  # 标准误差
      Length_bp = End - Start + 1,
      Total_BP = sum(BP_data),
      BP_density = Total_BP / Length_bp
    ) %>%
    ungroup() %>%
    mutate(ORF = factor(ORF, levels = ORF))   # 计算每bp的平均重组次数
}
orf_stats <- RegionStats(orf_regions)

# 初始化结果数据框
test_results <- data.frame(
  ORF = character(),
  p_value = numeric(),
  significant = logical(),
  stringsAsFactors = FALSE
)
average <- orf_stats[16,]
# 对每个ORF区域执行t检验
for (i in 1:(nrow(orf_stats))) {
  orf <- orf_stats[i, ]
  
  # 假设ORF区域内的重组事件服从正态分布
  # 执行双样本t检验（比较ORF区域与基因组平均水平）
  # 这里使用近似方法，因为我们知道均值、标准差和样本量
  wilcox_test <- wilcox.test(
    x = as.vector(orf$BP_data)[[1]],
    y = as.vector(average$BP_data)[[1]],
    var.equal = FALSE  # 使用Welch校正
  )
  
  test_results <- rbind(test_results, data.frame(
    ORF = as.character(orf$ORF),
    p_value = wilcox_test$p.value,
    significant = FALSE
  ))
}

# 应用Benjamini-Hochberg校正控制错误发现率
test_results$adjusted_p <- p.adjust(test_results$p_value, method = "BH")
test_results$significant <- test_results$adjusted_p < 0.05



# 图表4: 绘制带误差棒的柱状图 (使用标准误差)
ORFstatplot<-ggplot(orf_stats, aes(x = ORF, y = Mean_BP, fill = ORF)) +
  geom_col(fill = "#4DBBD5", color="black", alpha = 0.8, linewidth = 0.3) +
  geom_errorbar(aes(ymin = Mean_BP - SD_BP, ymax = Mean_BP + SD_BP),
                width = 0.3, color = "black", linewidth = 0.3) +
  labs(title = "",
       x = "",
       y = "Mean Recombination\n Events per base") +
  theme_bw() +
  theme(text = element_text(size=7),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

##########
Plot_BP<-plot_grid(
  plot_grid(BPcovplot,
            BPstatplot,
            align = 'hv',
            ncol = 1, 
            nrow = 2, 
            rel_heights = c(5, 0.6)
  ), 
  plot_grid(monthlyplot, 
            ORFstatplot,
            align = 'hv',
            ncol = 2, 
            nrow = 1, 
            rel_widths = c(1.2, 1)
            ), 
  align = 'hv',
  ncol = 1,
  rel_heights = c(5.6,1)
)



ggsave(paste0("Fig-recombinant_",DateLabel,".pdf"), 
       plot = Plot_BP, width = 17,height = 35, units = "cm")

# 定义新冠病毒Spike区域 (参考SARS-CoV-2基因组注释)
spike_regions <- read_xlsx("Data/SARS-CoV-2_Spike.xlsx")
spike_stats <- RegionStats(spike_regions)
spike_regions$ORF<-factor(spike_regions$ORF, levels = rev(spike_regions$ORF))
SpikeBPcovplot<-ggplot(spike_regions, 
       aes(x = Start, xend = End,
           y = ORF,   yend = ORF)) +
  geom_segment(color = "#0F6FC6", linewidth = 1) +
  labs(title = "",
       x = "",
       y = "Region") +
  theme_bw() +
  scale_x_continuous(limits = c(21563, 25384),
                     breaks = c(21563, 25384, 21500+seq(1:7)*500))+
  theme(text = element_text(size=7),
        panel.grid.major.x = element_line(color = "grey80"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank()
  )

SpikeBPstatplot<-ggplot(position_counts, aes(x = Position, y = BP_count)) +
  geom_line(linewidth = 0.5, color = "#0F6FC6") +
  #scale_color_gradient(low = "#91D1C299", high = "#00A087", name = "Count") +
  labs(title = "",
       x = "",
       y = "Breakpoint\n Coverage Count") +
  theme_bw() +
  scale_x_continuous(limits = c(21563, 25384),
                     breaks = c(21563, 25384, 21500+seq(1:7)*500))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # 避免y轴从0开始
  theme(text = element_text(size=7),
        legend.position = "none",
        panel.grid.minor = element_blank())

Spikestatplot<-ggplot(spike_stats, aes(x = ORF, y = Mean_BP, fill = ORF)) +
  geom_col(fill = "#0F6FC6", color="black", alpha = 0.8, linewidth = 0.3) +
  geom_errorbar(aes(ymin = Mean_BP - SD_BP, ymax = Mean_BP + SD_BP),
                width = 0.3, color = "black", linewidth = 0.3) +
  labs(title = "",
       x = "",
       y = "Mean Recombination\n Events per base") +
  theme_bw() +
  theme(text = element_text(size=7),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

Plot_Spike<-plot_grid(
  SpikeBPstatplot,
  SpikeBPcovplot,
  Spikestatplot,
  align = 'hv',
  ncol = 1,
  rel_heights = c(2,1,2)
)
ggsave(paste0("Fig-spike-recombinant_",DateLabel,".pdf"), 
       plot = Plot_Spike, width = 10,height = 10, units = "cm")
