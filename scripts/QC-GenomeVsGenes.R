
# Use Mahalanobis Distance to find those that deviate significantly from the mean (3 sd)
# Mahalanobis Distance generatlisation of Z-score (measure of the distance from mean)

# Pulled this from: https://stackoverflow.com/a/39051803/3856359
requiredPackages = c('dplyr','ggplot2','ggrepel', 'magrittr')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

process = function(path) {
    data <- read.csv(file = path, header = TRUE, row.names = 1)
    data$dist = mahalanobis(
        data,
        colMeans(data),
        cov(data)
    )
    data$pvalue = pchisq(
        data$dist,
        df=1,
        lower.tail=FALSE
    )
    data$id = row.names(data)

    data %>% mutate(color = case_when(
        dist > 3 ~ "Outliers",
        T ~ "Values"
    )) %>% ggplot(aes(x = size/10^6, y = count, color = factor(color))) +
        coord_cartesian() +
        geom_point() +
        geom_text_repel(data = . %>% filter(dist > 3), aes(label = id)) +
        labs(
            x = "Genome Size (Gbp)",
            y = "Frequency of Genes",
            color = "Legend"
        )
}

process("data/qc/count/archaea.csv")
ggsave(
    "plot/qc/count/archaea.png",
    width = 16,
    height = 9
)

process("data/qc/count/bacteria.csv")
ggsave(
    "plot/qc/count/bacteria.png",
    width = 16,
    height = 9
)
