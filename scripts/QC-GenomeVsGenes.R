
# Use Mahalanobis Distance to find those that deviate significantly from the mean (3 sd)
# Mahalanobis Distance generatlisation of Z-score (measure of the distance from mean)

archaea  <- read.csv(file = "data/qc/count/archaea.csv", header = TRUE, row.names = 1)
bacteria <- read.csv(file = "data/qc/count/archaea.csv", header = TRUE, row.names = 1)

archaea$dist = mahalanobis(
    archaea,
    colMeans(archaea),
    cov(archaea)
)
archaea$pvalue = pchisq(
    archaea$dist,
    df=1,
    lower.tail=FALSE
)
bacteria$dist = mahalanobis(
    bacteria,
    colMeans(bacteria),
    cov(bacteria)
)
bacteria$pvalue = pchisq(
    bacteria$dist,
    df=1,
    lower.tail=FALSE
)
