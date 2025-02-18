library(tidyverse)
library(cowplot)
library(growthcurver)
############################################################################
# ADAPTER THIS DEPENDING ON YOUR SOURCE DATA
dates <- read.table("./08-CoCultureExperiments/addingtimes.txt", header = TRUE, sep = "\t")

for (i in 1:nrow(dates)) {
    expdate <- dates[i, 1]
    adding <- (dates[i, 2] + 600) / 3600
############################################################################

data <- read.table(paste0("./08-CoCultureExperiments/",
                          expdate, "_CoCulture/",
                          expdate, "_GrowthCurvesCoCulture.txt"),
                          header = TRUE, sep = "\t")
meta <- read.table(paste0("./08-CoCultureExperiments/",
                          expdate, "_CoCulture/",
                          expdate, "_CoCultureLayout.txt"),
                          header = TRUE, sep = "\t")
counts <- read.table(paste0("./08-CoCultureExperiments/",
                            expdate, "_CoCulture/",
                            expdate, "_countdata_coculture.txt"),
                            header = TRUE, sep = "\t")
data$sec <- data$sec / 3600

# rename wells
meta$pos <- ifelse(meta$column<10, paste0(meta$row, "0", meta$column), paste0(meta$row, meta$column))
newcolnames <- meta$pos[order(meta$pos)]
colnames(data)[4:ncol(data)] <- newcolnames

# remove edge wells from data
pattern <- "01|12|A|H"
datafil <- data[, !grepl(pattern, colnames(data))]
meta <- meta[!grepl(pattern, meta$pos), ]
meta <- meta[meta$bac != "empty", ]
# check well with errors at the beginning (droplets, etc.)
startsit <- datafil[1, 4:ncol(datafil)]
hist(as.numeric(startsit), breaks = 100)

# remove columns with values > 0.15 in the first time point
dfp1 <- datafil[, 1:3]
dfp2 <- datafil[, -c(1:3)]
dfp2 <- dfp2[, dfp2[1, ] <= 0.15]
meta <- meta[meta$pos %in% colnames(dfp2), ]

# get the mean of the blank wells
blank <- dfp2 %>%
    select("G02", "G03")
blank$mean <- rowMeans(blank, na.rm = TRUE)

# substract blank from data
dfp2sub <- dfp2
dfp2sub <- sweep(dfp2sub, 1, blank$mean, "-")


# calculate change of growth in co-culture
counts$changebac <- counts$bac_cfu_co / counts$bac_cfu
counts$changefun <- counts$fun_cfu_co / counts$fun_cfu

# prepare data for plotting
baclist <- unique(meta$bac)
baclist <- baclist[baclist != "none"]
funlist <- unique(meta$fun)
funlist <- funlist[funlist != "none"]
batchcol <- c("#9E001B", "#002D41","#388E3C")

countsred <- counts %>%
    select(bac, fun, bac_cfu, fun_cfu, bac_cfu_co, fun_cfu_co)
countslong <- countsred %>% pivot_longer(cols = -c("bac", "fun"), names_to = "pos", values_to = "count")
countslong$combo <- paste0(countslong$bac, "_", countslong$fun)


countsred <- counts %>%
    select("bac", "bacname", "fun", "funname",  "changebac", "changefun") %>%
    pivot_longer(cols = c("changebac", "changefun"), names_to = "comparison", values_to = "change") %>%
    data.frame()

countsred$baclab <- paste0(countsred$bacname, "\n", countsred$bac)
countsred$funlab <- paste0(countsred$funname, "\n", countsred$fun)



plotlist <- list()
plotlist2 <- list()
plotlist3 <- list()
k = 0
for (i in seq_along(baclist)) {
    for (j in seq_along(funlist)) {
        print(paste0(baclist[i], " vs. ", funlist[j]))
        k = k + 1
        bacpos <- meta$pos[meta$bac == baclist[i] & meta$fun == "none"]
        funpos <- meta$pos[meta$fun == funlist[j] & meta$bac == "none"]
        copos <- meta$pos[meta$bac == baclist[i] & meta$fun == funlist[j]]
        if (length(bacpos) == 2) {
            bacmean <- rowMeans(dfp2sub[, bacpos], na.rm = TRUE)
            bacsd <- apply(dfp2sub[, bacpos], 1, sd, na.rm = TRUE)
        } else {
            bacmean <- dfp2sub[, bacpos]
            bacsd <- rep(0, nrow(dfp2sub))
        }
        if (length(funpos) == 2) {
            funmean <- rowMeans(dfp2sub[, funpos], na.rm = TRUE)
            funsd <- apply(dfp2sub[, funpos], 1, sd, na.rm = TRUE)
        } else {
            funmean <- dfp2sub[, funpos]
            funsd <- rep(0, nrow(dfp2sub))
        }
        if (length(copos) == 2) {
            comean <- rowMeans(dfp2sub[, copos], na.rm = TRUE)
            cosd <- apply(dfp2sub[, copos], 1, sd, na.rm = TRUE)
        } else {
            comean <- dfp2sub[, copos]
            cosd <- rep(0, nrow(dfp2sub))
        }
        comb <- cbind(dfp1, bacmean, funmean, comean)
        comb$additive <- comb$bac + comb$fun
        colnames(comb) <- c("cycle", "Time", "Temp", "bac", "fun", "co", "additive")

        comblong <- comb %>% pivot_longer(cols = c("bac", "fun", "co", "additive"), names_to = "comparison", values_to = "OD")
        comb2 <- cbind(dfp1, bacsd, funsd, cosd)
        colnames(comb2) <- c("cycle", "Time", "Temp", "bac", "fun", "co")
        comblong2 <- comb2 %>% pivot_longer(cols = c("bac", "fun", "co"), names_to = "comparison", values_to = "ODsd")
        combcomb <- left_join(comblong, comblong2, by = c("cycle", "Time", "Temp", "comparison"))
        combcomb <- data.frame(combcomb)
        combcomb$bac <- baclist[i]
        combcomb$fun <- funlist[j]
        maxOD <- max(combcomb$OD, na.rm = TRUE)
        Temp <- combcomb$Temp
        OD_min <- min(combcomb$OD, na.rm = TRUE)
        OD_max <- max(combcomb$OD, na.rm = TRUE)
        temp_min <- 20
        temp_max <- 30
        combcomb$rescaled_temp <- (Temp - temp_min) / (temp_max - temp_min) * (OD_max - OD_min) + OD_min
        combcomb$Time <- combcomb$Time + adding
        force(
        p1 <- ggplot(combcomb[combcomb$comparison!="additive", ], aes(x = Time, y = OD, color = comparison)) +
            geom_line(data = combcomb[combcomb$comparison=="additive", ], aes(x = Time, y = OD), color = "grey", size = 0.5, alpha = 0.6, linetype = "dashed") +
            geom_point(size = 1, alpha = 0.7) +
            geom_errorbar(aes(ymin = OD - ODsd, ymax = OD + ODsd), width = 0.1, alpha = 0.3) +
            #geom_smooth(method = "loess", se = FALSE) +
            theme_bw() +
            scale_color_manual(values = batchcol) +
            theme(legend.position = "none",
                    plot.title = element_text(size = 20),
                    plot.subtitle = element_text(size = 12),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14),
                    axis.text = element_text(size = 12)) +
            labs(title = paste0(baclist[i], " vs. ", funlist[j])) +
            geom_line(aes(x = Time, y = rescaled_temp), color = "grey", linetype = "dotted", alpha = 0.5) +
            scale_y_continuous(
                name = "OD",
                sec.axis = sec_axis(~ (. - OD_min) / (OD_max - OD_min) * (temp_max - temp_min) + temp_min, 
                                    name = "Temperature (°C)")
            ) +
            scale_x_continuous(breaks = seq(round(min(combcomb$Time),0), round(max(combcomb$Time),0), by = 5),
            limits = c(round(min(combcomb$Time),0), round(max(combcomb$Time),0)),
            minor_breaks = seq(round(min(combcomb$Time),0), round(max(combcomb$Time),0), by = 1))
            ) 
        
        plotlist[[k]] <- p1
        p2 <- ggplot(countsred[countsred$bac == baclist[i] & countsred$fun == funlist[j],], aes(x = comparison, y = log2(change), fill = comparison)) +
            geom_bar(stat = "identity") +
            scale_fill_manual(values = batchcol[c(1,3)]) +
            theme_bw() +
            theme(axis.text.x = element_blank(), legend.position = "none") +
            labs(x = "Comparison", y = "log2 fold change")
        plotlist2[[k]] <- p2
           p3 <- ggdraw() + draw_plot(p1) +
     draw_plot(p2, x = 0.1, y = 0.7, width = 0.2, height = 0.20, hjust = 0)
    plotlist3[[k]] <- p3
        ggsave(paste0("./08-CoCultureExperiments/", expdate, "_CoCulture/", expdate, "_GrowthCurvesCoCulture_", baclist[i], "_vs_", funlist[j], "_fp.svg"), p3, width = 12, height = 10)
    }
}

}