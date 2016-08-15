#read data
df = read.csv("...path_to_file/study_data.csv")

#helper functions            
mean_and_sd = function(list)
{
    mean = round(mean(list), 2)
    sd = round(sd(list), 2)
    output = paste(mean, ' (', sd, ')')
    return(output)
}

median_and_iqr = function(list)
{
    median = median(list)
    q1 = quantile(list)[2]
    q3 = quantile(list)[4]
    output = paste(median, ' (' , q1, ',', q3, ')')
    return(output)
}

#check data: there should be 180 observations of 12 variables
dim(df)
names(df)
str(df)

#convert "inoc_conc" variable from numeric to a factor variable; this variable should have 5 levels
df$inoc_conc = as.factor(df$inoc_conc)

#correct concentration notations
levels(df$inoc_conc) = c(levels(df$inoc_conc), c('3e+03'))
df$inoc_conc[df$inoc_conc %in% c('3000')] = '3e+03'
df = droplevels(df)

#identify and remove exclusions: 10 trials should be excluded
excluded = df[!df$comment %in% c(""), ]
excluded = droplevels(excluded)
table(excluded$comment)

#keep only included trials: 170 trials should be included
df = df[df$comment %in% c(""), ]
print("Included trials")
dim(df)

#microorganism distribution
table(df$org_name)
table(df$org_detail)

#inoculum concentration distribution
table(df$inoc_conc)

#actual colony counts from the penultimate dilution
col_counts = df[df$inoc_conc == "30", ]
col_counts$orig_2mcf_conc = col_counts$colony_count_raw * 10000 * 20 / 1000000 #calculate micro conc in millions per mL, accounting for two 100-fold dilutions and that 50 mcL was plated.

round(tapply(col_counts$orig_2mcf_conc, col_counts$org_name, mean), 2) #mean of the concentrations
round(tapply(col_counts$orig_2mcf_conc, col_counts$org_name, sd), 2) #sd of the concentrations

#identify trials available for time to detection comparison
ttd = df[!df$hours_diff %in% c(NA), ]
    
#create subsets of inoc conc for comparisons of time to detection
subset1 = ttd[ttd$inoc_conc %in% c("0.6"), ]
subset2 = ttd[ttd$inoc_conc %in% c("30"), ]
subset3 = ttd[ttd$inoc_conc %in% c("3e+03"), ]
subset4 = ttd[ttd$inoc_conc %in% c("3e+05"), ]

#generating table_2 basic
a = sapply(subset1[, c("csa_hours", "bact_hours", "hours_diff")], mean_and_sd)
b = sapply(subset2[, c("csa_hours", "bact_hours", "hours_diff")], mean_and_sd)
c = sapply(subset3[, c("csa_hours", "bact_hours", "hours_diff")], mean_and_sd)
d = sapply(subset4[, c("csa_hours", "bact_hours", "hours_diff")], mean_and_sd)
table_2 = as.data.frame(rbind(a, b, c, d))
rownames(table_2) = c('0.6', '30', '3e+3', '3e+5')
table_2

#adding p-values to table_2
p1 = t.test(subset1$csa_hours, subset1$bact_hours, paired = T)
p2 = t.test(subset2$csa_hours, subset2$bact_hours, paired = T)
p3 = t.test(subset3$csa_hours, subset3$bact_hours, paired = T)
p4 = t.test(subset4$csa_hours, subset4$bact_hours, paired = T)
p_value = c(round(p1$p.value, 4), round(p2$p.value, 4))
table_2 = cbind(table_2, p_value)
table_2

#csa and bact time to detection for subset1
round(tapply(subset1$csa_hours, subset1$org_name, mean), 1)
round(tapply(subset1$bact_hours, subset1$org_name, mean), 1)
