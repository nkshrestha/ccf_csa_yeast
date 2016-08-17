#read data
df = read.csv("C:/users/shrestn/dropbox/research/csa_yeast_pilot/data/study_data.csv")

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

actual_conc = function(orig_conc, dilution)
{
    if (dilution == "100"){
        conc = orig_conc / 100
    }
    else if (dilution == "10K"){
        conc = orig_conc / 10000
    }
    else if (dilution == "500K"){
        conc = orig_conc / 500000
    }
    else if (dilution == "1"){
        conc = orig_conc
    }
    else if (dilution == "0"){
        conc = NA
    }
        
    return(conc)
}

#check data: there should be 180 observations of 12 variables
dim(df)
names(df)
str(df)

#identify and remove exclusions: 10 trials should be excluded
excluded = df[!df$comment %in% c(""), ]
excluded = droplevels(excluded)
table(excluded$comment)

#assign a specimen no for all dilutions of the same specimen
df$spec_no = as.numeric(substr(df$study_no, 1, 4))

#keep only included trials: 170 trials should be included
df = df[df$comment %in% c(""), ]
print("Included trials")
dim(df)

#microorganism distribution
table(df$org_name)
table(df$org_detail)

#dilution distribution
table(df$dilution)
table(df$spec_no)

#determine actual colony counts from the penultimate dilution
col_counts = df[df$dilution == "10K", c("date", "spec_no", "org_name", "org_detail", "colony_count_raw")]
col_counts$orig_2mcf_conc = col_counts$colony_count_raw * 10000 * 20  #calculate micro conc per mL, accounting for two 100-fold dilutions and that 50 mcL was plated.

df = merge(df, col_counts[ , c("spec_no", "orig_2mcf_conc")], all = TRUE)

#calculate the mean colony count in the original 2mcf suspension for each yeast species
mean_cc = tapply(col_counts$orig_2mcf_conc, col_counts$org_name, mean)
mean_cc

#identify the specimens for which there is no colony count: there should be 3
no_conc = unique(df[df$orig_2mcf_conc %in% c(NA) & df$org_name != "none", c("spec_no", "org_name", "orig_2mcf_conc")])
no_conc

#impute the colony count for these missing values
for (i in 1:nrow(no_conc)){
    no_conc$orig_2mcf_conc[i] = mean_cc[row.names(mean_cc) == no_conc$org_name[i]]
}

no_conc

#merge these imputed values into the study data
temp = merge(df[df$spec_no %in% no_conc$spec_no, !names(df) %in% c("orig_2mcf_conc")], no_conc, by = c("spec_no", "org_name"))
names(temp)[names(temp) == "orig_2mcf_conc"] = "imputed"

df = merge(df, temp[ c("study_no", "imputed")], by = "study_no", all = TRUE)
df$orig_2mcf_conc[df$orig_2mcf_conc %in% c(NA)] = df$imputed[df$orig_2mcf_conc %in% c(NA)]

df = df[ , !names(df) %in% c("imputed")]

#calculate the actual concentration for each dilution of each specimen
for (i in 1:nrow(df)){
    df$actual_conc[i] = actual_conc(df$orig_2mcf_conc[i], df$dilution[i])
} 

#calculate conc in inoculated blood (1:20 dilution during preparation)
df$inoc_conc = df$actual_conc / 20

#calculate the mean and sd of inoc conc for each dilution
round(tapply(df$inoc_conc, df$dilution, mean), 2)
round(tapply(df$inoc_conc, df$dilution, sd), 2)

#create subsets of inoc conc for comparisons of time to detection
subset1 = df[df$dilution %in% c("500K"), ]
subset2 = df[df$dilution %in% c("10K"), ]
subset3 = df[df$dilution %in% c("100"), ]
subset4 = df[df$dilution %in% c("1"), ]

# actual meand and SD of yeast count in inoculated blood for the most diluted specimens
round(tapply(subset1$inoc_conc, subset1$org_name, mean), 2)
round(tapply(subset1$inoc_conc, subset1$org_name, sd), 2)

#identify trials available for time to detection comparison
ttd = df[!df$hours_diff %in% c(NA), ]
    
#create subsets of inoc conc for comparisons of time to detection
subset1 = ttd[ttd$dilution %in% c("500K"), ]
subset2 = ttd[ttd$dilution %in% c("10K"), ]
subset3 = ttd[ttd$dilution %in% c("100"), ]
subset4 = ttd[ttd$dilution %in% c("1"), ]

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
