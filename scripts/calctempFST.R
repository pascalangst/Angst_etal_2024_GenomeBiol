# caculate pairwiseFST between timepoints in slim simulations

# the calculation is based on the slim function calcFST

slim_muts_0 <- read.table("slim_muts_0.tmp", quote="\"", comment.char="")
slim_muts_1 <- read.table("slim_muts_1.tmp", quote="\"", comment.char="")
slim_muts_2 <- read.table("slim_muts_2.tmp", quote="\"", comment.char="")
slim_freq_0 <- read.table("slim_freq_0.tmp", quote="\"", comment.char="")
slim_freq_1 <- read.table("slim_freq_1.tmp", quote="\"", comment.char="")
slim_freq_2 <- read.table("slim_freq_2.tmp", quote="\"", comment.char="")

df_0 <- data.frame(mut0 = t(slim_muts_0), frq0 = t(slim_freq_0))
df_1 <- data.frame(mut1 = t(slim_muts_1), frq1 = t(slim_freq_1))
df_2 <- data.frame(mut2 = t(slim_muts_2), frq2 = t(slim_freq_2))

df <- merge(df_1, df_2, by.x = "mut1", by.y = "mut2", all = T)
df[is.na(df)] <- 0

mean_p = (df$frq1 + df$frq2) / 2.0
H_t = 2.0 * mean_p * (1.0 - mean_p)
H_s = df$frq1 * (1.0 - df$frq1) + df$frq2 * (1.0 - df$frq2)
fst = 1.0 - mean(H_s) / mean(H_t)

write(fst, "fst_1.tmp")

df <- merge(df_0, df_2, by.x = "mut0", by.y = "mut2", all = T)
df[is.na(df)] <- 0

mean_p = (df$frq0 + df$frq2) / 2.0
H_t = 2.0 * mean_p * (1.0 - mean_p)
H_s = df$frq0 * (1.0 - df$frq0) + df$frq2 * (1.0 - df$frq2)
fst = 1.0 - mean(H_s) / mean(H_t)

write(fst, "fst_0.tmp")

system("mv slim_muts_2.tmp slim_muts_1.tmp;
       mv slim_freq_2.tmp slim_freq_1.tmp")
