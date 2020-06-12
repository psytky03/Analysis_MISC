install.packages("Rmpfr")
library("Rmpfr")
# X.CHROM	POS	ID	REF	ALT	A1	A1_FREQ	MACH_R2	TEST	OBS_CT	BETA	SE	T_STAT
# 11	64361219	11:64361219:G:A	G	A	A	0.0225663	0.988891	ADD	104175	-1.24032	0.0140482	-88.2903

t.value <- abs(-88.2903)
extreme.P <- 2*Rmpfr::pnorm(mpfr(t.value, precBits=10), lower.tail=FALSE, log.p = FALSE)
extreme.P #7.6657e-1696

