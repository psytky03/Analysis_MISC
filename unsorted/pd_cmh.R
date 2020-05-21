require("dplyr")
require("coin") # use cmh_test function from coin package 
 
# Implementation of the method to test the population differnt genes proposed in 
# Ahn, Eunyong, and Taesung Park. "Analysis of population-specific pharmacogenomic variants using next-generation sequencing data." 
# Scientific reports 7.1 (2017): 1-11.

pd_cmh <- function(df){
  x <- df %>% 
    select(jm, jM, em, eM) #allele count of m, M for japanese (j) and European (e)
  k <- nrow(x)
  
  # make a 2*2*k table by flatening the data frame 
  x1 <- as.table(array(c(t(x)),dim=c(2,2,k)))
  print(x1)
  y = cmh_test(x1)
  
  # suppose 2 * 2 * k tables 
  # a b
  # c d
  # Cochran-Mantel-Haenszel OR: sum of (a*d/n) / sum of (b*c/n) of each 1:k 
  oeupper <- function(x){x[1,1]*x[2,2]/(x[1,1] + x[1,2] + x[2,1] + x[2,2])}
  oelower <- function(x){x[1,2]*x[2,1]/(x[1,1] + x[1,2] + x[2,1] + x[2,2])}                 
  or_baseline_e <- sum(unlist(apply(x1, 3, oeupper)))/sum(unlist(apply(x1, 3, oelower)))
  
  y1 = data.frame(p = pvalue(y), chisq = statistic(y), nsnp = k, or = or_baseline_e)
  
  return(y1)
}

# Test with gnomAD data
gnomad_essential <- gnomad_LoF_ALL_pickup_bi %>% 
  select(ID, Gene, SYMBOL, AC_nfe, AN_nfe, AF_nfe, AC_eas, AN_eas, AF_eas, LoF) %>% 
  filter(AC_nfe > 0 | AF_eas >0 ) %>% 
  mutate(eM = AN_nfe - AC_nfe ) %>% 
  rename(em = AC_nfe) %>% 
  mutate(jM = AN_eas - AC_eas) %>% 
  rename(jm = AC_eas) %>% 
  mutate(em = if_else(em == 0, 1, em)) %>% #add 1 for non-polymorphic site 
  mutate(jm = if_else(jm == 0, 1, jm)) %>% #add 1 for non-polymorphic site 
  unique() %>% 
  filter(AF_nfe > 0.005 | AF_eas > 0.005) %>% 
  filter(AF_nfe < 0.15 & AF_eas < 0.15) %>%  
  filter(LoF == 'HC')

gnomad_pd <- gnomad_essential %>% 
  group_by(Gene, SYMBOL) %>% 
  group_modify(~pd_cmh(.x)) %>% 
  ungroup()

gnomad_pd_jp <- gnomad_pd%>% 
  filter(or > 2) %>% 
  filter(p < 1e-5) %>% 
  filter(!str_detect(SYMBOL,"^OR") ) %>% 
  arrange(-chisq)



