require("org.Hs.eg.db")
require("GO.db")
require("topGO")
require("plyr")
require("readr")
require("dplyr")
load("Rdata_output/gencode_universe.Rdata")
# acceptable gene ID types for TopGo
#c("entrez", "genbank", "alias", "ensembl", "symbol", "genename", "unigene")
{
TopGO_over <- function(observedIDs, universeIDs, outdir, proj, gname='ensembl') 
{  
  
  geneFile = paste0(outdir, '/',proj, "_pathway_genes.txt")
  resFile = paste0(outdir, '/',proj, "_pathway.txt")
  
  geneList <- factor(as.integer((universeIDs %in% observedIDs)))
  names(geneList) <- universeIDs
 
  # Create topGO object
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList, 
                nodeSize = 10,
                annotation = annFUN.org, 
                mapping = "org.Hs.eg.db",
                ID= gname)
  
  # Run elimination test 
  resultTopGO.weight.fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher" )
  
  allRes <- GenTable(GOdata, 
                     classicfisher = resultTopGO.weight.fisher, 
                     numChar=120,
                     orderBy = "classicfisher", 
                     topNodes = length(usedGO(object = GOdata)))
  
  allRes2 <- allRes %>% 
  mutate(classicfisher = ifelse(classicfisher =="< 1e-30", 1e-30, classicfisher)) %>% 
  mutate(q.value=p.adjust(classicfisher, method="BH"))
    

  write_delim(allRes2, resFile, delim = "\t")
  
  #showSigOfNodes(GOdata, score(resultTopGO.weight.fisher), firstSigNodes = 5, useInfo = 'all')
  # printGraph(GOdata, resultTopGO.weight.fisher, firstSigNodes = 5, 
  #            fn.proj = pdfFile, useInfo = "all", pdfSW = TRUE)
  # 
  
  # Pull out genes of interest that belong to significantly enriched GO categories
  allGO <- genesInTerm(GOdata)
  sigGO_ids <- allRes[1:100,]$GO.ID
  sigGO <- allGO[sigGO_ids]
  sigInTerm <- lapply(sigGO, function(x) x[x %in% observedIDs] )
  sigInTerm <- plyr::ldply(sigInTerm, data.frame)
  
  names(sigInTerm) <- c("GO_id", "Gene")
  sigInTerm <- sigInTerm %>%
    mutate(Gene = as.character(Gene)) %>% 
    left_join(.,GeneTable) 
  
  write.table(sigInTerm, geneFile, sep = "\t", quote = F, row.names = F, col.names = T )
  #return(list(GOdata, resultTopGO.weight.fisher, resultTopGO.classic.fisher,allRes))
}


left_join0 <- function(x, y, fill = 0L,...){
  z <- left_join(x, y,...)
  tmp <- setdiff(names(z), names(x))
  z <- replace_na(z, setNames(as.list(rep(fill, length(tmp))), tmp))
  z
}

setGeneric("GOFisherUnder", function(object) 
  standardGeneric("GOFisherUnder") )

setMethod("GOFisherUnder", "classicCount", function(object){
            contMat <- contTable(object)
            if(all(contMat == 0))
              p.value <- 1
            else
              p.value <- fisher.test(contMat, alternative = "less")$p.value
            return(p.value)
          }
)



TopGO_under <- function(observedIDs, universeIDs, outdir, proj, gname='ensembl') 
{  
  
  geneFile = paste0(outdir, '/',proj, "_pathway_genes.txt")
  resFile = paste0(outdir, '/',proj, "_pathway.txt")
  
  geneList <- factor(as.integer((universeIDs %in% observedIDs)))
  names(geneList) <- universeIDs
  
  # Create topGO object
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList, 
                nodeSize = 10,
                annotation = annFUN.org, 
                mapping = "org.Hs.eg.db",
                ID= gname)
  #classicCount
  # Run elimination test and write the results file
  # test.stat <- new("weightCount", testStatistic = GOFisherUnder, 
  #                  name="Weight with Fisher's exact test for under-representation")
  test.stat <- new("classicCount", testStatistic = GOFisherUnder, 
                   name="Classic with Fisher's exact test for under-representation")
  resFisher.weight <- getSigGroups(GOdata, test.stat)
  

  allRes <- GenTable(GOdata, 
                     classicfisher = resFisher.weight, 
                     topNodes = length(usedGO(object = GOdata)) , 
                     numChar=120, 
                     orderBy = "classicfisher")
  
  allRes2 <- allRes %>% 
    mutate(classicfisher = ifelse(classicfisher =="< 1e-30", 1e-30, classicfisher)) %>% 
    mutate(q.value=p.adjust(classicfisher, method="BH")) %>% 
    mutate(fdr=p.adjust(classicfisher, method="fdr"))  
  
  write_delim(allRes2, resFile, delim = "\t")
  
  # Pull out genes of interest that belong to significantly enriched GO categories
  allGO <- genesInTerm(GOdata)
  sigGO_ids <- allRes[1:100,]$GO.ID
  sigGO <- allGO[sigGO_ids]
  sigInTerm <- lapply(sigGO, function(x) x[x %in% observedIDs] )
  sigInTerm <- plyr::ldply(sigInTerm, data.frame)
  
  names(sigInTerm) <- c("GO_id", "Gene")
  sigInTerm <- sigInTerm %>%
    mutate(Gene = as.character(Gene)) %>% 
    left_join(.,GeneTable)
  
  
  write.table(sigInTerm, geneFile, sep = "\t", quote = F, row.names = F, col.names = T )
}
}