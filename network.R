library(STRINGdb)

ciona <- STRINGdb$new(version='11',species=7719,score_threshold=200, input_directory="")
mouse <- STRINGdb$new(version='11',species=10090,score_threshold=200, input_directory="")

human <- STRINGdb$new(version='11',species=9606,score_threshold=200, input_directory="")


genes <- c('Arhgef8','Gna12/13','Farp','Depdc','Tyrosinase','Eph','Foxf','Gata','Rho',' Dync2h','NDE','EPB4.1','GNA-L/S','Rab5','Paxillin','Rabep','Rock','Col9a1','Ddr' ,'Daam','ADF/Cofilin')

khid <- c('KH.C8.840','KH.L37.10','KH.C8.441','KH.L108.56','KH.C12.469','KH.C1.404','KH.C3.170','KH.L20.1','KH.C2.651','KH.C14.432','KH.C3.451','KH.C9.562','KH.L121.2','KH.L22.59','KH.L108.1','KH.C11.355','KH.C7.324','KH.C8.248','KH.C9.371','KH.C10.209','KH.L37.33')


genes <- data.frame(khid=paste0("KH2013:",khid),name=genes)

genes <- merge(genes,gene.names,by.x=1,by.y=0)

sapply(genes, function(x) grep(x,mouse$get_proteins()[,2],ignore.case=T))


gene.names <- read.delim("gene_name.txt", row.names=1, stringsAsFactors=F)
ensembl <- read.delim('KH-ENS.blast',stringsAsFactors = F,header=F)
ensembl.gene <- data.frame(ensembl_transcript_id=ensembl[,1],KHID=sub('KH2012:(.*).v.*','\\1',ensembl[,2]))
ensembl.gene <- ensembl.gene[!duplicated(ensembl.gene),]
ensembl.gene <- merge(ensembl.gene, 
		      data.frame(KHID=sub('KH2013:','',row.names(gene.names)),KHname=gene.names[,1]))

sel <- sapply(khid,grep,ensembl[,2])

library(biomaRt)

genes <- read.csv('treatments.csv',row.names=1)

mart <- useMart('ensembl','cintestinalis_gene_ensembl')

features <- c(
	"ensembl_gene_id","ensembl_transcript_id","external_gene_name",
	#"refseq_peptide","refseq_peptide_predicted",
	"mmusculus_homolog_ensembl_peptide",
	"hsapiens_homolog_ensembl_peptide"
)

bm <- getBM(attributes=features,mart=mart)

genes <- merge(genes,bm)

lapply(sel, function(x){
	tpt <- ensembl[x,]
	merge(tpt,bm,by.x=1,by.y="ensembl_transcript_id")
})

cionaProtein <- ciona$get_proteins()
cionaRefseq <- cionaProtein$protein_external_id
cionaRefseq <- sub('^7719\\.(.*)\\.[0-9]+$','\\1',cionaRefseq)
cionaPeptide <- select(mart,cionaRefseq,c('ensembl_gene_id','ensembl_transcript_id',"external_gene_name","refseq_peptide","refseq_peptide_predicted"),"refseq_peptide")
cionaPredicted <- select(mart,cionaRefseq,c('ensembl_gene_id','ensembl_transcript_id',"external_gene_name","refseq_peptide","refseq_peptide_predicted"),"refseq_peptide_predicted")
cionaName <- select(mart,cionaProtein$preferred_name,c('ensembl_gene_id','ensembl_transcript_id',"external_gene_name","refseq_peptide","refseq_peptide_predicted"),"external_gene_name")

cionaEnsembl <- rbind(cionaPeptide,cionaPredicted,cionaName)
cionaEnsembl <- cionaEnsembl[!duplicated(cionaEnsembl),]

transcriptToGene <- select(mart,
			   ensembl.gene$ensembl_transcript_id,
			   c("ensembl_gene_id","ensembl_transcript_id","external_gene_name","refseq_peptide","refseq_peptide_predicted"),"ensembl_transcript_id")

ensembl.gene <- merge(ensembl.gene,transcriptToGene,'ensembl_transcript_id')
ensembl.gene <- ensembl.gene[!duplicated(ensembl.gene),]
write.csv(ensembl.gene,'cint_ensembl.csv')

genes <- merge(genes,ensembl.gene,by=1,all.x=T)

genes.hs <- select(mart,genes$ensembl_gene_id,c("ensembl_gene_id","hsapiens_homolog_ensembl_peptide"),'ensembl_gene_id')
genes.hs <- merge(genes,genes.hs,'ensembl_gene_id',all.x=T)

hsProt <- human$get_proteins()
hsProt$protein_external_id <- sub('9606\\.','',hsProt$protein_external_id)
genes.hs <- merge(genes.hs,hsProt,by.x='hsapiens_homolog_ensembl_peptide',by.y='protein_external_id',all.x=T)

genes.mm <- select(mart,genes$ensembl_gene_id,c("ensembl_gene_id","mmusculus_homolog_ensembl_peptide"),'ensembl_gene_id')
genes.mm <- merge(genes,genes.mm,'ensembl_gene_id',all.x=T)

mmProt <- mouse$get_proteins()
mmProt$protein_external_id <- sub('10090\\.','',mmProt$protein_external_id)
genes.mmProt <- merge(genes.mm,mmProt,by.x='mmusculus_homolog_ensembl_peptide',by.y='protein_external_id')

query <- setdiff(genes.mm$UniqueNAME,genes.mmProt$UniqueNAME)

regex <- paste0('^',gsub('\\(\\|','\\(',gsub(
  '([A-Z])|','\\1',
  gsub(
    '/','|',
    gsub(
      "([0-9/]+)","\\(\\1\\)",
      gsub('([0-9])([0-9][A-Za-z])','\\1\\.?\\2',query
    ))))))

regex <- paste0('^',gsub('\\(\\|','\\(',gsub(
  '([A-Z])|','\\1',
  gsub(
    '/','|',
    gsub(
      "([0-9/]+)","\\(\\1\\)",
      gsub('([0-9])([0-9][A-Za-z])','\\1\\.?\\2',genes$UniqueNAME
    ))))),'$')

prots <- lapply(
	setNames(regex,genes$UniqueNAME),
	function(x) mmProt[grep(x,mmProt$preferred_name,T),]
)
prots <- do.call(rbind,mapply(
      function(x,y) {
	      if(nrow(y)>0){
		      y$UniqueNAME <- x
		      merge(genes,y)
	      }
      },names(prots),prots
))

query <- setdiff(genes$ensembl_gene_id,prots$ensembl_gene_id)

genes.mm <- select(mart,query,c("ensembl_gene_id","mmusculus_homolog_ensembl_peptide"),'ensembl_gene_id')
genes.mm <- merge(genes,genes.mm,'ensembl_gene_id')

genes.mmProt <- merge(mmProt,genes.mm,by.y='mmusculus_homolog_ensembl_peptide',by.x='protein_external_id')

prots <- prots[,c(6,4,1:3,5,7:9)]
names(prots) <- names(genes.mmProt)
prots <- rbind(prots,genes.mmProt)

prots <- prots[,c(6,4,1:3,5,7:9)]
names(prots) <- names(genes.mmProt)
prots <- rbind(prots,genes.mmProt)
prots <- prots[!duplicated(prots[,c(1,3:9)]),]

library(igraph)
mmgraph <- mouse$get_graph()
E(mmgraph)

pdf('mmusculus_homologs.pdf')
mouse$plot_network(prots$protein_external_id)
dev.off()
