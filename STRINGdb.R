library(STRINGdb)

ciona <- STRINGdb$new(version='11',species=7719,score_threshold=200, input_directory="")
mouse <- STRINGdb$new(version='11',species=10090,score_threshold=200, input_directory="")
human <- STRINGdb$new(version='11',species=9606,score_threshold=200, input_directory="")

ensembl <- read.csv('ensembl/features.csv',row.names=1)

ortho <- read.csv('ensembl/orthologs.csv',row.names=1)
ortho <- ortho[apply(ortho[,4:5]!='',1,all),]
genes.ortho <- merge(ensembl,ortho,all.y=T)

write.csv(genes.ortho,'ensembl/crobusta_orthologs.csv')

cionaProtein <- ciona$get_proteins()
cionaRefseq <- cionaProtein$protein_external_id
cionaProtein$refseq_peptide <- sub('^7719\\.(.*)\\.[0-9]+$','\\1',cionaRefseq)

cionaPeptide <- merge(ensembl,cionaProtein,by='refseq_peptide')
cionaPred <- merge(ensembl,cionaProtein,by.x='refseq_peptide_predicted',by.y='refseq_peptide')

cionaName <- merge(ensembl,cionaProtein,by.x=c('refseq_peptide','external_gene_name'),by.y=c('refseq_peptide','preferred_name'),all.x=T,all.y=T)
cionaName <- cionaName[!apply(cionaName,1,anyNA),]
cionaName$preferred_name <- cionaName$external_gene_name

prots <- rbind(cionaPeptide,cionaPred[,names(cionaPeptide)],cionaName[,names(cionaPeptide)])
prots <- prots[!duplicated(prots),]

dir.create('STRINGdb')
write.csv(prots,'STRINGdb/cintestinalis.csv')

hsProt <- human$get_proteins()
hsProt$refseq_peptide <- sub('9606\\.','',hsProt$protein_external_id)
genes.hs <- merge(genes.ortho,hsProt,by.x='hsapiens_homolog_ensembl_peptide',by.y='refseq_peptide')

write.csv(genes.hs,'STRINGdb/hsapiens.csv')

mmProt <- mouse$get_proteins()
mmProt$refseq_peptide <- sub('10090\\.','',mmProt$protein_external_id)
genes.mm <- merge(genes.ortho,mmProt,by.x='mmusculus_homolog_ensembl_peptide',by.y='refseq_peptide')

write.csv(genes.mm,'STRINGdb/mmusculus.csv')

query <- setdiff(mmProt$preferred_name,genes.mm$external_gene_name)

gene.names <- read.delim("cionaGeneIDs/KH2013-UniqueNAME.txt")

gene.sub <- function(sdb,sdbq,genes,genesq){
	

match.names <- function(sdb,genes){
	regex <- paste0('^',gsub('\\(\\|','\\(',gsub(
	  '([A-Z])|','\\1',
	  gsub(
	    '/','|',
	    gsub(
	      "([0-9/]+)","\\(\\1\\)",
	      gsub('([0-9])([0-9][A-Za-z])','\\1\\.?\\2',genes
	    ))))),'$')

	out <- sapply(regex,grepl,sdb,T)
	colnames(out) <- genes
	row.names(out) <- sdb
	return(out)
}

ortho.ext <- function(sdb,genes,out){
	ext <- match.names(sdb$preferred_name,genes$UniqueNAME)

	res <- lapply(1:nrow(genes),function(x) {
			 tmp=sdb[ext[,x],]
			 if(nrow(tmp)>0) 
				 cbind(tmp,genes[x,,drop=F])
	})
	res <- do.call(rbind,res)
	write.csv(res,paste0("STRINGdb/",out,'.ext.csv'))
	return(res)
}

hs.ext <- ortho.ext(hsProt,gene.names,'hsapiens')
mm.ext <- ortho.ext(mmProt,gene.names,'mmusculus')

