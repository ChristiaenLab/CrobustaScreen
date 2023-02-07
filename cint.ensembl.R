library(biomaRt)

dir.create('ensembl')

gene.names <- read.delim("gene_name.txt", row.names=1, stringsAsFactors=F)
ensembl <- read.delim('KH-ENS.blast',stringsAsFactors = F,header=F)
ensembl.gene <- data.frame(ensembl_transcript_id=ensembl[,1],KHID=sub('KH2012:(.*).v.*','\\1',ensembl[,2]))
ensembl.gene <- ensembl.gene[!duplicated(ensembl.gene),]
ensembl.gene <- merge(ensembl.gene, 
		      data.frame(KHID=sub('KH2013:','',row.names(gene.names)),KHname=gene.names[,1]))

write.csv(ensembl.gene,'ensembl/khid.csv')

mart <- useMart('ensembl','cintestinalis_gene_ensembl')

features <- c(
	"ensembl_gene_id","ensembl_transcript_id","external_gene_name",
	"refseq_peptide","refseq_peptide_predicted"
)
ortho <- c(
	"ensembl_gene_id","ensembl_transcript_id","external_gene_name",
	"mmusculus_homolog_ensembl_peptide",
	"hsapiens_homolog_ensembl_peptide"
)

bmfeat <- getBM(attributes=features,mart=mart)
bmortho <- getBM(attributes=ortho,mart=mart)

write.csv(bmfeat,'ensembl/features.csv')
write.csv(bmortho,'ensembl/orthologs.csv')

transcriptToGene <- select(mart,
			   ensembl.gene$ensembl_transcript_id,
			   c("ensembl_gene_id","ensembl_transcript_id","external_gene_name","refseq_peptide","refseq_peptide_predicted"),"ensembl_transcript_id")

ensembl.gene <- merge(ensembl.gene,transcriptToGene,'ensembl_transcript_id')
ensembl.gene <- ensembl.gene[!duplicated(ensembl.gene),]
write.csv(ensembl.gene,'cint_ensembl.csv')

