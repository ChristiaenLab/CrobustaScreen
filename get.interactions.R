ciona <- STRINGdb$new(version='11',species=7719,score_threshold=200, input_directory="")
mouse <- STRINGdb$new(version='11',species=10090,score_threshold=200, input_directory="")
human <- STRINGdb$new(version='11',species=9606,score_threshold=200, input_directory="")

crgraph <- ciona$get_graph()
mmgraph <- mouse$get_graph()
hsgraph <- human$get_graph()

write.csv(do.call(rbind,E(crgraph)),'STRINGdb/cintestinalis.graph.csv')

cint <- read.csv("STRINGdb/cintestinalis.csv",row.names=1)
mmus <- read.csv("STRINGdb/mmusculus.csv",row.names=1)
hsap <- read.csv("STRINGdb/hsapiens.csv",row.names=1)

mm.ext <- read.csv("STRINGdb/mmusculus.ext.csv",row.names=1)
hs.ext <- read.csv("STRINGdb/hsapiens.ext.csv",row.names=1)

treat <- read.csv("treatments.csv",row.names=1)

cr.treat <- merge(treat,cint,by='ensembl_gene_id')
cr.int <- ciona$get_interactions(cr.treat$protein_external_id)

merge.sdb <- function(sdb,treat,prots,ext,out){
	ext.treat <- merge(treat,ext)

	treat <- merge(treat[treat$khid%in%setdiff(treat$khid,ext.treat$khid),],
		       prots,by='ensembl_gene_id')
	treat <- merge(treat,ext.treat,all.x=T,all.y=T)
	int <- sdb$get_interactions(treat$protein_external_id)
	int <- int[!duplicated(int),]

	kh <- treat[,c('protein_external_id','khid','name')]
	kh <- kh[!duplicated(kh),]
	int <- merge(int,kh,by.x='from',by.y='protein_external_id')
	int <- merge(int,kh,by.x='to',by.y='protein_external_id')
	int.kh <- int[,c('khid.x','khid.y','name.x','name.y')]
	int.kh <- int.kh[!duplicated(int.kh),]
	write.csv(int.kh,out)
	return(int.kh)
}

mmus.int <- merge.sdb(mouse,treat,mmus,mm.ext,'mmus.interactions.csv')

hsap.int <- merge.sdb(human,treat,hsap,hs.ext,'hsap.interactions.csv')

int <- rbind(mmus.int,hsap.int)
int <- int[!duplicated(int),]
write.csv(int,'interactions.csv')

