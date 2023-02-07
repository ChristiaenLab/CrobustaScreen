library(optparse)

parser <- OptionParser()
parser <- add_option(parser, '--params', action = 'store',default = 'out/params.csv')
parser <- add_option(parser, '--path', action = 'store',default = 'out/raw/2020_12_11/')
opts <- parse_args(parser)

dirs <- list.dirs(opts$path, recursive=F)

cmd <- paste0('runDewakss.R --params ', opts$params, ' --clusts ', dirs, ' --out ',opts$path, dirs)

sapply(cmd, system2, command='Rscript')
