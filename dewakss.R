runDewakss <- function(
	params='out/params.csv',
	groups='out/groups.csv',
	out='out/raw'
	#res=1
	){
	cmd <- paste(
		'python3 runDewakss.py --params', params,
		'--groups', groups,
		'--out', out
	#	'--res', as.character(res)
	)
	system(cmd)
}
