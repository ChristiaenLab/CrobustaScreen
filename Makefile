DATE=$(date +%Y-%m-%d)

out/raw/$DATE: out/params.csv
	python3 runDewakss.py

out/z/$DATE: out/z.csv
	python3 runDewakss.py --params out/z_dat.csv --out out/z

out/phenoSub/$DATE: out/phenoSub.csv
	python3 runDewakss.py --params out/phenoSubParam.csv --groups out/phenoSub.csv --out out/phenoSub

out/phenoSubZ/$DATE: out/phenoSubZ.csv
	python3 runDewakss.py --params out/phenoSubZ.csv --groups out/phenoSub.csv --out out/phenoSubZ

out/params.csv:
	Rscript readEmbryos.R
