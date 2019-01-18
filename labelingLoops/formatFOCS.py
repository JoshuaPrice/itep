import pandas as pd
from pybedtools.bedtool import BedTool
from subprocess import call

def getFocsPromotersMm10():
	# Returns dataframe in Bed6 format
	return pd.read_csv("/data2/josh/ep/focs_promoters_mm10.bed6",
						sep='\t', 
						header=None, 
						names=['chr','pos','pos2','gene','score','strand']
						)

def formatFocsPromoters(): 
	# note: not all sym_id imported successfully because some promoters had multiple corresponding genes
	raw_p_mm9 = pd.read_csv('/data2/josh/ep/focs_promoters_mm9.csv', 
	                    usecols=['name','score','thick.start','thick.width','sym_id']
	                   )

	# Clean mm9 dataframe to be in bed format
	raw_p_mm9['chr'] = raw_p_mm9['name']
	raw_p_mm9[['chr','raw_range']] = raw_p_mm9['name'].str.split(':',expand=True)
	# raw_p_mm9['chr'] = raw_p_mm9['chr'].str[3:]
	raw_p_mm9['pos2'] = raw_p_mm9['thick.start']
	raw_p_mm9['strand'] = '.'
	p_mm9 = raw_p_mm9.drop(columns=['name','raw_range','thick.width'])
	p_mm9.rename({'thick.start': 'pos', 'sym_id': 'gene'},
					 axis='columns', 
					 inplace=True
					 )
	p_bed_cols = ['chr','pos','pos2','gene','score','strand']
	p_mm9 = p_mm9[p_bed_cols]

	# Create bed file
	bed_name_mm9 = '/data2/josh/ep/focs_promoters_mm9_cleaned.bed6'
	p_mm9.to_csv(bed_name_mm9,sep='\t',index=False, header=False)
	
	# lift over to mm10
	chainfile = "/data2/josh/genome/mm9ToMm10.over.chain.gz"
	bed_name_mm10 = "/data2/josh/ep/focs_promoters_mm10.bed6"
	call(["CrossMap.py", "bed", chainfile, bed_name_mm9, bed_name_mm10])