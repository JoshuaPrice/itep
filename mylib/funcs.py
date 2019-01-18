### Functions

def mm9to10(fileDir, bedIn, bedOut):
	# Lifts mm9 records over to mm10 records
	from subprocess import call

	chainfile = '/data2/josh/genome/mm9ToMm10.over.chain.gz'
	bed_name_mm9 = fileDir + bedIn
	 #bed_name_mm9 = '/data2/josh/ep/focs_promoters_mm9_cleaned.bed6'
	bed_name_mm10 = fileDir + bedOut

	call(["CrossMap.py", "bed", chainfile, bed_name_mm9, bed_name_mm10])


# couldn't call this from freq folder (not sure why) - 11/20/18
def removeQuotes(fileIn, fileOut):
	# Removes all quotations from a csv file
	with open(fileIn, 'r') as f, open(fileOut, 'w') as fo:
	    for line in f:
	        fo.write(line.replace('"', '').replace("'", ""))
