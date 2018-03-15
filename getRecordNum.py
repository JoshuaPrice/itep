from Bio import SeqIO

def getRecordNum(filepath):
	num_records = 0
	for seq_record in SeqIO.parse(filepath, "fasta"):
		num_records += 1
	return num_records

def splitXsomeRecords(filepath):
	for seq_record in SeqIo.parse()