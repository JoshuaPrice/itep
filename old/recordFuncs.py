from Bio import SeqIO, Seq, SeqRecord, Alphabet
import csv

def getRecordNum(filepath):
	# get the number of records in the fasta file
	num_records = 0
	for seq_record in SeqIO.parse(filepath, "fasta"):
		num_records += 1
	return num_records

def splitRecords(filepath):
	# split the fasta file into many fasta files, each with one record
	xsome_files = []
	for seq_record in SeqIO.parse(filepath, "fasta"):
		record_name = seq_record.id
		new_filepath = "data/xsomes/xsome_" + str(record_name) + ".fa"
		xsome_files.append(new_filepath)
		SeqIO.write(seq_record, new_filepath, "fasta")

	with open('data/xsomes/xsome_files.csv', 'wb') as names_file:
		names_writer = csv.writer(names_file, delimiter=',')
		names_writer.writerow(xsome_files)


if __name__ == '__main__':
	splitRecords("data/Mus_musculus.GRCm38.dna.primary_assembly.fa")