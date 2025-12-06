#######################################################################################################
## Author	: M.L. Dubbelaar, adapted by A.Alsema. contact info a.m.alsema@umcg.nl
## Date		: 17-02-2018
## Name 	: UMItools_GEO.py
## Purpose	: Generates fastqc, various alignment log files, and countfiles of barcoded-smartseq2 single cell RNA sequencing data.
## Data structure: read 1 has UMI+barcode, read 2 cDNA of an adapted, barcoded-smartseq2 single cell RNA seq protocol described in PMID: 33192286
#######################################################################################################
import os
import sys
import glob
import commands
import argparse

knownBarcodeInfo = []
foundBarcodes = []

featureCounts = "/home/user/subread-1.6.0-Linux-x86_64/bin/featureCounts"
picard = "/home/user/picard/build/libs/picard.jar"
hisat2 = "/home/user/hisat2-2.1.0/hisat2"
humangtf = "/home/user/Homo_sapiens.GRCh38.91.gtf"
humanGenome = "/home/user/Human/GRCh38.91"
splicetxt = "/home/user/Homo_sapiens.GRCh38.91_spliceSites.txt"
threads = "5"

# declare arguments
dirToSCfiles = ""
knownBarcodeFile = ""
numberOfBarcodes = ""
pattern = ""
outputdir = "" 
summary = ""

# define arguments
def arguments ():
	'''
		This functions defines the different parameters used for this function and gives a help output if the parameters are not given.
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--fastqDir', help='Directory to the single cell fastq.gz files')
	parser.add_argument('-b', '--barcode', help="Give a text file with the known barcodes")
	parser.add_argument('-p', '--pattern', help="Define a UMI (C) and/or cell barcode (N) pattern")
	parser.add_argument('-o', '--outputDir', help='Pathways to output directory')
	parser.add_argument('-s', '--summaryDir', help='Pathways to countfile + multiQC summary directory')	

	args = parser.parse_args()

	global dirToSCfiles
	global knownBarcodeFile
	global pattern
	global outputdir
	global summary
	
	# if (not args.barcode or not args.fastqDir or not args.outputDir or not args.pattern or not args.splice):
	if (not args.barcode or not args.fastqDir or not args.outputDir or not args.pattern):
		print("Parameters: \n" +
		"-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n" + 
		"	-b/--barcode 					|		[Required] Path to a txt files that contains the known barcodes.\n" +
		"	-d/--fastqDir					|		[Required] Path to directory with all the fastq.gz files.\n" +
		"	-o/--outputDir 					|		[Required] Path to write the output to.\n" +
		"	-p/--pattern 					|		[Required] Define the pattern of the UMI (N) and/or barcode (C).\n"+
		"	-s/--summaryDir					|		[Required] Path to directory were countfile output and multiqc ends.\n"+
		"-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n\n" +

		"Use example:\n" +
		"-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n" +
		"	python pipeline_UMItools.py -b /home/user/readinCBC.csv -d /home/user/fastqDir/ -o /home/user/outputDir -p NNNNNNNCCCCCCCCCC - s /home/user/countfiles/ \n" +
		"-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n" +
		"	!IMPORTANT!\n" + 
		"	This files expects that there are 2 files: R1 that consists of the UMI and/or cell barcode, and R2 containing the sequences.\n" +
		"	Without these files the tool will not work, modification of the umi_tools extract method is necessary!\n" +
		"-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
		sys.exit()
	else:
		knownBarcodeFile = args.barcode
		dirToSCfiles = args.fastqDir
		outputdir = args.outputDir
		pattern = args.pattern
		summary = args.summaryDir
		
def generatefastqc():
	direct = os.getcwd()
	os.chdir(direct)
	myFileList = glob.glob('*.fastq.gz')
	print(myFileList)
	print("--------The number of files to analyze --------")
	print(len(myFileList))
	Nr = 0
	for myfile in myFileList:
		print(myfile)
		os.system("fastqc -o " + outputdir + " " + myfile)
		Nr+=1
		print("--------The number of files qc'ed " + str(Nr) + " --------")


def obtainKnownBarcodes():
	'''
		This functions obtaines the known in knownBarcodeFile from the defined file and the amount of barcodes is calculated and saved in numberOfBarcodes.
	'''
	global numberOfBarcodes
	# The known barcodes are obtained from the defined file and the amount of barcodes is calculated
	for knownBarcode in open(knownBarcodeFile, "r"):
		knownBarcodeInfo.append(knownBarcode.split()[0])
	numberOfBarcodes = len(knownBarcodeInfo)

def filterUnknownBarcodes(nameOfFile, nameOfAdjusted):
	'''
		This functions filters the unknown barcodes from the generated whitelists, saving the known barcodes into a new file (..._adj.txt).
	'''
	whitelistFile = open(nameOfFile, "r")
	lines = whitelistFile.readlines()
	whitelistFile.close()

	f = open(nameOfAdjusted,"w")
	for line in lines:
		barcode = line.strip().split("\t")[0]
		if [s for s in knownBarcodeInfo if barcode in s]:
			f.write(line.strip() + "\n")
	f.close()
	os.system("chmod -R oug+rwx " + outputdir)

def mappingMetrics():
	os.chdir(outputdir)
	#grep a list of sorted, indexed BAM files from directory
	myFileList = glob.glob('*sort_R2.bam')
	print("--------input files-----------")
	print(myFileList)
	print("--------The number of files to collect QC --------")
	print(len(myFileList))
	#use CollectRNASeqMetrics from Picard to calculate % mapped to exon, intron, UTR, and intergenic regions of the reference genome
	for files in myFileList:
		filesName = files.split("/")[-1]
		filesName = filesName.split("_") 
		filesName = "_".join(filesName[:5])  # MUST BE ADAPTED IF YOU USE A NON-STANDARD filename. Standard filename looks like iB6408-BA7_CTR_37-3_S10_sort_R2.bam
		print("\n")
		print("---------- " + filesName + " ----------")
		print("-------- calculate mapping statistics --------")
		os.system("java -jar " + picard + " CollectRnaSeqMetrics I=" + outputdir + "/" + files + " O=" + filesName + "_RNA_Metrics.txt REF_FLAT=" + outputdir + "/" + "refFlat.txt STRAND=NONE")

def main():
	'''
		All of the necessary steps are defined into the main.
			- Creating the whitelist
			- Filtering the unknown barcodes
			- Extracting the defined pattern
			- Alignment
			- FeatureCounts
			- Sorting and indexing
			- Creating the count file based on barcode
	'''
	arguments()
	generatefastqc()
	obtainKnownBarcodes()
	for files in glob.glob(dirToSCfiles + "*R1_001.fastq.gz"):
		filesName = files.split("/")[-1]
		filesName = filesName.split("_") 
		filesName = "_".join(filesName[:4])  # MUST BE ADAPTED IF YOU USE A NON-STANDARD filename. Standard filename looks like iB6408-BA7_CTR_37-3_S10_Rx_001.fastq.gz
		print("\n")
		print("---------- " + filesName + " ----------")
		print("### Creating Whitelist ###\n")
		# https://github.com/CGATOxford/UMI-tools/blob/master/doc/Single_cell_tutorial.md
		os.system("umi_tools whitelist -p " + pattern + " --extract-method=string --set-cell-number=" + str(numberOfBarcodes) + " -I " + files +
		" --log2stderr --plot-prefix=" +  outputdir + filesName + " > " + outputdir + "whitelist_" + filesName + ".txt")
		os.system("chmod -R  oug+rwx " + outputdir + "whitelist_" + filesName + ".txt")

		# Removes the unknown barcodes in the whitelist and writes this to an _adj file.
		filterUnknownBarcodes(outputdir + "whitelist_" + filesName + ".txt", outputdir + "whitelist_" + filesName + "_adj.txt")
		os.system("chmod -R oug+rwx " + outputdir)

		print("")
		print("### Extract ###\n")
		# Extract cellbarcode and put it in read1 header
		# https://github.com/CGATOxford/UMI-tools/blob/master/doc/Single_cell_tutorial.md
		#parallel -j 6 'umi_tools extract -I results/{} --read2-in=results/{}.MATEPAIR --bc-pattern=NNNNNNNNNN --log=processed.log --stdout=results_UMI/{}.read1.fastq --read2-out=results_UMI/{}.read2.fastq' ::: path/to/your/files/*.fileending
		os.system("umi_tools extract --whitelist " + outputdir + "whitelist_" + filesName + "_adj.txt --extract-method=string --error-correct-cell --bc-pattern=" + pattern + " --filter-cell-barcode -I " +
		files + " --stdout " + outputdir + filesName + "_R1.fastq.gz --read2-in " + dirToSCfiles + filesName + "_R2_001.fastq.gz --read2-out=" + outputdir + filesName + "_R2.fastq.gz -L " + outputdir + "ExtractLog.txt")
		os.system("chmod -R oug+rwx " + outputdir)

		print("")
		print("### Align ###\n")
		# Aligning with splice sites
		os.system(hisat2 + " --threads " + str(threads) +  " -x " + humanGenome + " -U " + outputdir + filesName + "_R2.fastq.gz --known-splicesite-infile " +  splicetxt + " 2>> " + outputdir + filesName + ".log | samtools view -Sbo " + outputdir + filesName + "_R2.bam")
		print(hisat2 + " --threads " + str(threads) + " -x " + humanGenome + " -U " + outputdir + filesName + "_R2.fastq --known-splicesite-infile " +  splicetxt + " 2>> " + outputdir + filesName + ".log | samtools view -Sbo " + outputdir + filesName + "_R2.bam")
		os.system("chmod -R oug+rwx " + outputdir)
		print("### Sort ###\n")
		os.system("samtools sort " + outputdir + filesName + "_R2.bam -o " + outputdir + filesName + "_sort_R2.bam -@ " + str(threads))
		os.system("chmod -R oug+rwx " + outputdir)
		print("### Index ###\n")
		os.system("samtools index -@ 5 " + outputdir + filesName + "_sort_R2.bam")
		os.system("chmod -R oug+rwx " + outputdir)
		
		
		print("")
		print("### Obtain featureCount bam ###\n") # -- primary means select only uniquely mapped reads
		os.system(featureCounts + " --primary -a " + humangtf + " -o " + outputdir + filesName + "_count.txt -R BAM " +  outputdir + filesName + "_sort_R2.bam -T 5")
		
		os.system("chmod -R oug+rwx " + outputdir)

		print("")
		print("### Generate Count Files ###\n")
		os.system("samtools sort " + outputdir + filesName + "_sort_R2.bam.featureCounts.bam -o " + outputdir + filesName + "_sortedFC_R2.bam -@ " + str(threads))
		
		os.system("chmod -R oug+rwx " + outputdir)
		os.system("samtools index " + outputdir + filesName + "_sortedFC_R2.bam -@ " + str(threads))
		
		os.system("umi_tools count --per-gene --per-cell --gene-tag XT --wide-format-cell-counts -L " + outputdir + filesName + "UMIcount_log.txt -I " + outputdir + filesName + "_sortedFC_R2.bam -S " + summary + filesName + "_completeCounts.txt")
		os.system("chmod -R oug+rwx " + outputdir)
		print("### Finished with CountFiles ###\n")


if __name__ == "__main__":
	main()
	mappingMetrics() # to do: move this to the main above.
	print("finished mapping and QC")
	os.system("sudo multiqc " + outputdir + " -o " + summary)
