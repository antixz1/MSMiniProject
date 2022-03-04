#Python Wrapper
import os
import subprocess
import csv


#SRA Accession to pull from NCBI
accession = 'SRR8185310'

#Make results directory and change working directory to it
os.chdir(os.path.expanduser("~"))
os.system('mkdir results')

#1. Set path to sratoolkit bin and prefetch desired accession in sra format. Use fasterq-dump to create fastq file from the sra file located in results folder.
def sra(accession):
    os.chdir(os.path.expanduser("~/sratoolkit.2.11.2-ubuntu64/bin")) #Change dir to sratoolkit bin
    print('Running sra prefetch for '+accession)
    os.system('./prefetch '+accession+' -O'+os.path.expanduser("~/results")) #Grabs run from ncbi in sra format and places in results dir
    print('Running sra fasterq-dump on '+accession+'.sra')
    os.system('./fasterq-dump '+os.path.expanduser("~/results/"+accession+"/"+accession+".sra")+' -O '+os.path.expanduser("~/results")) #Converts .sra file to .fastq file and places in results dir


sra(accession)


#2. Use SPAdes to assemble the genome using the fastq file generated from sratools fasterq-dump
def spades(fastq):
    os.chdir(os.path.expanduser("~/SPAdes-3.15.4-Linux/bin"))
    os.system('python3 spades.py -k 55,77,99,127 -t 2 --only-assembler -s '+os.path.expanduser("~/results/"+fastq)+' -o '+os.path.expanduser("~/results")) 
    os.chdir(os.path.expanduser("~/results"))
    os.system('echo python3 spades.py -k 55,77,99,127 -t 2 --only-assembler -s '+os.path.expanduser("~/results/"+fastq)+' -o '+os.path.expanduser("~/results")+' > miniproject.log') #Write SPAdes command to miniproject.log
    
spades(accession+'.fastq')



#3-4. Count the number of contigs longer than 1000 bp and count the total assembly length of those contigs containing > 1000 bp
def contig_count(contig_file):
    os.chdir(os.path.expanduser("~/results")) #Change directory to location of SPAdes output
    #Use SeqIO from BioPython to parse fasta file into list 
    with open(contig_file) as f_in:
        fasta = {}
    
        for line in f_in:
            if line.startswith('>'): #If line is a header, add it to fasta dict as a key and each subsequent line as the value until another header line is reached
                header = line.strip() 
                fasta[header] = '' 
            else:
                fasta[header]+=line.strip()
        
    count_over_1000 =0 #Hold number of contigs with lenghts > 1000 bp
    over_1000 = {} #Hold fasta sequences in new list if they are > 1000 bp long
    
    #Count contigs if they are longer than 1000 bp and add those contigs to a new list
    for j,i in fasta.items(): #Iterate through every contig header and strand
        if len(i) > 1000: #If the contig is > 1000 bp, increment the count and add the header and sequence strand to a list for contigs > 1000bp
            count_over_1000 += 1
            over_1000[j] = i
    #Add all contig lengths up if they're longer than 1000 bp
    assembly_length=0
    for i in over_1000.values():
        assembly_length+=len(i)

     #Write results to miniproject.log
    os.system('echo There are '+str(count_over_1000)+' contigs \> 1000 in the assembly. >> miniproject.log')
    os.system('echo There are '+str(assembly_length)+' bp in the assembly. >> miniproject.log')
    
     #Write contigs > 1000bp to a new fasta file and only use these contigs
    with open('contigs_over1000.fasta','w') as f_out:
        for j,i in over_1000.items():
            f_out.write(j+'\n')
            f_out.write(i+'\n')
            
contig_count('contigs.fasta')


#5. Use GeneMarkS-2 and predict coding regions from the contigs with lengths over 1000 bp. 
def gms(contigs_over_1000):
    os.chdir(os.path.expanduser("~/gms2_linux_64"))
    print('Running GeneMarkS-2 to predict protein sequences from contigs_over1000.fasta')
    os.system('perl gms2.pl --seq '+os.path.expanduser("~/results/contigs_over1000.fasta")+' --genome-type bacteria --faa '+os.path.expanduser("~/results/gms.fasta")) #Use contigs_over1000 for seq, specify that its for bacteria, and output the amino acid sequences in fasta format
    print('gms2.pl has finished running')
    
gms('contigs_over1000.fasta')


#6. Blastp the predicted protein sequences against the multi-FASTA Ecoli protein sequences (Ecoli.fasta)
def blastp(gms):
    os.chdir(os.path.expanduser("~"))
    print('Running blastp on predicted protein sequences against multi-FASTA Ecoli protein sequences')
    os.system('blastp -query '+os.path.expanduser("~/results/"+gms)+' -subject '+os.path.expanduser("~/Ecoli.fasta")+' -max_hsps 1 -max_target_seqs 1 -outfmt "10 qseqid sseqid pident qcovs" -out '+os.path.expanduser("~/results/predicted_functionality.csv")) #Queries the predicted coding sequences against a multi-fasta Ecoli amino acid sequence file from prokka, takes only the best alignment for each pair, and outputs to csv with the query and subject sequence ids, percent identity, and percent query coverage
    print('Blastp has finished running')

blastp('gms.fasta')

#7. The assembled genome in RefSeq for E.coli K-12 (NC_000913) has 4140 CDS annotated. Write any discrepancies found between the RefSeq and GeneMarkS-2 predictions
def predict_comp(prediction,headers):
    os.chdir(os.path.expanduser("~/results"))
    x=[]
    results=open(prediction,'r')
    rows=csv.DictReader(results,headers,delimiter=',')
    for row in rows:
        x.append(row)
    results.close()
    return x
   
x = predict_comp('predicted_functionality.csv',['qseqid','sseqid','pident','qcovs']) #Takes in annoted predicted CDS from csv file to list and outputs discrepancy if any between the prediction and the Refseq
if len(x) > 4140:
    os.system('echo GeneMarkS-2 found '+str(len(x)-4140)+' additional CDS than the RefSeq. >> miniproject.log')
elif len(x) < 4140:
    os.system('echo GeneMarkS-2 found '+str(4140-len(x))+' less CDS than the RefSeq. >> miniproject.log')
elif len(x) == 4140:
    os.system('echo GeneMarkS-2 found the same number of CDS as the RefSeq. >> miniproject.log')
    
    
#8. Use TopHat and Cufflinks to map the reads of the E.coli transcriptome project of a K-12 derivative BW38028 and quantify their expression, respectively. Map these reads to the complete annotated genome NC_000913.
sra("SRR1411276")
os.chdir(os.path.expanduser("~"))
os.system('wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/Escherichia_coli_K_12_substr__MG1655_uid57779/NC_000913.fna')#Download ref genome
os.chdir(os.path.expanduser("~/tophat-2.1.1.Linux_x86_64"))
os.system('./bowtie2-build '+os.path.expanduser("~/NC_000913.fna")+' ~/results/EcoliK12') #Build bowtie2 index with ecoli ref genome
os.system('./tophat2 --no-novel-juncs -o '+os.path.expanduser("~/results")+' '+os.path.expanduser("~/results/EcoliK12")+' '+os.path.expanduser("~/results/SRR1411276.fastq")) #Run tophat with rna-seq against the NC_000913 index created with bowtie2-build
os.chdir(os.path.expanduser("~/cufflinks-2.2.1.Linux_x86_64"))
os.system('./cufflinks -o '+os.path.expanduser("~/results")+' -p 2 '+os.path.expanduser("~/results/accepted_hits.bam")) #Run cufflinks with .bam file from tophat to quantify the expression






#9. Parse through Cufflinks output to create "transcriptome_data.fpkm" in csv format with seqname, start, end, strand, and FPKM for each record
def CuffParse(Top_out):
    os.chdir(os.path.expanduser("~/results"))
    #Read .gtf lines into list
    with open(Top_out,'r') as f_in:
        reads = []
        for line in f_in:
            reads.append(line.split())
    #Parse .gtf file for desired output        
    parsed = []
    for i in reads:
        parsed.append([i[11].strip('";'),i[3],i[4],i[2],i[13].strip('";')])
    
    header = ['seqname','start','end','strand','FPKM']
    with open('transcriptome_data.fpkm','w') as f_out:
            write = csv.writer(f_out)
            write.writerow(header)
            write.writerows(parsed)
        
CuffParse('transcripts.gtf')
