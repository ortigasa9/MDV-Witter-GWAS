####################
#Author(s):Alejandro Ortigas-Vasquez, David A. Kennedy
#Date:April 23, 2026
#Description:R script for generating variant table for GWAS
#Inputs: Multiple sequence alignment in .fasta format
#Outputs:Variant tables (SNP_Calls.txt,SNPLocations.txt)
####################

# Extract SNP table from alignment
FastaFiles=read.table("data/Witter_Masked_Trimmed_Alignment.fasta", sep="", skip=0, fill=TRUE, stringsAsFactors=FALSE)
FastaFiles=FastaFiles[,1]

GenomeStarts=grep(">", FastaFiles)
GenomeNames=FastaFiles[grep(">", FastaFiles)]
GenomeNames=gsub(">", "", GenomeNames)
GenomeNames=gsub("-", "_", GenomeNames)
NumGenomes=length(GenomeStarts)

NormalLength=nchar(FastaFiles[GenomeStarts[1]+1])
FinalLength=nchar(FastaFiles[GenomeStarts[2]-1])

FastaLines=GenomeStarts[2]-GenomeStarts[1]-1

GenomeLength=NormalLength*(FastaLines-1) + FinalLength

Genomes=matrix("", NumGenomes, GenomeLength)

for (j in 1:NumGenomes)
{
    for (i in (GenomeStarts[j]+1):(GenomeStarts[j]+FastaLines-1))
    {
        Genomes[j,(i-(GenomeStarts[j]+1))*NormalLength + (1:NormalLength)] = unlist(strsplit(FastaFiles[i], ""))
    }
    
    Genomes[j,(FastaLines-1)*NormalLength + 1:FinalLength] = unlist(strsplit(FastaFiles[i+1], ""))
}

PotentialSNPs=numeric()
for (i in 1:GenomeLength)
{
    Temp=names(table(Genomes[,i]))
    if(sum((Temp=="A") + (Temp=="T") + (Temp=="C") + (Temp=="G"))>1)
    {
        PotentialSNPs[length(PotentialSNPs)+1]=i
    }
}

Data=Genomes[,PotentialSNPs]
colnames(Data)=PotentialSNPs

QuestionableSNPs=numeric()
for(i in 1:ncol(Data))
{
    if (sum((Data[,i]=="A") +(Data[,i]=="T") + (Data[,i]=="C") +(Data[,i]=="G"))<(NumGenomes))
    {
        QuestionableSNPs[length(QuestionableSNPs)+1]=i
    }
}

Data=Data[,-QuestionableSNPs]

# Create repeat (132bp, Meq-PRR) and INDEL allele matrix
# Nucleotides shown are stand-ins for biallelic genotypes
Temp_132<-matrix(c("C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "C", "T", "T", "T", "T", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C"))
colnames(Temp_132)=114491
Temp_Meq<-matrix(c("A", "A", "A", "A", "G", "G", "G", "G", "G", "A", "G", "A", "A", "A", "A", "G", "A", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "G", "A", "A", "A", "G", "G", "A", "G"))
colnames(Temp_Meq)=120556
Temp_INDEL1<-matrix(c("A", "A", "A", "A",	"A", "A", "A", "A",	"A", "A", "A", "A",	"A", "A",	"A", "A",	"A", "T", "T", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "T", "T", "T", "T", "T", "T", "T", "T" ,"T" ,"T", "T", "T", "T", "T", "A", "A", "A", "T", "T", "A", "A"))
colnames(Temp_INDEL1)=5195
Temp_INDEL2<-matrix(c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "T", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A"))
colnames(Temp_INDEL2)=13268
Temp_INDEL3<-matrix(c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "T", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A"))
colnames(Temp_INDEL3)=16607
Temp_INDEL4<-matrix(c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "T", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A"))
colnames(Temp_INDEL4)=40923
Temp_INDEL5<-matrix(c("A", "T", "T", "T", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A"))
colnames(Temp_INDEL5)=52164
Temp_INDEL6<-matrix(c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "T", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A"))
colnames(Temp_INDEL6)=61127
Temp_INDEL7<-matrix(c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "T", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A"))
colnames(Temp_INDEL7)=62959
Temp_INDEL8<-matrix(c("A", "A", "A", "A", "A", "T", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "T"))
colnames(Temp_INDEL8)=65991
Temp_INDEL9<-matrix(c("A", "A", "A", "A", "A", "A", "A", "A", "A", "T", "A", "T", "A", "A", "A", "A", "T", "T", "T", "T", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "T", "A", "A", "A", "T", "A"))
colnames(Temp_INDEL9)=79045
Temp_INDEL10<-matrix(c("T", "T", "T", "T", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A","A", "A", "A", "A", "A", "A", "A", "T", "A", "A", "A", "A", "A", "A"))
colnames(Temp_INDEL10)=86489
Temp_INDEL11<-matrix(c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "T", "A", "A", "A"))
colnames(Temp_INDEL11)=96737
Temp_INDEL12<-matrix(c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "T", "A", "A", "A", "A", "A", "A", "A"))
colnames(Temp_INDEL12)=108503
Temp_INDEL13<-matrix(c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "T", "T", "A", "A", "A", "A"))
colnames(Temp_INDEL13)=117464
Temp_INDEL14<-matrix(c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "T", "A", "A", "A", "A", "A", "A", "A", "A", "T", "T", "A", "A", "A", "A"))
colnames(Temp_INDEL14)=122569
Temp_INDEL15<-matrix(c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "T", "A", "A", "A"))
colnames(Temp_INDEL15)=136323
Temp_INDEL16<-matrix(c("A", "T", "T", "T", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A"))
colnames(Temp_INDEL16)=140290
Temp_INDEL17<-matrix(c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "T", "A", "A", "A", "A", "A", "A"))
colnames(Temp_INDEL17)=144847 
                      
# Append repeat and INDEL alleles to SNP table 
Data<-cbind(Data, Temp_132, Temp_Meq, Temp_INDEL1, Temp_INDEL2, Temp_INDEL3, Temp_INDEL4, Temp_INDEL5, Temp_INDEL6, Temp_INDEL7, Temp_INDEL8, Temp_INDEL9, Temp_INDEL10, Temp_INDEL11, Temp_INDEL12, Temp_INDEL13, Temp_INDEL14, Temp_INDEL15, Temp_INDEL16, Temp_INDEL17)
SNPLocations=colnames(Data)
Temp=unlist(strsplit(GenomeNames, "_"))
StrainNames=Temp[seq(1,length(Temp),5)]
Pathotypes=Temp[seq(2,length(Temp),5)]
VirulenceRank=Temp[seq(3,length(Temp),5)]
Year=Temp[seq(4,length(Temp),5)]
Location=Temp[seq(5,length(Temp),5)]
OutputData=cbind(StrainNames, Pathotypes, VirulenceRank, Year, Location, Data)

# Remove multiallelic sites + hompolymer sites
# These are the genomic POSITIONS to remove, not column indices
positions_to_remove = c(10136, 110755, 114657, 114675, 119798, 120690, 135456, 135870, 135871, 135872, 135873, 135874, 136739)

# Find which SNPs match these positions
snps_to_keep = !(SNPLocations %in% positions_to_remove)

# Apply the filter
TempData = OutputData[, c(rep(TRUE, 5), snps_to_keep)]  # Keep all 5 metadata cols + filtered SNPs
TempLocations = SNPLocations[snps_to_keep]
TempData=TempData[order(TempData[,2]),]

# Generate output
write.table(TempData, file="data/SNPCalls.txt", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(TempLocations, file="data/SNPLocations.txt", sep=",", quote=FALSE, row.names=FALSE, col.names=FALSE)

