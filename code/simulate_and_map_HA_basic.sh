SEQS=$1
OUTFILE=$2
TAG=$3
TMP=tmp_HA_basic/tmp_"$TAG"
DB=$TMP/tmp_reference

# Make a directory to store the temporary files
mkdir $TMP

# Take two random sequences
python3 code/get_random_seq.py $SEQS > $TMP/tmp_query
python3 code/get_random_seq.py $SEQS > $TMP/tmp_reference
REFLEN=$(python3 code/get_ref_len.py "$TMP"/tmp_reference)
echo $TMP/tmp_reference

# Build databases
hisat2-build $DB $DB
bwa index -p $DB $DB
echo ./code/nextgenmap/bin/ngm-0.5.0/ngm
ngm -r $DB

echo "java -jar code/afg/ArtificialFastqGenerator.jar -O $TMP/afg_out_0 -R $TMP/tmp_query -F1 res/template_R1_001.fastq -F2 res/template_R2_001.fastq -RL 151 -S "$(grep ">" $TMP/tmp_query)" -URQS true -SE true -CMP 100"
# Simulate reads with AFG
java -jar code/afg/ArtificialFastqGenerator.jar -O $TMP/afg_out_0 -R $TMP/tmp_query -F1 res/template_R1_001.fastq -F2 res/template_R2_001.fastq -RL 151 -S "$(grep ">" $TMP/tmp_query)" -URQS true -SE true -CMP 100

# Mutate the reads additionally to account for experimental error and experimental noise
python3 code/read_mutator.py $TMP/afg_out_0.1.fastq > $TMP/afg_out.1.fastq
python3 code/read_mutator.py $TMP/afg_out_0.2.fastq > $TMP/afg_out.2.fastq

# Map with each program, retrieve number of mapped reads
#MINREADS=$(minimap2 -ax sr "$DB" "$TMP"/afg_out.1.fastq "$TMP"/afg_out.2.fastq | samtools view -c -F 260 -)
#NGMREADS=$(ngm -1 "$TMP"/afg_out.1.fastq -2 "$TMP"/afg_out.2.fastq -r "$DB" | samtools view -c -F 260 -)
#BWAREADSO=$(bwa mem "$DB" "$TMP"/afg_out.1.fastq "$TMP"/afg_out.2.fastq | samtools view -c -F 260 -)
#HI2READSO=$(hisat2 -x "$DB" -1 "$TMP"/afg_out.1.fastq -2 "$TMP"/afg_out.2.fastq | samtools view -c -F 260 -)

# Do the alignment
minimap2 -ax sr "$DB" "$TMP"/afg_out.1.fastq "$TMP"/afg_out.2.fastq > "$TMP"/tmp_min.bam
ngm -1 "$TMP"/afg_out.1.fastq -2 "$TMP"/afg_out.2.fastq -r "$DB" > "$TMP"/tmp_ngm.bam
bwa mem "$DB" "$TMP"/afg_out.1.fastq "$TMP"/afg_out.2.fastq > "$TMP"/tmp_bwa.bam
hisat2 -x "$DB" -1 "$TMP"/afg_out.1.fastq -2 "$TMP"/afg_out.2.fastq > "$TMP"/tmp_hi2.bam

# Retrieve the reads
samtools fastq -F 4 -0 "$TMP"/tmp_min.0.fq -1 "$TMP"/tmp_min.1.fq -2 "$TMP"/tmp_min.2.fq "$TMP"/tmp_min.bam
samtools fastq -F 4 -0 "$TMP"/tmp_ngm.0.fq -1 "$TMP"/tmp_ngm.1.fq -2 "$TMP"/tmp_ngm.2.fq "$TMP"/tmp_ngm.bam
samtools fastq -F 4 -0 "$TMP"/tmp_bwa.0.fq -1 "$TMP"/tmp_bwa.1.fq -2 "$TMP"/tmp_bwa.2.fq "$TMP"/tmp_bwa.bam
samtools fastq -F 4 -0 "$TMP"/tmp_hi2.0.fq -1 "$TMP"/tmp_hi2.1.fq -2 "$TMP"/tmp_hi2.2.fq "$TMP"/tmp_hi2.bam

# New counting approach
MINL1=$(wc -l "$TMP"/tmp_min.1.fq | cut -f 1 -d' ')
NGML1=$(wc -l "$TMP"/tmp_ngm.1.fq | cut -f 1 -d' ')
BWAL1=$(wc -l "$TMP"/tmp_bwa.1.fq | cut -f 1 -d' ')
HI2L1=$(wc -l "$TMP"/tmp_hi2.1.fq | cut -f 1 -d' ')

MINL2=$(wc -l "$TMP"/tmp_min.2.fq | cut -f 1 -d' ')
NGML2=$(wc -l "$TMP"/tmp_ngm.2.fq | cut -f 1 -d' ')
BWAL2=$(wc -l "$TMP"/tmp_bwa.2.fq | cut -f 1 -d' ')
HI2L2=$(wc -l "$TMP"/tmp_hi2.2.fq | cut -f 1 -d' ')

MINL0=$(wc -l "$TMP"/tmp_min.0.fq | cut -f 1 -d' ')
NGML0=$(wc -l "$TMP"/tmp_ngm.0.fq | cut -f 1 -d' ')
BWAL0=$(wc -l "$TMP"/tmp_bwa.0.fq | cut -f 1 -d' ')
HI2L0=$(wc -l "$TMP"/tmp_hi2.0.fq | cut -f 1 -d' ')

MINCOUNT=$(expr "$MINL1" + "$MINL2" + "$MINL0")
NGMCOUNT=$(expr "$NGML1" + "$NGML2" + "$NGML0")
BWACOUNT=$(expr "$BWAL1" + "$BWAL2" + "$BWAL0")
HI2COUNT=$(expr "$HI2L1" + "$HI2L2" + "$HI2L0")

MINREADS=$(expr "$MINCOUNT" / 4)
NGMREADS=$(expr "$NGMCOUNT" / 4)
BWAREADS=$(expr "$BWACOUNT" / 4)
HI2READS=$(expr "$HI2COUNT" / 4)

# Old counting approach
MINREADSO=$(samtools view -c -F 260 "$TMP"/tmp_min.bam)
NGMREADSO=$(samtools view -c -F 260 "$TMP"/tmp_ngm.bam)
BWAREADSO=$(samtools view -c -F 260 "$TMP"/tmp_bwa.bam)
HI2READSO=$(samtools view -c -F 260 "$TMP"/tmp_hi2.bam)

# Get the total number of reads, lines/4
R1L=$(wc -l $TMP/afg_out.1.fastq | cut -f 1 -d' ')
R1C=$(expr "$R1L" / 4)
R2L=$(wc -l $TMP/afg_out.2.fastq | cut -f 1 -d' ')
R2C=$(expr "$R2L" / 4)
NREADS=$(expr "$R1C" + "$R2C")

# Get true percentage IDs of sequences
SEQ2REFPID=$(python3 code/get_pid_between_2_seqs.py  "$TMP"/tmp_query  "$TMP"/tmp_reference)

# Get the original sequence names
REFSEQ=$(grep ">" $TMP/tmp_reference)
QUERYSEQ=$(grep ">" $TMP/tmp_query)

if [[ $R1C -ne $R2C ]]; then
    echo "NERROR! N" $R1C $R2C
fi

if [[ $MINREADSO -ne $MINREADS ]]; then
    echo "ERROR! Min" $NREADS $MINREADSO $MINREADS
fi
if [[ $BWAREADSO -ne $BWAREADS ]]; then
    echo "ERROR! BWA" $NREADS $BWAREADSO $BWAREADS
fi
if [[ $NGMREADSO -ne $NGMREADS ]]; then
    echo "ERROR! NGM" $NREADS $NGMREADSO $NGMREADS
fi
if [[ $HI2READSO -ne $HI2READS ]]; then
    echo "ERROR! HI2" $NREADS $HI2READSO $HI2READS
fi

echo "DIFF! Min" $NREADS $MINREADSO $MINREADS
echo "DIFF! BWA" $NREADS $BWAREADSO $BWAREADS
echo "DIFF! NGM" $NREADS $NGMREADSO $NGMREADS
echo "DIFF! HI2" $NREADS $HI2READSO $HI2READS


# Output
echo "RESULTS"
echo $DB,$REFSEQ,$QUERYSEQ,$SEQ2REFPID,$REFLEN,$NREADS,$MINREADS,$NGMREADS,$BWAREADS,$HI2READS
echo $DB,$REFSEQ,$QUERYSEQ,$SEQ2REFPID,$REFLEN,$NREADS,$MINREADS,$NGMREADS,$BWAREADS,$HI2READS >> $OUTFILE

# Finally remove the temporary files, don't necessarily have space, but also to stay tidy
rm -r $TMP
