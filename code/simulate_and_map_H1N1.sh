DB=$1
SEQS=$2
N=$3
OUTFILE=$4
TAG=$5
TMP=tmp_H1N1/tmp_"$TAG"

# Make a directory to store the temporary files
mkdir $TMP

# Fetch the Nth sequence
echo "python3 code/get_seq_by_index.py $SEQS $N > $TMP/tmp"
python3 code/get_seq_by_index.py $SEQS $N > $TMP/tmp
REFLEN=$(python3 code/get_ref_len.py "$TMP"/tmp)

# Simulate reads with AFG
java -jar code/afg/ArtificialFastqGenerator.jar -O $TMP/afg_out_0 -R $TMP/tmp -F1 res/template_R1_001.fastq -F2 res/template_R2_001.fastq -RL 151 -S "$(grep ">" $TMP/tmp)" -URQS true -SE true -CMP 100

# Mutate the reads additionally to account for experimental error and experimental noise
python3 code/read_mutator.py $TMP/afg_out_0.1.fastq > $TMP/afg_out.1.fastq
python3 code/read_mutator.py $TMP/afg_out_0.2.fastq > $TMP/afg_out.2.fastq

# Map with each program, retrieve number of mapped reads

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

# Get the total number of reads, lines/4
R1L=$(wc -l $TMP/afg_out.1.fastq | cut -f 1 -d' ')
R2L=$(wc -l $TMP/afg_out.2.fastq | cut -f 1 -d' ')
RL=$(expr "$R1L" + "$R2L")
NREADS=$(expr "$RL" / 4)

# Get the original sequence name
SEQ=$(grep ">" $TMP/tmp)

# Get the pid
SEQ2DBPID=$(python3 code/get_pid_between_2_seqs.py  "$DB"  "$TMP"/tmp)

# Output
echo $DB,$SEQ,$SEQ2DBPID,$REFLEN,$NREADS,$MINREADS,$NGMREADS,$BWAREADS,$HI2READS
echo $DB,$SEQ,$SEQ2DBPID,$REFLEN,$NREADS,$MINREADS,$NGMREADS,$BWAREADS,$HI2READS >> $OUTFILE

# Finally remove the temporary files, don't necessarily have space, but also to stay tidy
rm -r $TMP
