#!/bin/bash

# Check arguments
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 /path/to/fasta/files /path/to/output/directory"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="${2%/}"
METRICS_FILE=$OUTPUT_DIR/"performance_metrics.tsv"

echo -e "Prefix\tTool\tRuntime (seconds)\tCPU (%)\tMemory (KB)" > "$METRICS_FILE"

# Check if cd-hit is already installed
# Uncommment the next 11 lines (highlight text and press cmd+/) only if you wish to use CD-Hit

# if ! command -v cd-hit &> /dev/null; then
#     echo "CD-HIT is not installed. Proceeding with installation..."
#     curl -L -O https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz
#     tar -xvf cd-hit-v4.8.1-2019-0228.tar.gz
#     cd cd-hit-v4.8.1-2019-0228
#     make openmp=no
#     sudo make install
#     echo "CD-HIT installation complete."
# else
#     echo "CD-HIT is already installed."
# fi
#     cd -

# Verify input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Directory $INPUT_DIR does not exist"
    exit 1
fi

# Create output directories
mkdir -p "$OUTPUT_DIR"/{ref,Prokka,Barrnap,GeMoMa,Prodigal,AMRFinder_Prodigal,AMRFinder_GeMoMa,AMRFinder_Merged,Merged,Shared_GFFs}


# Install cpanimus if it isn't already there
if ! command -v cpanm &> /dev/null; then
    echo "cpanminus not found, installing..."
    brew install cpanminus
else
    echo "cpanminus is already installed."
fi

# Install XML::Simple (Perl module)
cpan XML::Simple

# Download reference genome and annotation if missing
cd "$OUTPUT_DIR/ref"
if [ ! -f "GCF_000008805.1_ASM880v1_genomic.fna" ] || [ ! -f "GCF_000008805.1_ASM880v1_genomic.gff" ]; then
    curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/245/025/GCF_008245025.1_ASM824502v1/GCF_008245025.1_ASM824502v1_genomic.fna.gz
    curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/245/025/GCF_008245025.1_ASM824502v1/GCF_008245025.1_ASM824502v1_genomic.gff.gz
    gunzip *.gz
else
    echo "Reference genome already exists."
fi
cd -

# Check if the conda environment 'gpga_64' exists
if ! conda info --envs | grep -q "gpga_64"; then
    echo "Creating conda environment 'gpga_64'..."
    CONDA_SUBDIR=osx-64 conda create -n gpga_64 openjdk=11 -y
else
    echo "Conda environment 'gpga_64' already exists."
fi

# Activate the environment
conda activate gpga_64

# Check if 'prodigal' is installed
if ! conda list | grep -q "prodigal"; then
    echo "Installing prodigal..."
    conda install -c bioconda prodigal -y
else
    echo "Prodigal is already installed."
fi

# Check if 'prokka' is installed
if ! conda list | grep -q "prokka"; then
    echo "Installing prokka..."
    conda install -c bioconda prokka -y
else
    echo "Prokka is already installed."
fi

# Check if 'gemoma' is installed
if ! conda list | grep -q "gemoma"; then
    echo "Installing gemoma..."
    conda install -c bioconda gemoma -y
else
    echo "Gemoma is already installed."
fi

# Check if 'barrnap' is installed
if ! conda list | grep -q "barrnap"; then
    echo "Installing barrnap..."
    conda install -c bioconda -c conda-forge barrnap -y
else
    echo "Barrnap is already installed."
fi

# Check if 'bedtools' is installed
if ! conda list | grep -q "bedtools"; then
    echo "Installing bedtools..."
    conda install -c bioconda -c conda-forge bedtools -y
else
    echo "Bedtools is already installed."
fi

# Check if 'amrfinderplus' is installed
if ! conda list | grep -q "amrfinderplus"; then
    echo "Installing amrfinderplus..."
    conda install -y -c conda-forge -c bioconda --strict-channel-priority ncbi-amrfinderplus
else
    echo "amrfinderplus is already installed."
fi

amrfinder -u

# Loop through FASTA files in input directory
for file in "$INPUT_DIR"/*_BBDUK_SPADES.fasta; do
    BASE_NAME=$(basename "$file" .fasta)
    prefix=${BASE_NAME%%_*}

    # Create tool-specific directories per prefix
    mkdir -p "$OUTPUT_DIR"/{Prokka,Barrnap,GeMoMa,Prodigal,AMRFinder_Prodigal,AMRFinder_GeMoMa,AMRFinder_Merged,Merged,Shared_GFFs}/"$prefix"

    ###################################
    ############ PRODIGAL #############
    ###################################

    # Run Prodigal (macOS compatible timing)
    echo -e "Running Prodigal for $prefix.fasta"
    LOG_FILE="$OUTPUT_DIR/Prodigal/$prefix/${prefix}_prodigal_log.txt"
    /usr/bin/time -l prodigal \
        -i "$file" \
        -c \
        -m \
        -f gff \
        -o "$OUTPUT_DIR/Prodigal/$prefix/${prefix}_cds.gff" \
        -d "$OUTPUT_DIR/Prodigal/$prefix/${prefix}_cds_nucleotide.fasta" \
        -a "$OUTPUT_DIR/Prodigal/$prefix/${prefix}_cds_protein.fasta" \
        >"$LOG_FILE" 2>&1

    # Copying GFFs to shared folder
    cp "$OUTPUT_DIR/Prodigal/$prefix/${prefix}_cds.gff" "$OUTPUT_DIR/Shared_GFFs/$prefix/Prodigal_$prefix.gff"
    gzip "$OUTPUT_DIR/Shared_GFFs/$prefix/Prodigal_$prefix.gff"

    # Initialize variables
    TOTAL_CDS=0
    TOTAL_LENGTH=0

    while read -r line; do
        # Skip comments
        if [[ $line == \#* ]]; then
            continue
        fi
        # Extract start and end positions (columns 4 and 5 in GFF format)
        START=$(echo "$line" | awk '{print $4}')
        END=$(echo "$line" | awk '{print $5}')
        
        # Calculate CDS length
        CDS_LENGTH=$((END - START + 1))
        
        # Update totals
        TOTAL_CDS=$((TOTAL_CDS + 1))
        TOTAL_LENGTH=$((TOTAL_LENGTH + CDS_LENGTH))
    done < "$OUTPUT_DIR/Prodigal/$prefix/${prefix}_cds.gff"

    # Compute average CDS length
    if [[ $TOTAL_CDS -gt 0 ]]; then
        AVG_CDS_LENGTH=$(echo "$TOTAL_LENGTH / $TOTAL_CDS" | bc)
    else
        AVG_CDS_LENGTH=0
    fi

    echo -e "$prefix\t$TOTAL_CDS\t$AVG_CDS_LENGTH" >> $OUTPUT_DIR/Prodigal/$prefix/${prefix}_cds_metrics.tsv
    echo -e "Finished running Prodigal for $prefix.fasta"

    # Run CD-HIT on Prodigal's protein sequences
    # Uncommment the next 8 lines (highlight text and press cmd+/) only if you wish to use CD-Hit
    # LOG_CDHIT="$OUTPUT_DIR/Prodigal/$prefix/${prefix}_cdhit_log.txt"
    # /usr/bin/time -l cd-hit \
    #     -i "$OUTPUT_DIR/Prodigal/$prefix/${prefix}_cds_protein.fasta" \
    #     -o "$OUTPUT_DIR/Prodigal/$prefix/${prefix}_cds_protein_cdhit.fasta" \
    #     -c 0.9 \
    #     -n 5 \
    #     -d 0 \
    #     -T 10 \
    #     > "$LOG_CDHIT" 2>&1

    # Extract metrics (macOS compatible parsing)
    PRODIGAL_RUNTIME=$(grep -E "user" "$LOG_FILE" | tail -n 1 | awk '{print $3}')
    PRODIGAL_MEM_KB=$(grep -E "maximum resident set size" "$LOG_FILE" | awk '{print $1}' | tr -d '\r' | awk '{print $1/1048576}')
    USER_TIME=$(grep -E "user" "$LOG_FILE" | tail -n 1 | awk '{print $3}')
    SYSTEM_TIME=$(grep "sys" "$LOG_FILE" | awk '{print $5}' | head -n 1)
    # Calculate CPU usage
    TOTAL_CPU_TIME=$(echo "$USER_TIME + $SYSTEM_TIME" | bc -l)
    CPU_USAGE=$(echo "scale=2; ($TOTAL_CPU_TIME / $PRODIGAL_RUNTIME) * 100" | bc -l)
    echo -e "$prefix\tProdigal\t$PRODIGAL_RUNTIME\t${CPU_USAGE}\t$PRODIGAL_MEM_KB" >>"$METRICS_FILE"

    
    echo -e "Running AMRFinder for Prodigal results of $prefix.fasta"
    # Run AMRFinderPlus
    AMR_LOG_FILE="$OUTPUT_DIR/AMRFinder_Prodigal/$prefix/${prefix}_amrfinder_log.txt"
    /usr/bin/time -l amrfinder -n "$OUTPUT_DIR/Prodigal/$prefix/${prefix}_cds_nucleotide.fasta" -O Campylobacter -o "$OUTPUT_DIR/AMRFinder_Prodigal/$prefix/${prefix}_amr_results.tsv" > "$AMR_LOG_FILE" 2>&1
    awk -v prefix="$prefix" 'BEGIN {OFS="\t"} NR==1 {print "SampleID", $0; next} {print prefix, $0}' \
    "$OUTPUT_DIR/AMRFinder_Prodigal/$prefix/${prefix}_amr_results.tsv" > "$OUTPUT_DIR/AMRFinder_Prodigal/$prefix/${prefix}_amr_results.tmp" && \
    mv "$OUTPUT_DIR/AMRFinder_Prodigal/$prefix/${prefix}_amr_results.tmp" "$OUTPUT_DIR/AMRFinder_Prodigal/$prefix/${prefix}_amr_results.tsv"

    AMR_RUNTIME=$(grep -E "user" "$LOG_FILE" | tail -n 1 | awk '{print $3}')
    AMR_MEM_KB=$(grep -E "maximum resident set size" "$AMR_LOG_FILE" | awk '{print $1}' | tr -d '\r' | awk '{print $1/1048576}')
    USER_TIME=$(grep -E "user" "$LOG_FILE" | tail -n 1 | awk '{print $3}')
    SYSTEM_TIME=$(grep -E "sys" "$AMR_LOG_FILE" | awk '{print $5}' | head -n 1)
    # Calculate CPU usage
    TOTAL_CPU_TIME=$(echo "$USER_TIME + $SYSTEM_TIME" | bc -l)
    CPU_USAGE=$(echo "scale=2; ($TOTAL_CPU_TIME / $AMR_RUNTIME) * 100" | bc -l)
    echo -e "$prefix\tAMRFinder_Merged\t$AMR_RUNTIME\t$CPU_USAGE\t$AMR_MEM_KB" >>"$METRICS_FILE"
    echo -e "Finished running AMRFinder for Prodigal results of $prefix.fasta"


    ###################################
    ############# PROKKA ##############
    ###################################


    echo -e "Running Prokka for $prefix.fasta"
    # Run Prokka (assuming prokka is installed via conda or brew)
    LOG_FILE="$OUTPUT_DIR/Prokka/$prefix/${prefix}_prokka_log.txt"
    /usr/bin/time -l prokka --outdir "$OUTPUT_DIR/Prokka/$prefix" --prefix "$prefix" --force --kingdom Bacteria --genus Campylobacter --cpus 10 "$file" >"$LOG_FILE" 2>&1
    echo -e "Finished running Prokka for $prefix.fasta"

    # Copying GFFs to shared folder
    cp "$OUTPUT_DIR/Prokka/$prefix/${prefix}.gff" "$OUTPUT_DIR/Shared_GFFs/$prefix/Prokka_$prefix.gff"
    gzip "$OUTPUT_DIR/Shared_GFFs/$prefix/Prokka_$prefix.gff"

    PROKKA_RUNTIME=$(grep -E "user" "$LOG_FILE" | tail -n 1 | awk '{print $3}')
    PROKKA_MEM_KB=$(grep -E "maximum resident set size" "$LOG_FILE" | awk '{print $1}' | tr -d '\r' | awk '{print $1/1048576}')
    USER_TIME=$(grep -E "user" "$LOG_FILE" | tail -n 1 | awk '{print $3}')
    SYSTEM_TIME=$(grep -E "sys" "$LOG_FILE" | awk '{print $5}' | head -n 1)
    # Calculate CPU usage
    TOTAL_CPU_TIME=$(echo "$USER_TIME + $SYSTEM_TIME" | bc -l)
    CPU_USAGE=$(echo "scale=2; ($TOTAL_CPU_TIME / $PROKKA_RUNTIME) * 100" | bc -l)
    echo -e "$prefix\tProkka\t$PROKKA_RUNTIME\t$CPU_USAGE\t$PROKKA_MEM_KB" >>"$METRICS_FILE"


    ####################################
    ############# Barrnap ##############
    ####################################

    # Run Barrnap (assuming installed via conda or brew)
    echo -e "Running Barrnap for $prefix.fasta"
    LOG_FILE="$OUTPUT_DIR/Barrnap/$prefix/barrnap.log"
    GFF_OUTPUT="$OUTPUT_DIR/Barrnap/$prefix/${prefix}.gff"
    FASTA_OUTPUT="$OUTPUT_DIR/Barrnap/$prefix/${prefix}_16S.fasta"
    mkdir -p "$(dirname "$FASTA_OUTPUT")"

    /usr/bin/time -l barrnap --kingdom bac --threads 10 --outseq "$FASTA_OUTPUT" "$file" >"$GFF_OUTPUT" 2>"$OUTPUT_DIR/Barrnap/$prefix/barrnap.log"

    echo -e "Finished running Barrnap for $prefix.fasta"

    # Copying GFFs to shared folder
    cp "$OUTPUT_DIR/Barrnap/$prefix/${prefix}.gff" "$OUTPUT_DIR/Shared_GFFs/$prefix/Barrnap_$prefix.gff"
    gzip "$OUTPUT_DIR/Shared_GFFs/$prefix/Barrnap_$prefix.gff"

    # Initialize variables
    TOTAL_CDS=0
    TOTAL_LENGTH=0

    while read -r line; do
        # Skip comments
        if [[ $line == \#* ]]; then
            continue
        fi
        # Extract start and end positions (columns 4 and 5 in GFF format)
        START=$(echo "$line" | awk '{print $4}')
        END=$(echo "$line" | awk '{print $5}')
        
        # Calculate CDS length
        CDS_LENGTH=$((END - START + 1))
        
        # Update totals
        TOTAL_CDS=$((TOTAL_CDS + 1))
        TOTAL_LENGTH=$((TOTAL_LENGTH + CDS_LENGTH))
    done < "$GFF_OUTPUT"

    # Compute average CDS length
    if [[ $TOTAL_CDS -gt 0 ]]; then
        AVG_CDS_LENGTH=$(echo "$TOTAL_LENGTH / $TOTAL_CDS" | bc)
    else
        AVG_CDS_LENGTH=0
    fi

    echo -e "$prefix\t$TOTAL_CDS\t$AVG_CDS_LENGTH" >> $OUTPUT_DIR/Barrnap/$prefix/${prefix}_rna_metrics.tsv

    BARRNAP_RUNTIME=$(grep -E "user" "$LOG_FILE" | tail -n 1 | awk '{print $3}')
    BARRNAP_MEM_KB=$(grep -E "maximum resident set size" "$LOG_FILE" | awk '{print $1}' | tr -d '\r' | awk '{print $1/1048576}')
    USER_TIME=$(grep -E "user" "$LOG_FILE" | tail -n 1 | awk '{print $3}')
    SYSTEM_TIME=$(grep -E "sys" "$LOG_FILE" | awk '{print $5}' | head -n 1)
    # Calculate CPU usage
    TOTAL_CPU_TIME=$(echo "$USER_TIME + $SYSTEM_TIME" | bc -l)
    CPU_USAGE=$(echo "scale=2; ($TOTAL_CPU_TIME / $BARRNAP_RUNTIME) * 100" | bc -l)
    echo -e "$prefix\tBarrnap\t$BARRNAP_RUNTIME\t$CPU_USAGE\t$BARRNAP_MEM_KB" >>"$METRICS_FILE"

    # Run Bedtools 
    echo -e "Converting Barrnap GFF to FASTA for $prefix.fasta"
    bedtools getfasta -fi "$file" -bed "$GFF_OUTPUT" >"$FASTA_OUTPUT" 2> "$LOG_FILE" 2>&1
    echo -e "Done!"

    # Check Bedtools execution status
    if [ $? -ne 0 ]; then
        echo "Bedtools failed for '$GFF_OUTPUT'. Check logs."
        exit 1
    fi

    ###################################
    ############# GeMoMa ##############
    ###################################

    echo -e "Running GeMoMa for $prefix.fasta, this may take a while!"
    # Run GeMoMa (macOS compatible timing)
    GEMOMA_OUT="$OUTPUT_DIR/GeMoMa/$prefix"
    mkdir -p "$GEMOMA_OUT"
    LOG_FILE="$OUTPUT_DIR/GeMoMa/$prefix/${prefix}_gemoma_log.txt"
    /usr/bin/time -l gemoma GeMoMaPipeline \
        g="$OUTPUT_DIR/ref/GCF_008245025.1_ASM824502v1_genomic.fna" \
        a="$OUTPUT_DIR/ref/GCF_008245025.1_ASM824502v1_genomic.gff" \
        t="$file" \
        tblastn=true \
        threads=10 \
        o=true \
        i=1 \
        outdir="$GEMOMA_OUT" \
        AnnotationFinalizer.r=NO > "$LOG_FILE" 2>&1

    # Copying GFFs to shared folder
    cp "$OUTPUT_DIR/GeMoMa/$prefix/final_annotation.gff" "$OUTPUT_DIR/Shared_GFFs/$prefix/GeMoMa_$prefix.gff"
    gzip "$OUTPUT_DIR/Shared_GFFs/$prefix/GeMoMa_$prefix.gff"

    # Initialize variables
    TOTAL_CDS=0
    TOTAL_LENGTH=0

    while read -r line; do
        # Skip comments
        if [[ $line == \#* ]]; then
            continue
        fi
        # Extract start and end positions (columns 4 and 5 in GFF format)
        START=$(echo "$line" | awk '{print $4}')
        END=$(echo "$line" | awk '{print $5}')
        
        # Calculate CDS length
        CDS_LENGTH=$((END - START + 1))
        
        # Update totals
        TOTAL_CDS=$((TOTAL_CDS + 1))
        TOTAL_LENGTH=$((TOTAL_LENGTH + CDS_LENGTH))
    done < "$OUTPUT_DIR/GeMoMa/$prefix/final_annotation.gff"

    # Compute average CDS length
    if [[ $TOTAL_CDS -gt 0 ]]; then
        AVG_CDS_LENGTH=$(echo "$TOTAL_LENGTH / $TOTAL_CDS" | bc)
    else
        AVG_CDS_LENGTH=0
    fi

    echo -e "$prefix\t$TOTAL_CDS\t$AVG_CDS_LENGTH" >> $OUTPUT_DIR/GeMoMa/$prefix/${prefix}_cds_metrics.tsv


    # Extract metrics (macOS compatible parsing)
    GEMOMA_RUNTIME=$(grep -E "user" "$LOG_FILE" | tail -n 1 | awk '{print $3}')
    GEMOMA_MEM_GB=$(grep -E "maximum resident" "$LOG_FILE" | awk '{print $1}' | tail -n 1 | awk '{print $1}' | tr -d '\r' | awk '{print $1/1048576}')
    USER_TIME=$(grep -E "user" "$LOG_FILE" | tail -n 1 | awk '{print $3}')
    SYSTEM_TIME=$(grep -E "sys" "$LOG_FILE" | awk '{print $5}' | head -n 1)
    # Calculate CPU usage
    TOTAL_CPU_TIME=$(echo "$USER_TIME + $SYSTEM_TIME" | bc -l)
    CPU_USAGE=$(echo "scale=2; ($TOTAL_CPU_TIME / $GEMOMA_RUNTIME) * 100" | bc -l)
    echo -e "$prefix\tGeMoMa\t$GEMOMA_RUNTIME\t$CPU_USAGE\t$GEMOMA_MEM_GB" >> "$METRICS_FILE"

    echo "GeMoMa completed for $prefix. Results saved in $OUTPUT_DIR/GeMoMa/$prefix."

    echo -e "Running AMRFinder for GeMoMa result of $prefix.fasta!"
    # Run AMRFinderPlus
    AMR_LOG_FILE="$OUTPUT_DIR/AMRFinder_GeMoMa/$prefix/${prefix}_amrfinder_log.txt"
    /usr/bin/time -l amrfinder -p "$OUTPUT_DIR/GeMoMa/$prefix/predicted_proteins.fasta" -O Campylobacter -o "$OUTPUT_DIR/AMRFinder_GeMoMa/$prefix/${prefix}_amr_results.tsv" > "$AMR_LOG_FILE" 2>&1
    awk -v prefix="$prefix" 'BEGIN {OFS="\t"} NR==1 {print "SampleID", $0; next} {print prefix, $0}' \
    "$OUTPUT_DIR/AMRFinder_GeMoMa/$prefix/${prefix}_amr_results.tsv" > "$OUTPUT_DIR/AMRFinder_GeMoMa/$prefix/${prefix}_amr_results.tmp" && \
    mv "$OUTPUT_DIR/AMRFinder_GeMoMa/$prefix/${prefix}_amr_results.tmp" "$OUTPUT_DIR/AMRFinder_GeMoMa/$prefix/${prefix}_amr_results.tsv"


    echo -e "Finished running AMRFinder for GeMoMa result of $prefix.fasta!"

    # BEDTOOLS MERGE

    echo -e "Now merging your Prodigal and GeMoMa results for $prefix.fasta!"
    # Merge Prodigal and GeMoMa GFF files using bedtools merge
    LOG_FILE="$OUTPUT_DIR/Merged/$prefix/${prefix}_merged_log.txt"
    MERGED_GFF="$OUTPUT_DIR/Merged/$prefix/${prefix}_merged.gff"
    bedtools merge -i <(cat "$OUTPUT_DIR/Prodigal/$prefix/${prefix}_cds.gff" "$OUTPUT_DIR/GeMoMa/$prefix/final_annotation.gff") > "$MERGED_GFF"

    # Convert merged GFF to FASTA using bedtools getfasta
    MERGED_FASTA="$OUTPUT_DIR/Merged/$prefix/${prefix}_merged_cds.fasta"
    bedtools getfasta -fi "$file" -bed "$MERGED_GFF" > "$MERGED_FASTA" 

    echo "Merged GFF and FASTA files generated for $prefix. Results saved in $OUTPUT_DIR/Merged/$prefix."

    cp "$OUTPUT_DIR/Merged/$prefix/${prefix}_merged.gff" "$OUTPUT_DIR/Shared_GFFs/$prefix/Merged_$prefix.gff"
    gzip "$OUTPUT_DIR/Shared_GFFs/$prefix/Merged_$prefix.gff"

    # Initialize variables
    TOTAL_CDS=0
    TOTAL_LENGTH=0

    while read -r line; do
        # Skip comments
        if [[ $line == \#* ]]; then
            continue
        fi
        # Extract start and end positions (columns 2 and 3 in GFF format)
        START=$(echo "$line" | awk '{print $2}')
        END=$(echo "$line" | awk '{print $3}')
        
        # Calculate CDS length
        CDS_LENGTH=$((END - START + 1))
        
        # Update totals
        TOTAL_CDS=$((TOTAL_CDS + 1))
        TOTAL_LENGTH=$((TOTAL_LENGTH + CDS_LENGTH))
    done < "$MERGED_GFF"

    # Compute average CDS length
    if [[ $TOTAL_CDS -gt 0 ]]; then
        AVG_CDS_LENGTH=$(echo "$TOTAL_LENGTH / $TOTAL_CDS" | bc)
    else
        AVG_CDS_LENGTH=0
    fi

    echo -e "$prefix\t$TOTAL_CDS\t$AVG_CDS_LENGTH" >> $OUTPUT_DIR/Merged/$prefix/${prefix}_cds_metrics.tsv

    # Run AMRFinderPlus
    AMR_LOG_FILE="$OUTPUT_DIR/AMRFinder_Merged/$prefix/${prefix}_amrfinder_log.txt"
    /usr/bin/time -l amrfinder -n "$OUTPUT_DIR/Merged/$prefix/${prefix}_merged_cds.fasta" -O Campylobacter -o "$OUTPUT_DIR/AMRFinder_Merged/$prefix/${prefix}_amr_results.tsv" > "$AMR_LOG_FILE" 2>&1
    awk -v prefix="$prefix" 'BEGIN {OFS="\t"} NR==1 {print "SampleID", $0; next} {print prefix, $0}' \
    "$OUTPUT_DIR/AMRFinder_Merged/$prefix/${prefix}_amr_results.tsv" > "$OUTPUT_DIR/AMRFinder_Merged/$prefix/${prefix}_amr_results.tmp" && \
    mv "$OUTPUT_DIR/AMRFinder_Merged/$prefix/${prefix}_amr_results.tmp" "$OUTPUT_DIR/AMRFinder_Merged/$prefix/${prefix}_amr_results.tsv"

done

echo -e "Preparing your final TSVs!"


# PRODIGAL CDS CALCULATION
Prodigal_MERGED_TSV="$OUTPUT_DIR/Prodigal/Prodigal_CDS_Metrics.tsv"

# Loop through each CDS output TSV and append its contents 
for tsv in "$OUTPUT_DIR"/Prodigal/*/*_cds_metrics.tsv; do
    cat "$tsv" >> "$Prodigal_MERGED_TSV"
done

echo "Merged CDS metrics saved to $Prodigal_MERGED_TSV"

# Barrnap 16s rRNA CALCULATION
Barrnap_MERGED_TSV="$OUTPUT_DIR/Barrnap/Barrnap_rna_Metrics.tsv"

# Loop through each CDS output TSV and append its contents 
for tsv in "$OUTPUT_DIR"/Barrnap/*/*_rna_metrics.tsv; do
    cat "$tsv" >> "$Barrnap_MERGED_TSV"
done

echo "Merged rRNA metrics saved to $Barrnap_MERGED_TSV"

# GeMoMA CDS Calculation
GeMoMa_MERGED_TSV="$OUTPUT_DIR/GeMoMa/GeMoMa_CDS_Metrics.tsv"

# Loop through each CDS output TSV and append its contents 
for tsv in "$OUTPUT_DIR"/GeMoMa/*/*_cds_metrics.tsv; do
    cat "$tsv" >> "$GeMoMa_MERGED_TSV"
done

echo "Merged CDS metrics saved to $GeMoMa_MERGED_TSV"


# Merged CDS Calculation
Merged_MERGED_TSV="$OUTPUT_DIR/Merged/Merged_CDS_Metrics.tsv"

# Loop through each CDS output TSV and append its contents 
for tsv in "$OUTPUT_DIR"/Merged/*/*_cds_metrics.tsv; do
    cat "$tsv" >> "$Merged_MERGED_TSV"
done

echo "Merged CDS metrics saved to $Merged_MERGED_TSV"

## Prodigal Annotation
# Define the merged output file
PAnn_MERGED_TSV="$OUTPUT_DIR/AMRFinder_Prodigal/MergedAnnotations.tsv"

# Get the header from the first TSV file (if available)
first_file=$(find "$OUTPUT_DIR/AMRFinder_Prodigal" -name "*_amr_results.tsv" | head -n 1)
if [[ -f "$first_file" ]]; then
    head -n 1 "$first_file" > "$PAnn_MERGED_TSV"
fi

# Loop through each AMRFinder output TSV and append its contents (excluding the header)
for tsv in "$OUTPUT_DIR"/AMRFinder_Prodigal/*/*_amr_results.tsv; do
    tail -n +2 "$tsv" >> "$PAnn_MERGED_TSV"
done

echo "Merged AMRFinder annotations saved to $PAnn_MERGED_TSV"

## GeMoMa Annotation
# Define the merged output file
GAnn_MERGED_TSV="$OUTPUT_DIR/AMRFinder_GeMoMa/MergedAnnotations.tsv"

# Get the header from the first TSV file (if available)
first_file=$(find "$OUTPUT_DIR/AMRFinder_GeMoMa" -name "*_amr_results.tsv" | head -n 1)
if [[ -f "$first_file" ]]; then
    head -n 1 "$first_file" > "$GAnn_MERGED_TSV"
fi

# Loop through each AMRFinder output TSV and append its contents (excluding the header)
for tsv in "$OUTPUT_DIR"/AMRFinder_GeMoMa/*/*_amr_results.tsv; do
    tail -n +2 "$tsv" >> "$GAnn_MERGED_TSV"
done

echo "Merged AMRFinder annotations saved to $GAnn_MERGED_TSV"

## Merged Annotation
# Define the merged output file
MAnn_MERGED_TSV="$OUTPUT_DIR/AMRFinder_Merged/MergedAnnotations.tsv"

# Get the header from the first TSV file (if available)
first_file=$(find "$OUTPUT_DIR/AMRFinder_Merged" -name "*_amr_results.tsv" | head -n 1)
if [[ -f "$first_file" ]]; then
    head -n 1 "$first_file" > "$MAnn_MERGED_TSV"
fi

# Loop through each AMRFinder output TSV and append its contents (excluding the header)
for tsv in "$OUTPUT_DIR"/AMRFinder_Merged/*/*_amr_results.tsv; do
    tail -n +2 "$tsv" >> "$MAnn_MERGED_TSV"
done

echo "Merged AMRFinder annotations saved to $MAnn_MERGED_TSV"

echo -e "Done preparing your final TSVs, that was fun!"

echo "End of script\n\nFun Fact: Did you know that some bacteria can have genomes as small as 160,000 base pairs, while others can have over 13 million base pairs? The largest bacterial genomes are often found in symbiotic or free-living bacteria, highlighting the diversity of bacterial life!"
