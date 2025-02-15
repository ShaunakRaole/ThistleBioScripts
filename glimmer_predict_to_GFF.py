import sys

def glimmer_predict_to_gff(predict_file, output_file):
    with open(predict_file, "r") as in_handle, open(output_file, "w") as out_handle:
        # Write GFF header
        out_handle.write("##gff-version 3\n")
        
        current_contig = None
        
        # Parse predict file
        for line in in_handle:
            if line.startswith(">"):
                # Extract contig information
                parts = line.strip().split()
                current_contig = parts[0]
            elif line.startswith("orf"):
                # Extract gene information
                parts = line.strip().split()
                gene_id = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                strand = "+" if int(parts[3]) > 0 else "-"
                
                # Ensure start <= end
                if start > end:
                    start, end = end, start
                    strand = "-"  # Correct strand if necessary
                
                # Write GFF line
                out_handle.write(f"{current_contig}\tGlimmer\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene_id};Name={gene_id}\n")

# Usage
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.predict output.gff")
        sys.exit(1)
    
    predict_file = sys.argv[1]
    output_file = sys.argv[2]
    
    glimmer_predict_to_gff(predict_file, output_file)
