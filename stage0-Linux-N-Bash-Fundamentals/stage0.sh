# Project 1
echo "StMark"                     # Prints the string "StMark" to the terminal
mkdir StMark                     # Creates a new directory named StMark in the current location
mkdir biocomputing && cd biocomputing  # Creates a directory named biocomputing and changes to it in one command
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna  # Downloads the wildtype.fna file from the specified GitHub URL
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk  # Downloads the wildtype.gbk file from the specified GitHub URL
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk  # Downloads the wildtype.gbk file again, creating a duplicate (wildtype.gbk.1)
mv wildtype.fna StMark/         # Moves the wildtype.fna file to the StMark directory
rm wildtype.gbk                 # Deletes the original wildtype.gbk file, leaving the duplicate (wildtype.gbk.1)
grep -c "tatatata" ../StMark/wildtype.fna  # Counts the number of lines containing "tatatata" in the wildtype.fna file
grep -c "tata" ../StMark/wildtype.fna      # Counts the number of lines containing "tata" in the wildtype.fna file
grep "tatatata" ../StMark/wildtype.fna > mutant_lines.txt  # Finds all lines with "tatatata" and saves them to mutant_lines.txt
grep -v "^LOCUS" wildtype.gbk.1 | wc -l    # Counts lines in wildtype.gbk.1, excluding the header line starting with "LOCUS"
head -n 1 wildtype.gbk.1 | grep -oP '(?<=LOCUS\s+)[^\s]+' | awk '{print $2}'  # Extracts the sequence length from the LOCUS tag in the first line of wildtype.gbk.1 (may fail due to regex compatibility)
head -n 1 wildtype.gbk.1 | awk '{for(i=2;i<=NF;i++) if ($i ~ /^[0-9]+$/) {print $i; exit}}'  # Extracts the sequence length as a number from the LOCUS tag in the first line of wildtype.gbk.1
grep -m 1 "^SOURCE" wildtype.gbk.1 | cut -d' ' -f2-  # Extracts the source organism from the first SOURCE line in wildtype.gbk.1
grep -oP '(?<=/gene=")[^"]*' wildtype.gbk.1  # Extracts all gene names between /gene=" and the closing quote in wildtype.gbk.1
clear                          # Clears the terminal screen
history                        # Displays the command history of the current session
ls -l ../StMark/               # Lists the files in the StMark directory in long format
ls -l .                        # Lists the files in the current directory (biocomputing) in long format

# Project 2
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh  # Downloads the Miniconda installer for Linux
bash Miniconda3-latest-Linux-x86_64.sh  # Runs the Miniconda installer script with interactive prompts
conda activate base               # Activates the base Conda environment
conda create -n funtools          # Creates a new Conda environment named funtools
conda activate funtools           # Activates the funtools Conda environment
sudo apt-get update               # Updates the package lists for the system 
sudo apt-get install figlet -y    # Installs the figlet package with automatic yes to prompts 
figlet StMark                     # Runs figlet to display "StMark" in ASCII art
conda install -c bioconda bwa     # Installs the bwa tool from the bioconda channel
conda install -c bioconda blast   # Installs the blast tool from the bioconda channel
conda install -c bioconda samtools  # Installs the samtools tool from the bioconda channel
conda install -c bioconda bedtools  # Installs the bedtools tool from the bioconda channel
conda install -c bioconda spades   # Installs the spades tool from the bioconda channel
conda install -c bioconda bcftools  # Installs the bcftools tool from the bioconda channel
conda install -c bioconda fastp    # Installs the fastp tool from the bioconda channel
conda install -c bioconda multiqc  # Installs the multiqc tool from the bioconda channel 


#Professional Profile
Github repo submission - https://github.com/StMarkFx/HackBioNGS/blob/main/stage0.sh
LinkedIn post - https://www.linkedin.com/posts/stmarkadebayo_bioinformatics-hackbio-genomics-activity-7369740282684182531-Etm6?utm_source=share&utm_medium=member_desktop&rcm=ACoAAEPqdXgB8Q0PUcaOrdpjrOIcdNQeM4UehLA