# SHAPE-MaP Processing Pipeline
The analysis pipeline for selective 2′-hydroxyl acylation analyzed by primer extension and mutational profiling (SHAPE-MaP) assays, processes raw FASTQ data through sequential steps such as initial basecall quality trimming, paired read merging, genome alignment, mutation handling, calculation of mutation rates, and reactivity profile calculation and normalization, and it supports multiple samples. Also, we provide a fully containerized Singularity environment that bundles all required tools and dependencies, and with a single command, the entire workflow can be executed reproducibly on any compatible system, supporting multiple samples.

# Part I Workflow
Here stands an throughout workflow of data analysis.
<img width="1731" height="655" alt="1 workflow" src="https://github.com/user-attachments/assets/f2de1227-d413-4c6c-9d90-3e769f39cd13" />

# Part II Requirements
1.  **Recommended System Configuration**:

      * 8-core CPU
      * 24 GB RAM

2.  **Singularity**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. For other operating systems, please refer to the official installation guide: [https://docs.sylabs.io/guides/3.0/user-guide/installation.html](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

      * **Step 1: Install System Dependencies**

        ```bash
        # Update package lists and install dependencies
        sudo apt-get update
        sudo apt-get install -y \
            build-essential \
            libseccomp-dev \
			libfuse3-dev \
            pkg-config \
            squashfs-tools \
            cryptsetup \
            curl wget git
        ```

      * **Step 2: Install Go Language**

        ```bash
        # Download and install Go
        wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
        sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
        rm go1.21.3.linux-amd64.tar.gz

        # Configure Go environment variables and apply them
        echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
        echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
        source ~/.bashrc
        ```

      * **Step 3: Download, Build, and Install Singularity**

        ```bash
        # Note: The script navigates to /mnt/share/software. 
        # You can change this to your preferred directory for source code.
        cd /mnt/share/software

        # Download the Singularity CE source code
        wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz

        # Extract the archive and clean up
        tar -xvzf singularity-ce-4.0.1.tar.gz
        rm singularity-ce-4.0.1.tar.gz
        cd singularity-ce-4.0.1

        # Configure the build
        ./mconfig

        # Build Singularity (this can be time-consuming)
        cd builddir
        make

        # Install Singularity to the system
        sudo make install
        ```

      * **Step 4: Verify the Installation**

        ```bash
        # Check the installed version
        singularity --version

        # Display help information
        singularity -h
        ```


3.  **snakemake**: Snakemake must be installed on your system and requires a Python 3 distribution.

      ```bash
      pip install snakemake
      ```

4.  **Data Preparation**: The data run by this pipeline is from SRR30943151 and SRR30943152 in the SRA database.The specific processing method is as follows

    ```bash
    # Download the test sra data
    mkdir -p data/samples
    cd data/samples
    prefetch SRR30943151
    prefetch SRR30943152

    # Convert sra data to fastq data
    fastq-dump --gzip --split-files SRR30943151\SRR30943151.sra
    fastq-dump --gzip --split-files SRR30943152\SRR30943152.sra

    # Download the target fasta file
    cd ../
    wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279192/suppl/GSE279192%5FHs%5FDRAIC.fa.gz
    # Download the primer file that may exists
    wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279192/suppl/GSE279192%5FHs%5FDRAIC%5FAmp%5FPrimers.fa.gz
    ```

5.   **Required File Structure**

      ```bash
        root/
            ├── config.yaml
            ├── data
            │   ├── GSE279192_Hs_DRAIC_Amp_Primers.fa
            │   ├── GSE279192_Hs_DRAIC.fa
            │   └── samples
            │       ├── SRR30943151_1.fastq.gz
            │       ├── SRR30943151_2.fastq.gz
            │       ├── SRR30943152_1.fastq.gz
            │       └── SRR30943152_2.fastq.gz
            ├── shapemap.sif
            ├── split_reference.py
            └── SHAPE-MaP.smk
      ```
      
      - **SHAPE-MaP.smk** — The main Snakemake workflow script.  
      - **config.yaml** — Configuration file containing paths, parameters, and sample information.  
        ⚠️ Must be located in the same directory as `SHAPE-MaP.smk`.
      - **shapemap.sif** — Singularity container image with all required software and dependencies pre-installed.
      - **split_reference.py** — A python script to split target fasta into single-sequence files.

# Part III Running

   * **Example code**

      * **Step 1: Edit `config.yaml`**

        ```bash
        # Please use ABSOLUTE PATH for all file paths
        # Sample configuration file for SHAPE-MaP pipeline
        samples:
          GSE279192: 
            modified:
              R1 : "/mnt/zhangam/SHAPE_MaP/data/samples/SRR30943151_1.fastq.gz"
              R2 : "/mnt/zhangam/SHAPE_MaP/data/samples/SRR30943151_2.fastq.gz"
            untreated:
              R1 : "/mnt/zhangam/SHAPE_MaP/data/samples/SRR30943152_1.fastq.gz"
              R2 : "/mnt/zhangam/SHAPE_MaP/data/samples/SRR30943152_2.fastq.gz"
        # Output directory no "/" at the end
        output_dir: "/mnt/zhangam/SHAPE_MaP/output"
        container: "/mnt/zhangam/SHAPE_MaP/shapemap.sif"

        # Target reference file
        target: "/mnt/zhangam/SHAPE_MaP/data/GSE279192_Hs_DRAIC.fa"

        # Script settings
        scripts:
        	split: "/mnt/zhangam/SHAPE_MaP/split_reference.py"

        # ShapeMapper2 parameters
        denatured: false
        overwrite: true
        indiv_norm: false

        # Primer settings
        amplicon: true
        certain_primer: true
        primers_file: /mnt/zhangam/SHAPE_MaP/data/GSE279192_Hs_DRAIC_Amp_Primers.fa
        random_primer: false
        random_primer_length: 0
        max_primer_offset: 10

        # ShapeMapper2 parameters (delete to set default values)
        min_depth: 1000
        max_bg: 0.05
        min_qual_to_trim: 20
        window_to_trim: 5
        min_length_to_trim: 25
        min_qual_to_count: 30

        # RNAstructure parameters
        temperature: 37
        maximum: 20
        loop: 30
        ```

		Please note the primer settings. If `certain_primer` is used, please set as follows:

	  	```bash
   		# Primer settings
   		amplicon: (set based on facts)
        certain_primer: true
        primers_file: /absolute/path/to/primers/file
        random_primer: false
        random_primer_length: (will not be read)
        max_primer_offset:  (set based on facts)
    	```

		If `random_primer` is used, please set as follows:

	  	```bash
   		# Primer settings
   		amplicon: false
        certain_primer: false
        primers_file: (will not be read)
        random_primer: true
        random_primer_length: (set based on facts)
        max_primer_offset: (will not be read) 
    	```

	  * **Step 2: Dry-run and dag-make**
      Here /mnt/zhangam/SHAPE_MaP/ represents the root directory.

		```bash
    	# Dry-run
        snakemake -np \
          -s SHAPE-MaP.smk \
          --use-singularity \
          --singularity-args "--bind /mnt/zhangam/SHAPE_MaP/"

		# Dag-make
		snakemake -s SHAPE-MaP.smk \
		          --use-singularity \
		          --singularity-args "--bind /mnt/zhangam/SHAPE_MaP/" \
		          --dag  | \
		dot -Tsvg > dag.svg
        ```
		Please try dry-run and dag-make first to check pipeline usability and generate flow diagram.

      * **Step 3: run snakemake**
      Here /mnt/zhangam/SHAPE_MaP/ represents the root directory.

        ```bash
        cd /path/to/snakemake/file/
        snakemake -s SHAPE-MaP.smk \
		      --cores 8 \
          --use-singularity \
          --singularity-args "--bind /mnt/zhangam/SHAPE_MaP/"
        ```
      
   * **Command Parameters**

      **edit `config.yaml`**
      - `samples`:(required) A list describing all input samples, including their names and raw FASTQ file paths. Each sample entry must contain: `sample`: the sample name; `modified`: experimental group, in which RNA molecules are covalently modified at structurally flexible, unpaired nucleotides using SHAPE chemicals; `untreated`: negative control, where no SHAPE chemical reagent is added;`R1`: path to the Read 1 FASTQ file; `R2`: path to the Read 2 FASTQ file (for paired-end data). (optional) `denatured`: positive control, where RNA is fully denatured using chemical denaturants or high temperature prior to or concurrently with the addition of the SHAPE reagent. If your dataset includes a denatured group, please set ShapeMapper2 parameters `denatured` as "true" and add sample information following the format described above.
      
      - `output_dir`:(required) Path to the output directory where all results will be saved.

      - `container`:(required) Path to the Singularity container image (`shapemap.sif`) containing all required software and dependencies.

      - `target`:(required) FASTA file containing the target DNA or RNA sequences. Lowercase positions will be excluded from reactivity profile, and should be used to indicate primer-binding sites if using directed primers. If multiple primer pairs were used, set `certain_primer` as "true" and provide the primer sequences in a separate file with `primers_file`.
    
      - `scripts`:(required) The absolute path of the script "split_reference.py".The single-sequence results will be saved in `output_dir`/reference_split.
        
      - `overwrite`:(optional) Overwrite existing files in output and temporary file folders without warning. Default: false.

      - `indiv_norm`:(optional) Normalize multiple reactivity profiles individually, instead of as a group. Default: false.

      - `amplicon`:(optional) Require reads to align near expected primer pair locations, and intelligently trim primer sites. If multiple pairs or internal locations are needed, set `certain_primer` as "true" and specify primers with a `primers_file`.

      - `primers_file`:(optional) Amplicon primer pair sequences. Each line should contain a pair of primer sequences: the forward primer first followed by the reverse primer, separated by whitespace. To specify primers for multiple RNAs, add a line with each RNA name preceded by '>' before each group of primer pairs. RNA names must match those in any provided .fa files. If a primer file is needed, `certain_primer` must be set to "true"; otherwise, the primer file will not be read.

      - `random_primer_length`:(optional) Length of random primers used (if any). Mutations within (length+1) of the 3-prime end of a read will be excluded, as will read depths over this region. Unused if `amplicon` and/or `certain_primer` are provided. If random primers are used, `random_primer` must be set to "true"; otherwise, the random primer length will not be read. Default: 0.

      - `max_primer_offset`:If `amplicon` and/or `certain_primer` used, require read ends to align to within +/- this many nucleotides of expected amplicon primer pairs. Default: 10.

      - `min_depth`:(optional) Minimum effective sequencing depth for including data (threshold must be met for all provided samples). Default: 5000.
      
      - `max_bg`:(optional) Maximum allowed mutation frequency in untreated sample. Default: 0.05.
      
      - `min_qual_to_trim`:(optional) Minimum phred score in initial basecall quality trimming. Default: 20.
      
      - `window_to_trim`:(optional) Window size in initial basecall quality trimming. Default: 5.
      
      - `min_length_to_trim`:(optional) Minimum trimmed read length in initial basecall quality trimming.Default: 25.
      
      - `min_qual_to_count`:(optional) Only count mutations with all basecall quality scores meeting this minimum score (including the immediate upstream and downstream basecalls). This threshold is also used when computing the effective read depth. Default: 30.

      - `temperature`:(optional) Specify the temperature at which calculation takes place in Kelvin. Default：37 degrees C, which is 310.15 K.
      
      - `maximum`:(optional) Specify a maximum number of structures. Note that suboptimal structures are generated until the maximum number of structures are reached. Default: 20.
      
      - `loop`: (optional) Specify the maximum number of unpaired nucleotides in an internal or bulge loop.Default: 30.

      **run snakemake**
      - `--use-singularity`: Enables execution of rules within a Singularity container to ensure a fully reproducible environment.
      - `--singularity-args`: Allows passing additional arguments to the Singularity runtime (e.g., `--bind`, `--nv`, or custom options).
      - `--cores`: Specifies the maximum number of CPU cores (threads) that Snakemake can use in parallel when executing workflow rules.
      - `--bind`: Specifies the directories to be mounted within the Singularity container. Include all required paths such as raw data, scripts, container images, and references. The format is `/project_directory:/project_directory`. Multiple directories can be mounted by separating them with commas, for example: `/path1:/path1,/path2:/path2`. If you're setting `-- bind /project_directory/`, note that the path in the config file must be set as ABSOLUTE PATH and make sure that all the files you need are included in `/project_directory`. (required)

# Part IV Output

   * **Output Structure**
     ```bash
		output_dir/
		    ├── multiqc_data
		    ├── multiqc_report.html
		    ├── qc
		    │   ├── modified
		    │   │   └── (modified group raw qc)
		    │   └── untreated
		    │       └── (untreated group raw qc)
		    ├── reference_split
		    │   ├── DNA
		    │   │   └── (DNA single-sequence fasta files)
		    │   └── RNA
		    │       └── (RNA single-sequence fasta files)
		    └── GSE279192
		        ├── Hs_DRAIC_ncRNA
		        │   ├── structure
		        │   │   ├── GSE279192_Hs_DRAIC_ncRNA.ct
		        │   │   ├── GSE279192_Hs_DRAIC_ncRNA.dbn
		        │   │   └── GSE279192_Hs_DRAIC_ncRNA_folding.svg
		        │   ├── GSE279192_Hs_DRAIC_ncRNA_histograms.pdf
		        │   ├── GSE279192_Hs_DRAIC_ncRNA.map
		        │   ├── GSE279192_Hs_DRAIC_ncRNA_mapped_depths.pdf
		        │   ├── GSE279192_Hs_DRAIC_ncRNA_per-amplicon_abundance.txt
		        │   ├── GSE279192_Hs_DRAIC_ncRNA_profiles.pdf
		        │   ├── GSE279192_Hs_DRAIC_ncRNA_profile.txt
		        │   ├── GSE279192_Hs_DRAIC_ncRNA_ribosketch_colors.txt
		        │   ├── GSE279192_Hs_DRAIC_ncRNA.shape
		        │   ├── GSE279192_Hs_DRAIC_ncRNA_varna_colors.txt
		        │   ├── GSE279192_Modified_Hs_DRAIC_ncRNA_mutation_counts.txt
		        │   ├── GSE279192_Modified_Hs_DRAIC_ncRNA_parsed.mut
		        │   ├── GSE279192_Hs_DRAIC_ncRNA_shapemapper_log.txt
		        │   ├── GSE279192_Untreated_Hs_DRAIC_ncRNA_mutation_counts.txt
		        │   └── GSE279192_Untreated_Hs_DRAIC_ncRNA_parsed.mut
		        └── (Other target sequences, if any exist)
      ```
    
   * **Output Interpretation**
     
	  - **`multiqc_report.html`**: Open multiqc_report.html in a web browser to explore all sections interactively.

        - **General Statistics**: A combined table summarizing important metrics for each sample.

			

        - **FastQC**: Quality-control metrics on raw and trimmed reads.
     
     	

      - **`*_shapemapper_log.txt`**
        
		<img width="642" height="751" alt="2 log" src="https://github.com/user-attachments/assets/e9603eb4-2cca-4d98-ab5d-d5cb556df11d" />

        - **Content**: Run progress and summary outputs. Includes mate pair merging stats, read alignment stats, reactivity profile quality control checks, and amplicon primer pair read depths.
        - **Application**: This is the first file you should check to assess the overall quality of your data.

      - **`*_mapped_depths.pdf`**
        
		<img width="1194" height="654" alt="3 mapped depth" src="https://github.com/user-attachments/assets/ce5f972d-45d4-48f2-8cb6-41410986c680" />

        - **Content**: Figures showing simple mapped read depths. Shows reads excluded due to low aligner-reported MAPQ (mapping quality score), and shows off-target reads excluded due to not aligning near expected amplicon primer pair locations. 
        - **Application**: Reads included in analysis are further broken down by primer pair.

      - **`*_profile.txt`**

        - **Content**: The description of the analysis results including read depths, mutation rates, reactivity profile and so on. The detailed format and contents are as follows:

        | Column name               |  Content                                          |
        |---------------------------|---------------------------------------------------|
        |`Nucleotide`               |  Nucleotide number (1-based)                      |
        |`Sequence`                 |  Nucleotide (AUGCaugc)                            |
        |`<Sample>_mutations`       |  Mutation counts                                  |
        |`<Sample>_read_depth`      |  Read depth                                       |
        |`<Sample>_effective_depth` |  see [Effective read depth](analysis_steps.md#effective-read-depth) |
        |`<Sample>_rate`            |  Effective mutation rate calculated as <br> <tt>(mutation count / effective read depth)</tt> |
        |`<Sample>_off_target_mapped_depth` | Simple mapped read depths for reads <br> not meeting <kbd>--amplicon</kbd>/<kbd>--primer</kbd> location <br> requirements |
        |`<Sample>_low_mapq_mapped_depth` | Simple mapped read depths for reads <br> not meeting <kbd>--min-mapq</kbd>  |
        |`<Sample>_mapped_depth` <br> or `<Sample>_primer_pair_<n>_mapped_depth` | Simple mapped read depths for included <br> reads, broken down by primer pair if <br> applicable |
        |`Reactivity_profile`       |  Calculated reactivity profile                    |
        |`Std_err`                  |  Standard error                                   |
        |`HQ_profile`               |  Reactivity profile with high-background and <br> low-depth positions excluded (set to nan)     |
        |`HQ_stderr`                |  Standard error with high-background and <br> low-depth positions excluded |
        |`Norm_profile`             |  Reactivity profile after normalization <br> (scaling) | 
        |`Norm_stderr`              |  Standard error after normalization               |

        - **Application**: Show the analysis results, where visualizations of read depths, mutation rates, and reactivity profile can be viewed at `*_profiles.pdf` and `*_histograms.pdf`.
       
          <img width="1722" height="234" alt="4 profile" src="https://github.com/user-attachments/assets/44d38d16-f1e1-416b-9580-e0f0d3ff27d7" />

		  <img width="1131" height="680" alt="5 profile histogram" src="https://github.com/user-attachments/assets/f5c2e342-57da-451b-8fcb-213509a76903" />

      - **`*_mutation_counts.txt`**

        - **Content**: This file reports mutation counts and sequencing depth metrics for each nucleotide position in the reference sequence.    Columns include counts for each mutation type, listed 5' to 3', along with depth metrics such as effective_depth, depths of filtered reads, and the total mapped depth of analyzed reads.
        - **Application**: Used for subsequent analysis of mutation patterns and calculation of mutation rates. The provided depth metrics are critical for evaluating data quality and ensuring the reliability of the results.

      - **`*_parsed.mut`**

        - **Content**: This is a log where each line represents a single mapped read. It catalogs essential alignment information, including read type, mapping coordinates, and classification (INCLUDED, LOW_MAPQ, OFF_TARGET). It encodes per-read nucleotide-level data through binary arrays for mapped coverage, effective depth, and identified mutation sites. Detailed mutation information and file format are as follows：

	        | Field                     |  Content                                          |
	        |---------------------------|---------------------------------------------------|
	        | 1 | read type (see below) |
	        | 2 | read name |
	        | 3 | 0-based leftmost mapping position (inclusive) |
	        | 4 | 0-based rightmost mapping position (inclusive) |
	        | 5 | read mapping category (see below)|
	        | 6 | primer pair index (0-based), <br> or -999 if no associated primers |
	        | 7 | mapped depth array (see below) |
	        | 8 | effective depth array (see below) |
	        | 9 | mutation count array (see below) |
	        | 10 | mutations (see below) |

        - **Application**: This file serves as the primary data source for calculating mutation rates and profiling mutation patterns.  The detailed mutation records enable advanced analysis of mutation spectra, while the mapping categories facilitate quality control by distinguishing reads used in analysis from those excluded due to low quality or off-target alignment.

      - **`*.shape` and `*.map`**

        - **Content**: These files contain the final nucleotide-resolution SHAPE reactivity profiles, which quantify RNA structure flexibility. The shape file provides the essential data in a standard two-column format (position, reactivity), with unreactive positions marked as -999. The map file contains two additional columns: the standard error of the reactivity measurement and the corresponding nucleotide sequence.
        - **Application**: These files are the direct input for RNAstructure. The reactivity values serve as experimental constraints to guide and validate the computational folding algorithms, enabling more accurate determination of the RNA secondary structure.

      - **`*_varna_colors.txt` and `*_ribosketch_colors.txt`**

        - **Content**: These files provide a simplified, normalized reactivity profile optimized for visualization software.  Each file contains a single column of normalized reactivity values: values exceeding 0.85 are capped at 0.85, negative values are set to 0, and positions with missing data are assigned a value of 0.
        - **Application**: The files serve as direct inputs for VARNA and RiboSketch. The adjusted reactivity scale ensures consistent and interpretable color mapping or styling of nucleotides based on their structural flexibility in the resulting secondary structure diagrams.

      - **`*.ct` and `*.dbn`** : Open multiqc_report.html in a web browser to explore all sections interactively.

        - **Content**: CT and Dot Bracket (.dbn) files encapsulate an RNA sequence and its predicted secondary structure. The CT file uses a table format, detailing the nucleotide sequence and the pairing partner for each position. The Dot Bracket format represents the structure more compactly using a string of symbols, where matching brackets denote base pairs.
        - **Application**: Served as the direct input for visualization tools that generate structure diagrams. The `*.svg` file can be viewed to see the graphical representation of the RNA secondary structure.​

		![6 structure](https://github.com/user-attachments/assets/3e393cb8-85a1-4c20-ad9d-e6462a47cc4e)
	  
# Part Ⅴ Reference
[1] Busan S, Weeks KM. Accurate detection of chemical modifications in RNA by mutational profiling (MaP) with ShapeMapper 2. RNA. 2018, 24(2):143-148.

[2] Reuter, J.S., Mathews, D.H. RNAstructure: software for RNA secondary structure prediction and analysis. BMC Bioinformatics 11, 129 (2010).
