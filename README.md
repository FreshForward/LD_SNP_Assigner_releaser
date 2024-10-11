# LD_SNP_Assigner

Python code that computes all LD values of a set of variants with itself or another set of variants, to address whether the variant is on the right chromosome and position.
This is particularly useful in allopolyploid organisms that have high sequence similarity among different genomes, and as such sequencing reads can be aligned on the biologically wrong position and eventually result in spurious variant calls. This progam tries to estimate whether this is the case for variants or not.

Please cite the following article when using this tool:
How to handle high subgenome sequence similarity in allopolyploid Fragaria x ananassa: Linkage Disequilibrium Based Variant Filtering: https://www.researchsquare.com/article/rs-4686718/v1

# Installation
```bash
conda create --name LD_assignment python=3.11.2
conda activate LD_assignment
pip install scipy==1.10.1
pip install pandas==1.5.3
pip install numpy==1.24.3
pip install tables==3.8.0
pip install PyYAML==6.0
pip install seaborn==0.12.2
pip install h5py==3.10.0
```

Ensure that you call python correctly: this program calls 'python' when calling python files. 
If your computer uses bash and requires 'python3'; enter:
```bash
alias python="python3"
```

# Input
This program currently requires CSV files as input which need to be defined in the config file. To obtain these quick from .vcf files the +dosage plugin from bcftools can be used in combination with bash cut, awk and sed commands:
```bash
sbatch --wrap="bcftools +dosage /Path/To/Input/input_file.vcf -- -t GT | cut -f 1-3,5- | awk -F'\t' 'BEGIN{OFS=\"\t\"} NR==1 {\$1=\"Marker\"; \$2=\"Chrom\"; \$3=\"Position\"; print} NR>1 {\$3=\$2; \$2=\$1; \$1=\$2\"_\"\$3; print}' | sed -e '1 s/\[[^][]*\]//g' | sed 's/\t/,/g' > /Path/To/LD_assigner/Datafiles/input_file.csv"
```
Note: for above command it is important that Chrom+Position is a unique identifier for each variant.  

# Run pipeline

There are a few steps to take before running the pipeline. These are:
- Ensure your terminal is in the correct LD_assignment folder.
- Ensure the correct settings are applied in the "config_example.yaml" file. 

To run the program locally, use the command:
```python
python3 calculate_chromosome.py config_example.yaml
```

To run a long running program over the weekend or night, use the command:
```python
nohup python3 -u calculate_chromosome.py config_example.yaml &
```

The nohup command ensures that breaking the SSH connection does not stop the program from running, where the -u command ensures the print statements and errors get printed to the "nohup.out" file. 

# Run with screen command so you can exit the terminal

```python
screen
python3 calculate_chromosome.py config_example.yaml
```

Re- connect to screen after it is done by:
```python
screen -ls
```

Find the screen you are interesed in- which ran the process, and type screen with the number code:
```python
screen xxxxxx
exit
```

If the pipeline only finished one script, please ensure 'python' is a valid arg.
Or run the scripts in order manually. 
```python
python calculate_chromosome.py config_example.yaml
python calculate_position.py config_example.yaml
python summaryfile_creation.py config_example.yaml
```

# Run on SLURM server
If using this script on a SLURM server program can be run with:
```bash
sbatch --output=LD_XX_output.txt --error=LD_XX_error.txt --wrap="python3 calculate_chromosome.py config_example.yaml"
```