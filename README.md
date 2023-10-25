# ppmFixer
This script takes a pGlyco3.0 output file (e.g. pGlycoDB-GP-FDR-Pro-Quant-Site.txt)
and outputs a separate .txt file that adds additional columns with corrected
N-glycan compositions (not linkages).

Usage: python ppmFixer <input_file> <output_file>

Example: python ppmFixer.py pGlycoDB-GP-FDR-Pro-Quant-Site.txt output.txt

Required packages:
- numpy
- pandas

To install the required packages, I recommend the use of miniconda:
https://docs.conda.io/projects/miniconda/en/latest/

Once miniconda is installed and activated, use the following commands to install the needed packages:
```
conda install numpy
conda install pandas
```
For more information please see associated publication:
