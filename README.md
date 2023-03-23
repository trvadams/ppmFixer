# ppmFixer
This script takes a pGlyco3.0 output file (e.g. pGlycoDB-GP-FDR-Pro-Quant-Site.txt)
and outputs a separate .txt file that adds additional columns with corrected
N-glycan compositions (not linkages).

Usage: python ppmFixer <input_file> <output_file>

Example: python ppmFixer.py pGlycoDB-GP-FDR-Pro-Quant-Site.txt output.txt

Required packages:
- numpy
- pandas

For more information please see associated publication:
