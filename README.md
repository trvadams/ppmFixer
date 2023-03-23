# ppmFixer
This script takes a pGlyco3.0 output file (e.g. pGlycoDB-GP-FDR-Pro-Quant-Site.txt)
and outputs a separate .txt file that adds additional columns with corrected
glycan compositions (not linkages).

Usage: python ppmFixer <input_file> <output_file>
e.g.: python ppmFixer pGlycoDB-GP-FDR-Pro-Quant-Site.txt output.txt

Required packages:
- numpy
- pandas
For more information please see associated publication:
