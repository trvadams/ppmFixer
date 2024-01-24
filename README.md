# ppmFixer
This script takes a pGlyco3.0 output file (e.g. pGlycoDB-GP-FDR-Pro-Quant-Site.txt)
and outputs a separate .txt file that adds additional columns with corrected
N-glycan compositions (not linkages).

Usage: 

```
python ppmFixer <input_file> <output_file>
```

Example: 
```
python ppmFixer.py pGlycoDB-GP-FDR-Pro-Quant-Site.txt output.txt
```

Required packages:
- numpy
- pandas

To install the required packages, I recommend the use of miniconda:
https://docs.conda.io/projects/miniconda/en/latest/

Once miniconda is installed and activated, use the following commands to install the required packages:
```
conda install numpy
conda install pandas
```
For more information please see associated publication:
https://doi.org/10.1093/glycob/cwae006

FOR WINDOWS USERS:
Using Python on a Windows machine can be a little tricky, so here are additional instructions for Windows users who might be having trouble installing the proper packages.

1. Download and install Python 3.x (latest version) from here: https://www.python.org/downloads/
2. Download and install miniconda from here: https://docs.conda.io/projects/miniconda/en/latest/index.html
3. Open the newly installed Anaconda Powershell Prompt (Miniconda3).
4. Install the numpy and pandas packages via the powershell. To do this, simply enter:
	conda install numpy
	conda install pandas
5. Use the powershell to run the script. The usage goes:
	python ppmFixer.py <input_file> <output_file>
	with your input file being the "pGlycoDB-GP-FDR-Pro-Quant-Site.txt” output from pGlyco, and the output file being whatever you want (e.g. output.txt).
	Note: If you’re unfamiliar with navigating directories in powershell (which uses a bash-like interface), this link has some tips.
	https://www.sharepointdiary.com/2021/04/change-directory-in-powershell.html
