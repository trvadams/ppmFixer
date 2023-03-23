"""pGlyco PPM Fixer

This script takes a pGlyco output file (e.g. pGlycoDB-GP-FDR-Pro-Quant-Site.txt)
and outputs a separate .txt file that adds additional columns with corrected
glycan compositions (not linkages).
For more information please see associated publication:

"""
# Import required packages
import sys
import numpy as np
import pandas as pd

# Get command line arguments
input_file_path = str(sys.argv[1])
output_file_path = str(sys.argv[2])

## Specify input and output file paths
#input_file_path = 'input_file.txt'
#output_file_path = 'output_file.csv'

# Read in data file
df = pd.read_csv(input_file_path, delimiter='\t')

# Make useful dataset-wide variables
ppm_mean = df['PPM'].mean()
ppm_stdDev = df['PPM'].std()

hexoseMass = 162.0528
hexnacMass = 203.0794
fucoseMass = 146.0579
neu5acMass = 291.0954
protonMass = 1.0073

numSugars = ['numHexose', 'numHexNAc', 'numSialic', 'numFucose']

## Add column that removes spaces from peptides
#df['Peptide_noSpace'] = df['Peptide'].str.replace(" ", "")

# Add number of each type of sugar.
df[['numHexose', 'numHexNAc', 'numSialic', 'numFucose']] = \
    df['Glycan(H,N,A,F)'].str.split(' ', expand=True)

# Convert number of sugars to numeric for math.
for x in numSugars:
    df[x] = pd.to_numeric(df[x])

# Identify whether GPSMs are less (Pass) or greater (fail)
# than one SD away from the mean.
df['ppmCheck'] = np.where(
    df['PPM'] < (ppm_mean - ppm_stdDev),
    'Fail', 'Pass'
)

# Flag GPSMs that fail the ppmCheck and have multiple fucoses.
df['ppmFlag'] = np.where(
    (df['ppmCheck'] == "Fail") & (df['numFucose'] >= 2),
    True, False
)

# Where GPSMs are flagged, recalculate mass error as if 
# two fucoses are converted to one sialic acid.
df['ppmFixPrecursorMZ'] = np.where(df['ppmFlag'] == True, \
    (df['PeptideMH'] + \
        (df['numHexose'] * hexoseMass) + \
        (df['numHexNAc'] * hexnacMass) + \
        ((df['numSialic'] + 1) * neu5acMass) + \
        ((df['numFucose'] - 2) * fucoseMass) + \
        (df['Charge'] * protonMass)) / \
        df['Charge'],
    0
)

# Calculate new ppm error, otherwise set it to 100 
# (high enough where it will always fail the ppmFixSuccess check).
df['ppmFixError'] = np.where(df['ppmFixPrecursorMZ'] != 0,
    (df['PrecursorMZ'] - df['ppmFixPrecursorMZ'])/df['ppmFixPrecursorMZ']*1000000,
    100)

# Flag GPSMs where the new mass error is less than the original mass error.
#df['ppmFixSuccess'] = np.where(abs(df['ppmFixError']) < abs(df['PPM']), True, False)
df['ppmFixSuccess'] = np.where(
    (abs(df['ppmFixError'] - ppm_mean)) < abs(df['PPM'] - ppm_mean),
    True, False
)

# Create new columns with adjusted fucoses/sialic acids
for x in numSugars:
    df[x + '_ppmFixer'] = pd.to_numeric(df[x])
df['numFucose_ppmFixer'] = df['numFucose_ppmFixer'] - 2
df['numSialic_ppmFixer'] = df['numSialic_ppmFixer'] + 1

# Turn ppmFixer columns into strings.
for x in ['numHexose_ppmFixer', 'numHexNAc_ppmFixer', 'numSialic_ppmFixer', 'numFucose_ppmFixer']:
    df[x] = df[x].astype('str')


# Make new column 'Glycan(H,N,A,F)_corrected_ppmFixer' where either the initial composition, built-in
# pGlyco correction, or the ppmFixer correction is used (ppmFixSuccess = True).
# Where both corrections are applicable, the built-in pGlyco correction takes precedence.

# Assign original assignments to column.
df['Glycan(H,N,A,F)_ppmFixer'] = df['Glycan(H,N,A,F)']


## Overwrite where ppmFixSuccess = True
df['Glycan(H,N,A,F)_ppmFixer'].mask(
    df['ppmFixSuccess'] == True,
    df['numHexose_ppmFixer'] + ' ' + df['numHexNAc_ppmFixer'] + ' ' \
        + df['numSialic_ppmFixer'] + ' ' + df['numFucose_ppmFixer'],
        inplace=True
)

# Overwrite where built-in pGlyco correction exists
df['Glycan(H,N,A,F)_ppmFixer'].mask(
    df['CorrectedGlycan(H,N,A,F)'].notna(),
    df['CorrectedGlycan(H,N,A,F)'],
    inplace=True
)

# Make new columns to store these final values
df[['numHexose_final', 'numHexNAc_final', 'numSialic_final', 'numFucose_final']] = \
    df['Glycan(H,N,A,F)_ppmFixer'].str.split(' ', expand=True)

# Make new column 'GlycanComposition_corrected_ppmFixer',
# Apologies for this being pretty confusing, couldn't think of a better way to do it.
df['GlycanComposition_ppmFixer'] = np.where(
    df['numHexose_final'].astype('int') > 0,
    'H(' + df['numHexose_final'] + ')' + np.where(
        df['numHexNAc_final'].astype('int') > 0,
        'N(' + df['numHexNAc_final'] + ')' + np.where(
            df['numSialic_final'].astype('int') > 0,
            'A(' + df['numSialic_final'] + ')' + np.where(
                df['numFucose_final'].astype('int') > 0,
                'F(' + df['numFucose_final'] + ')',
                ''),
            np.where(df['numFucose_final'].astype('int') > 0,
                     'F(' + df['numFucose_final'] + ')',
                     '')),
        ''),
    ''
    )

# Remove all generated columns except for the final _corrected_ppmFixer columns (for tidyness)
df.drop(columns=['numHexose', 'numHexNAc', 'numSialic', 'numFucose'],
        inplace=True)
df.drop(columns=['numHexose_ppmFixer', 'numHexNAc_ppmFixer', 'numSialic_ppmFixer', 'numFucose_ppmFixer'],
        inplace=True) 
df.drop(columns=['numHexose_final', 'numHexNAc_final', 'numSialic_final', 'numFucose_final'],
        inplace=True)


# Save resulting dataframe to tsv (same as original input)
df.to_csv(output_file_path, sep = '\t')