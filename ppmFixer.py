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

# Specify data types for each column
dataTypes = {'GlySpec': str,
             'PepSpec': str,
             'Scan': float,
             'RT': float,
             'PrecursorMH': float,
             'PrecursorMZ': float,
             'Charge': int,
             'Rank': int,
             'Peptide': str,
             'Mod': str,
             'PeptideMH': float,
             'Glycan(H,N,A,F)': str,
             'GlycanComposition': str,
             'PlausibleStruct': str,
             'GlyID': int,
             'GlyFrag': str,
             'GlyMass': float,
             'GlySite': int,
             'TotalScore': float,
             'PepScore': float,
             'GlyScore': float,
             'CoreMatched': int,
             'MassDeviation': float,
             'PPM': float,
             'GlyIonRatio': float,
             'byIonRatio': float,
             'czIonRatio': float,
             'GlyDecoy': float,
             'PepDecoy': float,
             'IsSmallGlycan': int,
             'GlycanPEP': float,
             'GlycanFDR': float,
             'PeptidePEP': float,
             'PeptideFDR': float,
             'TotalFDR': float,
             'Proteins': str,
             'Genes': str,
             'ProSites': str,
             'MonoArea': float,
             'IsotopeArea': float,
             'CorrectedGlycan(H,N,A,F)': str,
             'CorrectedComposition': str,
             'CorrectedMonoArea': float,
             'CorrectedIsotopeArea': float,
             'ETDScan': int,
             'LocalizedSiteGroups': float,
             'LocalizedScore': float,
             'LocalizedIonRatio': float,
             'PreLocalizedScore': float
             }

## Specify input and output file paths (for dummy runs)
#input_file_path = 'input_file.txt'
#output_file_path = 'output_file.csv'

print('Beginning ppmFixer...')

# Read in data file
df = pd.read_csv(input_file_path,
                 delimiter='\t',
                 dtype=dataTypes)

# Create useful dataset-wide constants
ppm_mean = df['PPM'].mean()
ppm_stdDev = df['PPM'].std()

print('Mean mass error is ' + str(ppm_mean)[:5] + ' PPM.')
print('Mass error standard deviation is ' + str(ppm_stdDev)[:5] + ' PPM.')

hexoseMass = 162.052824
hexnacMass = 203.079373
fucoseMass = 146.057909
neu5acMass = 291.095417
protonMass = 1.007276
c13_deltaMass = 1.00335
numSugars = ['numHexose', 'numHexNAc', 'numSialic', 'numFucose']

## Add column that removes spaces from peptides
#df['Peptide_noSpace'] = df['Peptide'].str.replace(" ", "")

# Create some columns used for storing values
df['ppmFixerType'] = 'none'
df['ppmFixPrecursorMZ_temp'] = 0
df['Glycan(H,N,A,F)_ppmFixer'] = df['Glycan(H,N,A,F)']

# Create column to store newly computed mass errors
df['massErrorPPM_ppmfixer'] = df['PPM']

# Overwrite where built-in pGlyco correction exists, note that this took place.
df['Glycan(H,N,A,F)_ppmFixer'].mask(
    df['CorrectedGlycan(H,N,A,F)'].notna(),
    df['CorrectedGlycan(H,N,A,F)'],
    inplace=True
)

df['ppmFixerType'].mask(
    df['CorrectedGlycan(H,N,A,F)'].notna(),
    'native pGlyco correction',
    inplace=True
)

# Add number of each type of sugar.
df[['numHexose', 'numHexNAc', 'numSialic', 'numFucose']] = \
    df['Glycan(H,N,A,F)'].str.split(' ', expand=True)

# Convert number of sugars to numeric for math.
for x in numSugars:
    df[x] = pd.to_numeric(df[x])

# Calculate new mass error where native pGlyco correction took place and save value to temporary column.
df['ppmFixPrecursorMZ_temp'] = np.where(df['CorrectedGlycan(H,N,A,F)'].notna(), \
    (df['PeptideMH'] - protonMass + \
        (df['numHexose'] * hexoseMass) + \
        (df['numHexNAc'] * hexnacMass) + \
        ((df['numSialic'] + 1) * neu5acMass) + \
        ((df['numFucose'] - 2) * fucoseMass) + \
        (df['Charge'] * protonMass)) / \
        df['Charge'],
    0
)

# Calculate new temporary ppm error, otherwise set it to 100 
# (high enough where it will always fail the ppmFixSuccess check).
# For this correction, we are moving down an isotope for the theoretical mass (so subtract this value: (1*protonMass/charge))
df['ppmFixError_temp'] = np.where(df['ppmFixPrecursorMZ_temp'] != 0,
    ((df['PrecursorMZ']-(c13_deltaMass/df['Charge'])) - df['ppmFixPrecursorMZ_temp'])/df['ppmFixPrecursorMZ_temp']*1000000,
    100)

# Set new mass error for native pGlyco correction
df['massErrorPPM_ppmfixer'].mask(
        df['CorrectedGlycan(H,N,A,F)'].notna() == True,
        df['ppmFixError_temp'],
        inplace=True
    )

#########################################
###        Function Definitions       ###
#########################################

# Identify whether GPSMs are less (Pass) or greater (fail)
# than one SD away from the mean.
def ppmCheck():
    df['ppmCheck'] = np.where(
        (df['PPM'] < (ppm_mean - ppm_stdDev)) | (df['PPM'] > (ppm_mean + ppm_stdDev)),
        'Fail', 'Pass'
    )


# Correct double fucoses to single sialic acids.
def ppmFix_2F_A1():
    # Flag GPSMs that fail the ppmCheck and have multiple fucoses.
    df['ppmFucoseFlag'] = np.where(
        (df['ppmCheck'] == "Fail") & (df['numFucose'] >= 2),
        True, False
    )

    # Where GPSMs are flagged, recalculate mass error as if 
    # two fucoses are converted to one sialic acid, store in temporary column.
    df['ppmFixPrecursorMZ_temp'] = np.where(df['ppmFucoseFlag'] == True, \
        (df['PeptideMH'] - protonMass + \
            (df['numHexose'] * hexoseMass) + \
            (df['numHexNAc'] * hexnacMass) + \
            ((df['numSialic'] + 1) * neu5acMass) + \
            ((df['numFucose'] - 2) * fucoseMass) + \
            (df['Charge'] * protonMass)) / \
            df['Charge'],
        0
    )

    # Calculate new temporary ppm error, otherwise set it to 100 
    # (high enough where it will always fail the ppmFixSuccess check).
    # For this correction, we are moving down an isotope for the theoretical mass (so subtract this value: (1*protonMass/charge))
    df['ppmFixError_temp'] = np.where(df['ppmFixPrecursorMZ_temp'] != 0,
        ((df['PrecursorMZ']-(c13_deltaMass/df['Charge'])) - df['ppmFixPrecursorMZ_temp'])/df['ppmFixPrecursorMZ_temp']*1000000,
        100)

    # Flag GPSMs where the new mass error is less than the original mass error.
    #df['ppmFixSuccess'] = np.where(abs(df['ppmFixError']) < abs(df['PPM']), True, False)
    df['ppmFixSuccess'] = np.where(
        (((abs(df['ppmFixError_temp'] - ppm_mean)) < abs(df['PPM'] - ppm_mean)) & \
            (abs(df['ppmFixError_temp'] - ppm_mean) < abs(df['massErrorPPM_ppmfixer'] - ppm_mean))),
        True, False
    )


    # Create new temporary columns with adjusted fucoses/sialic acids
    for x in numSugars:
        df[x + '_ppmFixer'] = pd.to_numeric(df[x])
    df['numFucose_ppmFixer'] = df['numFucose_ppmFixer'] - 2
    df['numSialic_ppmFixer'] = df['numSialic_ppmFixer'] + 1

    # Turn ppmFixer columns into strings.
    for x in ['numHexose_ppmFixer', 'numHexNAc_ppmFixer', 'numSialic_ppmFixer', 'numFucose_ppmFixer']:
        df[x] = df[x].astype('str')

    # Overwrite glycans where ppmFixSuccess = True (native pGlyco correction takes precedence)
    df['Glycan(H,N,A,F)_ppmFixer'].mask(
        (df['ppmFixSuccess'] == True) & (df['CorrectedGlycan(H,N,A,F)'].notna() == False),
        df['numHexose_ppmFixer'] + ' ' + df['numHexNAc_ppmFixer'] + ' ' \
            + df['numSialic_ppmFixer'] + ' ' + df['numFucose_ppmFixer'],
            inplace=True
    )

    # Note the type of correction made
    df['ppmFixerType'].mask(
        (df['ppmFixSuccess'] == True) & (df['CorrectedGlycan(H,N,A,F)'].notna() == False),
        'F2 -> A1',
        inplace=True
    )

    # Where changes are made, save new mass error.
    df['massErrorPPM_ppmfixer'].mask(
        (df['ppmFixSuccess'] == True) & (df['CorrectedGlycan(H,N,A,F)'].notna() == False),
        df['ppmFixError_temp'],
        inplace=True
    )

    # Clean up table
    df.drop(columns='ppmFucoseFlag',
            inplace=True)       

# Correct quadruple fucoses to double sialic acids.
def ppmFix_4F_A2():
    # Flag GPSMs that fail the ppmCheck and have multiple fucoses.
    df['ppmFucoseFlag'] = np.where(
        (df['ppmCheck'] == "Fail") & (df['numFucose'] >= 4),
        True, False
    )

    # Where GPSMs are flagged, recalculate mass error as if 
    # two fucoses are converted to one sialic acid, store in temporary column.
    df['ppmFixPrecursorMZ_temp'] = np.where(df['ppmFucoseFlag'] == True, \
        (df['PeptideMH'] - protonMass + \
            (df['numHexose'] * hexoseMass) + \
            (df['numHexNAc'] * hexnacMass) + \
            ((df['numSialic'] + 2) * neu5acMass) + \
            ((df['numFucose'] - 4) * fucoseMass) + \
            (df['Charge'] * protonMass)) / \
            df['Charge'],
        0
    )

    # Calculate new temporary ppm error, otherwise set it to 100 
    # (high enough where it will always fail the ppmFixSuccess check).
    # For this correction, we are moving down an isotope for the theoretical mass (so subtract this value: (1*protonMass/charge))
    df['ppmFixError_temp'] = np.where(df['ppmFixPrecursorMZ_temp'] != 0,
        ((df['PrecursorMZ']-(2*c13_deltaMass/df['Charge'])) - df['ppmFixPrecursorMZ_temp'])/df['ppmFixPrecursorMZ_temp']*1000000,
        100)

    # Flag GPSMs where the new mass error is less than the original mass error.
    #df['ppmFixSuccess'] = np.where(abs(df['ppmFixError']) < abs(df['PPM']), True, False)
    df['ppmFixSuccess'] = np.where(
        (((abs(df['ppmFixError_temp'] - ppm_mean)) < abs(df['PPM'] - ppm_mean)) & \
            (abs(df['ppmFixError_temp'] - ppm_mean) < abs(df['massErrorPPM_ppmfixer'] - ppm_mean))),
        True, False
    )


    # Create new temporary columns with adjusted fucoses/sialic acids
    for x in numSugars:
        df[x + '_ppmFixer'] = pd.to_numeric(df[x])
    df['numFucose_ppmFixer'] = df['numFucose_ppmFixer'] - 4
    df['numSialic_ppmFixer'] = df['numSialic_ppmFixer'] + 2

    # Turn ppmFixer columns into strings.
    for x in ['numHexose_ppmFixer', 'numHexNAc_ppmFixer', 'numSialic_ppmFixer', 'numFucose_ppmFixer']:
        df[x] = df[x].astype('str')

    # Overwrite glycans where ppmFixSuccess = True (native pGlyco correction takes precedence)
    df['Glycan(H,N,A,F)_ppmFixer'].mask(
        (df['ppmFixSuccess'] == True) & (df['CorrectedGlycan(H,N,A,F)'].notna() == False),
        df['numHexose_ppmFixer'] + ' ' + df['numHexNAc_ppmFixer'] + ' ' \
            + df['numSialic_ppmFixer'] + ' ' + df['numFucose_ppmFixer'],
            inplace=True
    )

    # Note the type of correction made
    df['ppmFixerType'].mask(
        (df['ppmFixSuccess'] == True) & (df['CorrectedGlycan(H,N,A,F)'].notna() == False),
        'F4 -> A2',
        inplace=True
    )

    # Where changes are made, save new mass error.
    df['massErrorPPM_ppmfixer'].mask(
        (df['ppmFixSuccess'] == True) & (df['CorrectedGlycan(H,N,A,F)'].notna() == False),
        df['ppmFixError_temp'],
        inplace=True
    )

    # Clean up table
    df.drop(columns='ppmFucoseFlag',
            inplace=True)       

# Correct N4HxFxAx -> N2H(+7)
def ppmFix_N4HxFxAx_N2Hp7():
    # Flag GPSMs that fail the ppmCheck and have multiple fucoses.
    df['ppmN4Hx+7Flag'] = np.where(
        (df['ppmCheck'] == "Fail") & (df['numHexNAc'] == 4),
        True, False
    )

    # Where GPSMs are flagged, recalculate mass error as if 
    # a complex structure is converted to high mannose, store in temporary column.
    df['ppmFixPrecursorMZ_temp'] = np.where(df['ppmN4Hx+7Flag'] == True, \
        (df['PeptideMH'] - protonMass + \
            ((df['numHexose'] + 7) * hexoseMass) + \
            ((df['numHexNAc'] - 2) * hexnacMass) + \
            ((df['numSialic'] - df['numSialic']) * neu5acMass) + \
            ((df['numFucose'] - df['numFucose']) * fucoseMass) + \
            (df['Charge'] * protonMass)) / \
            df['Charge'],
        0
    )

    # Calculate new temporary ppm error, otherwise set it to 100 
    # (high enough where it will always fail the ppmFixSuccess check).
    # For this correction, we are moving down an isotope for the theoretical mass (so subtract this value: (1*protonMass/charge))
    df['ppmFixError_temp'] = np.where(df['ppmFixPrecursorMZ_temp'] != 0,
        ((df['PrecursorMZ']-(c13_deltaMass/df['Charge'])) - df['ppmFixPrecursorMZ_temp'])/df['ppmFixPrecursorMZ_temp']*1000000,
        100)

    # Flag GPSMs where the new mass error is less than the original mass error.
    # Additionally, check against other current corrections and only flag as success if error is less than previous correction.
    #df['ppmFixSuccess'] = np.where(abs(df['ppmFixError']) < abs(df['PPM']), True, False)
    df['ppmFixSuccess'] = np.where(
        (((abs(df['ppmFixError_temp'] - ppm_mean)) < abs(df['PPM'] - ppm_mean)) & \
            (abs(df['ppmFixError_temp'] - ppm_mean) < abs(df['massErrorPPM_ppmfixer'] - ppm_mean))),
        True, False
    )


    # Create new temporary columns with adjusted fucoses/sialic acids
    for x in numSugars:
        df[x + '_ppmFixer'] = pd.to_numeric(df[x])
    df['numHexose_ppmFixer'] = df['numHexose_ppmFixer'] + 7
    df['numHexNAc_ppmFixer'] = df['numHexNAc_ppmFixer'] - 2
    df['numFucose_ppmFixer'] = df['numFucose_ppmFixer'] - df['numFucose_ppmFixer']
    df['numSialic_ppmFixer'] = df['numSialic_ppmFixer'] - df['numSialic_ppmFixer']

    # Turn ppmFixer columns into strings.
    for x in ['numHexose_ppmFixer', 'numHexNAc_ppmFixer', 'numSialic_ppmFixer', 'numFucose_ppmFixer']:
        df[x] = df[x].astype('str')

    # Overwrite glycans where ppmFixSuccess = True (native pGlyco correction takes precedence)
    df['Glycan(H,N,A,F)_ppmFixer'].mask(
        (df['ppmFixSuccess'] == True) & (df['CorrectedGlycan(H,N,A,F)'].notna() == False),
        df['numHexose_ppmFixer'] + ' ' + df['numHexNAc_ppmFixer'] + ' ' \
            + df['numSialic_ppmFixer'] + ' ' + df['numFucose_ppmFixer'],
            inplace=True
    )

    # Note the type of correction made
    df['ppmFixerType'].mask(
        (df['ppmFixSuccess'] == True) & (df['CorrectedGlycan(H,N,A,F)'].notna() == False),
        'N4HxAxFx -> N2H(x+7)',
        inplace=True
    )

    # Where changes are made, save new mass error.
    df['massErrorPPM_ppmfixer'].mask(
        (df['ppmFixSuccess'] == True) & (df['CorrectedGlycan(H,N,A,F)'].notna() == False),
        df['ppmFixError_temp'],
        inplace=True
    )

    # Clean up table
    df.drop(columns='ppmN4Hx+7Flag',
            inplace=True)       

# Correct NxHxAx -> N(x-4)H(x+5)
def ppmFix_NxHx_Nxm4Hxp5():
    # Flag GPSMs that fail the ppmCheck and have composition of N3H6.
    df['ppmN6HxFlag'] = np.where(
        (df['ppmCheck'] == "Fail") & (5 < df['numHexNAc']),
        True, False
    )

    # Where GPSMs are flagged, recalculate mass error as if 
    # Hexoses and HexNAcs are adjusted.
    df['ppmFixPrecursorMZ_temp'] = np.where(df['ppmN6HxFlag'] == True, \
        (df['PeptideMH'] - protonMass + \
            ((df['numHexose'] + 5) * hexoseMass) + \
            ((df['numHexNAc'] - 4) * hexnacMass) + \
            (df['numSialic'] * neu5acMass) + \
            (df['numFucose'] * fucoseMass) + \
            (df['Charge'] * protonMass)) / \
            df['Charge'],
        0
    )

    # Calculate new temporary ppm error for GPSMs with fitting criteria, otherwise set it to 100 
    # (high enough where it will always fail the ppmFixSuccess check).
    # For this correction, we are moving down two isotopes for the theoretical mass (so subtract this value: (2*protonMass/charge))
    df['ppmFixError_temp'] = np.where(df['ppmFixPrecursorMZ_temp'] != 0,
        ((df['PrecursorMZ']-(2*c13_deltaMass/df['Charge'])) - df['ppmFixPrecursorMZ_temp'])/df['ppmFixPrecursorMZ_temp']*1000000,
        100)

    # Flag GPSMs where the new mass error is less than the original mass error.
    #df['ppmFixSuccess'] = np.where(abs(df['ppmFixError']) < abs(df['PPM']), True, False)
    df['ppmFixSuccess'] = np.where(
        (abs(df['ppmFixError_temp'] - ppm_mean)) < abs(df['PPM'] - ppm_mean),
        True, False
    )


    # Create new temporary columns with adjusted Hexoses/HexNAcs
    for x in numSugars:
        df[x + '_ppmFixer'] = pd.to_numeric(df[x])
    df['numHexose_ppmFixer'] = df['numHexose_ppmFixer'] + 5
    df['numHexNAc_ppmFixer'] = df['numHexNAc_ppmFixer'] - 4

    # Turn ppmFixer columns into strings.
    for x in ['numHexose_ppmFixer', 'numHexNAc_ppmFixer', 'numSialic_ppmFixer', 'numFucose_ppmFixer']:
        df[x] = df[x].astype('str')

    # Overwrite glycans where ppmFixSuccess = True
    df['Glycan(H,N,A,F)_ppmFixer'].mask(
        df['ppmFixSuccess'] == True,
        df['numHexose_ppmFixer'] + ' ' + df['numHexNAc_ppmFixer'] + ' ' \
            + df['numSialic_ppmFixer'] + ' ' + df['numFucose_ppmFixer'],
            inplace=True
    )

    # Note the type of correction made
    df['ppmFixerType'].mask(
        df['ppmFixSuccess'] == True,
        'NxHx -> N(x-4)H(x+5)',
        inplace=True
    )

    # Where changes are made, save new mass error.
    df['massErrorPPM_ppmfixer'].mask(
        (df['ppmFixSuccess'] == True) & (df['CorrectedGlycan(H,N,A,F)'].notna() == False),
        df['ppmFixError_temp'],
        inplace=True
    )

    # Clean up table
    df.drop(columns='ppmN6HxFlag',
            inplace=True)


#########################################
###          Script Execution         ###
#########################################
numFixes = None
newNumFixes = 0
numCycles = 1

while newNumFixes != numFixes:
    print('Beginning round ' + str(numCycles) + ' of adjustments...')
    numFixes = newNumFixes
    ppmCheck()
    ppmFix_2F_A1()
    ppmFix_4F_A2()
    ppmFix_N4HxFxAx_N2Hp7()
    ppmFix_NxHx_Nxm4Hxp5()
    newNumFixes = len(df[df['Glycan(H,N,A,F)_ppmFixer'] != df['Glycan(H,N,A,F)']])
    ppm_mean = df['massErrorPPM_ppmfixer'].mean()
    ppm_stdDev = df['massErrorPPM_ppmfixer'].std()
    print('Completed round ' + str(numCycles) + ' of adjustments.')
    numCycles += 1
    print('Performed ' + str(newNumFixes) + ' total adjustments, an addition of ' + str((newNumFixes - numFixes)) + ' adjustments since the last round.')
    print('New mean mass error is ' + str(ppm_mean)[:5] + ' PPM.')
    print('New mass error standard deviation is ' + str(ppm_stdDev)[:5] + ' PPM.')

print('Completed. Adjusted ' + str(numFixes) + ' of ' + str(len(df)) + ' GPSMs.')

#########################################
###           Table Cleanup           ###
#########################################

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
df.drop(columns=['ppmFixError_temp', 'ppmFixPrecursorMZ_temp', 'ppmCheck', 'ppmFixSuccess'],
        inplace=True)
df.drop(columns=['numHexose', 'numHexNAc', 'numSialic', 'numFucose'],
        inplace=True)
df.drop(columns=['numHexose_ppmFixer', 'numHexNAc_ppmFixer', 'numSialic_ppmFixer', 'numFucose_ppmFixer'],
        inplace=True) 
df.drop(columns=['numHexose_final', 'numHexNAc_final', 'numSialic_final', 'numFucose_final'],
        inplace=True)


# Save resulting dataframe to tsv (same as original input)
df.to_csv(output_file_path, sep = '\t')