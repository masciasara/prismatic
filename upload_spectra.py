import os
import re
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.stats import mode as scipy_mode
import gdown
import zipfile
import tkinter as tk
from tkinter import filedialog

pd.options.mode.chained_assignment = None

def download_and_extract_zip(url, output_folder):
    """
    Download a zip file from Google Drive and extract it.
    
    Args:
        url: Google Drive URL
        output_folder: Folder to extract the contents
    """
    # Download the zip file
    zip_path = os.path.join(output_folder, "data.zip")
    gdown.download(url, zip_path, quiet=False)

    # Check if the downloaded file is a valid zip file
    if not zipfile.is_zipfile(zip_path):
        raise zipfile.BadZipFile("File is not a valid zip file")

    # Extract the zip file
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(output_folder)
    
    # Remove the zip file after extraction
    os.remove(zip_path)

def get_spectra_file_paths_and_ids(folder, pointing, field):
    """
    Collect spectra file paths and organize them by source ID.
    Prioritizes x1d_optext files over regular x1d files when available.
    
    Args:
        folder: Directory containing the spectra files
        pointing: Pointing identifier used in the filename
        field: Field name (UDS, EGS, COSMOS)
        
    Returns:
        Dictionary with source IDs as keys and file paths for x1d and s2d files
    """
    spectra_dict = {}
    
    for filename in os.listdir(folder):
        # Skip hidden files
        if filename.startswith('.'):
            continue
            
        # Define the filename pattern to match based on field
        init_name = f'CAPERS_{field}_{pointing}_s'
        match = re.search(init_name + r'(\d{9})_(x1d|x1d_optext|s2d).fits', filename)
        
        if match:
            source_id = match.group(1)
            file_type = match.group(2)
            
            # Initialize entry for this source ID if it doesn't exist
            if source_id not in spectra_dict:
                spectra_dict[source_id] = {'x1d': None, 's2d': None}
                
            # Handle different file types
            if file_type == 'x1d_optext':
                # Always prefer x1d_optext over regular x1d
                spectra_dict[source_id]['x1d'] = os.path.join(folder, filename)
            elif file_type == 'x1d' and spectra_dict[source_id]['x1d'] is None:
                # Use regular x1d only if x1d_optext isn't available
                spectra_dict[source_id]['x1d'] = os.path.join(folder, filename)
            elif file_type == 's2d':
                spectra_dict[source_id]['s2d'] = os.path.join(folder, filename)
    
    return spectra_dict

def create_pointing_table(pointing, folder_path, field):
    """
    Create a table of sources for a specific pointing.
    
    Args:
        pointing: Pointing identifier
        folder_path: Path to the data folder
        field: Field name (UDS, EGS, COSMOS)
        
    Returns:
        tuple: Lists of source IDs and corresponding file paths
    """
    if field == "UDS":
        folder = os.path.join(folder_path, f'CAPERS_{field}_V0.1/{pointing}')
    else:
        folder = os.path.join(folder_path, f'Spectra/CAPERS_{field}_{pointing}')
    spectra_files = get_spectra_file_paths_and_ids(folder, pointing, field)

    source_ids = []
    file1d = []
    file2d = []

    for source_id, files in spectra_files.items():
        source_ids.append(source_id)
        file1d.append(files['x1d'])
        file2d.append(files['s2d'])

    return source_ids, file1d, file2d

def read_redshift_catalogs(pointing, field):
    """
    Read redshift catalogs for a specific pointing and field.
    
    Args:
        pointing: Pointing identifier
        field: Field name (UDS, EGS, COSMOS)
        
    Returns:
        tuple: Various redshift catalog data
    """
    # Get catalog path once
    script_dir = os.path.dirname(os.path.abspath(__file__))
    catalog_path = os.path.join(script_dir, "solutions", field)
    
    # Precompile regex patterns for better performance
    s_pattern = re.compile(r's(\d+)_')
    s_pattern_extended = re.compile(r's(\d{1,9})_')
    
    def extract_id_vectorized(series, pattern):
        """Vectorized ID extraction using pandas string methods"""
        return series.str.extract(pattern)[0]
    
    # BAGPIPES catalog
    if field == "UDS" or field == "COSMOS":
        bagpipes_file = os.path.join(catalog_path, f'capers_{field}_bagpipes_{pointing}.cat')
        BAGPIPES = pd.read_csv(bagpipes_file, sep=r'\s+')
    elif field == "EGS": 
        bagpipes_file = os.path.join(catalog_path, 'bagpipes_zfits.cat')
        BAGPIPES = pd.read_csv(bagpipes_file, sep=r'\s+')
    else:
        BAGPIPES = None
    
    # AT catalog
    at_file = os.path.join(catalog_path, f'CAPERS_{pointing}_zspec.txt')
    AT = pd.read_csv(at_file, sep=',')
    # Use vectorized string operations instead of loop
    ID_AT = AT['file'].str.extract(r's(\d+)_')[0]
    
    # MSAEXP catalog
    if field == "UDS":
        msaexp_file = os.path.join(catalog_path, f'capers_{pointing.lower()}_msaexp_zcat.csv')
    elif field == "EGS":
        msaexp_file = os.path.join(catalog_path, f'capers_{field.lower()}_msaexp_zcat.csv')
    else:
        msaexp_file = None
    
    msaexp = pd.read_csv(msaexp_file) if msaexp_file else None
    
    # LiMe catalog
    if field == "UDS":
        lime_file = os.path.join(catalog_path, f'CAPERS_{field}_V01_ASPECT_redshifts_v2.csv')
        lime = pd.read_csv(lime_file)
        ID_lime = lime['file'].astype(str).str.extract(r's(\d{1,9})_')[0]
    elif field == "EGS":
        lime_file = os.path.join(catalog_path, f'CAPERS_{field}_V0.2_lime_redshifts_v2.txt')
        lime = pd.read_csv(lime_file, sep='\s+')
        ID_lime = lime['optext'].astype(str).str.extract(r's(\d{1,9})_')[0]
    else:
        lime = None
        ID_lime = None
    
    # MARZ catalog
    if field == "UDS":
        marz_file = os.path.join(catalog_path, f'capers_{field.lower()}_{pointing.lower()}_redshift.fits')
    elif field == "EGS":
        marz_file = os.path.join(catalog_path, f'capers_{field.lower()}_{pointing.lower()}_marz_result.fits')
    else:
        marz_file = None
    
    if marz_file:
        with fits.open(marz_file) as marz:
            marz_data = marz[1].data
            MARZ = pd.DataFrame(marz_data)
        # Use vectorized string operations
        file_MARZ = MARZ['Name'].str.extract(r's(\d+)')[0]
    else:
        MARZ = None
        file_MARZ = None
    
    # CIGALE catalog
    if field == "UDS":
        cigale_file = os.path.join(catalog_path, f'CAPERS_v0.1_cigale_redshift_with_photometry.csv')
    elif field == "EGS":
        cigale_file = os.path.join(catalog_path, f'redshift-cigale-v2.csv')
    else:
        cigale_file = None
    
    CIGALE = pd.read_csv(cigale_file) if cigale_file else None
    CIGALE_ID = CIGALE['id'] if CIGALE is not None else None
    
    # Photometric catalog
    if field == "UDS":
        photo_file = os.path.join(catalog_path, f'CAPERS_{field}_master_yield_v2.1_actual.csv')
    elif field == "EGS" or field == "COSMOS":
        photo_file = os.path.join(catalog_path, f'CAPERS_{field}_master_yield_v2.csv')
    else:
        photo_file = None
    
    if photo_file:
        photometric_catalog = pd.read_csv(photo_file)
        photo_z = photometric_catalog[['ID', 'z_UNICORN', 'ra', 'dec']]
        photo_z_ID = photo_z['ID']
    else:
        photo_z = None
        photo_z_ID = None
    
    return BAGPIPES, AT, ID_AT, msaexp, lime, ID_lime, MARZ, file_MARZ, CIGALE, CIGALE_ID, photo_z, photo_z_ID

def process_field(field, google_drive_url, output_folder, pointings):
    """
    Process data for a specific field - OPTIMIZED VERSION.
    
    Args:
        field: Field name (UDS, EGS, COSMOS)
        google_drive_url: Google Drive URL for this field
        output_folder: Output folder for this field
        pointings: List of pointings to process
        
    Returns:
        DataFrame: Processed catalog for this field
    """
    print(f"\n{'='*50}")
    print(f"Processing {field} field")
    print(f"{'='*50}")
    
    # Check if we need to download data
    version_file = os.path.join(output_folder, "version.txt")
    expected_version = "v0.1"
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if os.path.exists(version_file):
        with open(version_file, 'r') as f:
            current_version = f.read().strip()
        if current_version == expected_version:
            print(f"Data for {field} is already downloaded and up-to-date.")
        else:
            print(f"Data version mismatch for {field}. Downloading the latest version.")
            download_and_extract_zip(google_drive_url, output_folder)
            with open(version_file, 'w') as f:
                f.write(expected_version)
    else:
        print(f"Downloading data for {field}...")
        download_and_extract_zip(google_drive_url, output_folder)
        with open(version_file, 'w') as f:
            f.write(expected_version)
    
    # Use lists to collect data efficiently 
    source_data = []
    catalog_data = []

    # Process each pointing for this field
    for pointing in pointings:
        print(f"Processing pointing {pointing} for {field}")
        
        # Get source IDs and file paths
        source_ids, file1d, file2d = create_pointing_table(pointing, output_folder, field)
        
        # Store source data with pointing info
        for sid, f1d, f2d in zip(source_ids, file1d, file2d):
            source_data.append({
                'galaxy_id': float(sid),
                'file1d': f1d,
                'file2d': f2d,
                'pointing': pointing
            })

        # Read redshift catalogs
        catalogs = read_redshift_catalogs(pointing, field)
        BAGPIPES, AT, ID_AT, msaexp, lime, ID_lime, MARZ, file_MARZ, CIGALE, CIGALE_ID, photo_z, photo_z_ID = catalogs
        
        # Store catalog data with pointing info for merging
        catalog_data.append({
            'pointing': pointing,
            'BAGPIPES': BAGPIPES,
            'AT': AT,
            'ID_AT': ID_AT,
            'msaexp': msaexp,
            'lime': lime,
            'ID_lime': ID_lime,
            'MARZ': MARZ,
            'file_MARZ': file_MARZ,
            'CIGALE': CIGALE,
            'CIGALE_ID': CIGALE_ID,
            'photo_z': photo_z,
            'photo_z_ID': photo_z_ID
        })

    # Convert source data to DataFrame for efficient operations
    source_df = pd.DataFrame(source_data)
    
    # Create consolidated catalog DataFrames with proper ID conversion
    def create_consolidated_catalogs(catalog_data):
        """Create consolidated DataFrames for each catalog type"""
        consolidated = {}
        
        # Initialize lists for each catalog type
        bagpipes_list = []
        at_list = []
        msaexp_list = []
        lime_list = []
        marz_list = []
        cigale_list = []
        photo_z_list = []
        
        for cat_info in catalog_data:
            pointing = cat_info['pointing']
            
            # BAGPIPES
            if cat_info['BAGPIPES'] is not None:
                bp_df = cat_info['BAGPIPES'].copy()
                bp_df['galaxy_id'] = bp_df['ID'].astype(float)
                bp_df['pointing'] = pointing
                bagpipes_list.append(bp_df[['galaxy_id', 'z_bp', 'pointing']])
            
            # AT
            if cat_info['AT'] is not None and cat_info['ID_AT'] is not None:
                at_df = cat_info['AT'].copy()
                at_df['galaxy_id'] = pd.to_numeric(cat_info['ID_AT'], errors='coerce')
                at_df['pointing'] = pointing
                at_list.append(at_df[['galaxy_id', 'z', 'pointing']].dropna(subset=['galaxy_id']))
            
            # MSAEXP
            if cat_info['msaexp'] is not None:
                msaexp_df = cat_info['msaexp'].copy()
                msaexp_df['galaxy_id'] = msaexp_df['msa_id'].astype(float)
                msaexp_df['pointing'] = pointing
                msaexp_list.append(msaexp_df[['galaxy_id', 'z_spec', 'pointing']])
            
            # LIME
            if cat_info['lime'] is not None and cat_info['ID_lime'] is not None:
                lime_df = cat_info['lime'].copy()
                lime_df['galaxy_id'] = pd.to_numeric(cat_info['ID_lime'], errors='coerce')
                lime_df['pointing'] = pointing
                lime_list.append(lime_df[['galaxy_id', 'z_manual', 'pointing']].dropna(subset=['galaxy_id']))
            
            # MARZ
            if cat_info['MARZ'] is not None and cat_info['file_MARZ'] is not None:
                marz_df = cat_info['MARZ'].copy()
                marz_df['galaxy_id'] = pd.to_numeric(cat_info['file_MARZ'], errors='coerce')
                marz_df['pointing'] = pointing
                marz_list.append(marz_df[['galaxy_id', 'AutoZ', 'pointing']].dropna(subset=['galaxy_id']))
            
            # CIGALE
            if cat_info['CIGALE'] is not None and cat_info['CIGALE_ID'] is not None:
                cigale_df = cat_info['CIGALE'].copy()
                cigale_df['galaxy_id'] = cat_info['CIGALE_ID'].astype(float)
                cigale_df['pointing'] = pointing
                cigale_list.append(cigale_df[['galaxy_id', 'best.universe.redshift', 'pointing']])
            
            # Photo-z
            if cat_info['photo_z'] is not None and cat_info['photo_z_ID'] is not None:
                photo_df = cat_info['photo_z'].copy()
                photo_df['galaxy_id'] = cat_info['photo_z_ID'].astype(float)
                photo_df['pointing'] = pointing
                photo_z_list.append(photo_df[['galaxy_id', 'z_UNICORN', 'ra', 'dec', 'pointing']])
        
        # Concatenate all DataFrames
        consolidated['bagpipes'] = pd.concat(bagpipes_list, ignore_index=True) if bagpipes_list else pd.DataFrame()
        consolidated['at'] = pd.concat(at_list, ignore_index=True) if at_list else pd.DataFrame()
        consolidated['msaexp'] = pd.concat(msaexp_list, ignore_index=True) if msaexp_list else pd.DataFrame()
        consolidated['lime'] = pd.concat(lime_list, ignore_index=True) if lime_list else pd.DataFrame()
        consolidated['marz'] = pd.concat(marz_list, ignore_index=True) if marz_list else pd.DataFrame()
        consolidated['cigale'] = pd.concat(cigale_list, ignore_index=True) if cigale_list else pd.DataFrame()
        consolidated['photo_z'] = pd.concat(photo_z_list, ignore_index=True) if photo_z_list else pd.DataFrame()
        
        return consolidated
    
    # Create consolidated catalogs
    consolidated_cats = create_consolidated_catalogs(catalog_data)
    
    # Merge all redshift data with source data using vectorized operations
    result_df = source_df.copy()
    
    # Initialize redshift columns with zeros
    redshift_columns = ['redshift_bp', 'redshift_at', 'redshift_msaexp', 
                       'redshift_lime', 'redshift_marz', 'redshift_cigale', 
                       'photo_z_values', 'ra_values', 'dec_values']
    
    for col in redshift_columns:
        result_df[col] = 0.0
    
    # Efficient merging using pandas merge operations
    def safe_merge(df, cat_df, cat_name, value_col, result_col):
        """Safely merge catalog data"""
        if not cat_df.empty:
            # For each pointing, merge the data
            merged = df.merge(
                cat_df[['galaxy_id', 'pointing', value_col]], 
                on=['galaxy_id', 'pointing'], 
                how='left'
            )
            df[result_col] = merged[value_col].fillna(0)
        return df
    
    # Merge each catalog type
    result_df = safe_merge(result_df, consolidated_cats['bagpipes'], 'bagpipes', 'z_bp', 'redshift_bp')
    result_df = safe_merge(result_df, consolidated_cats['at'], 'at', 'z', 'redshift_at')
    result_df = safe_merge(result_df, consolidated_cats['msaexp'], 'msaexp', 'z_spec', 'redshift_msaexp')
    result_df = safe_merge(result_df, consolidated_cats['lime'], 'lime', 'z_manual', 'redshift_lime')
    result_df = safe_merge(result_df, consolidated_cats['marz'], 'marz', 'AutoZ', 'redshift_marz')
    result_df = safe_merge(result_df, consolidated_cats['cigale'], 'cigale', 'best.universe.redshift', 'redshift_cigale')
    
    # Special handling for photo_z (multiple columns)
    if not consolidated_cats['photo_z'].empty:
        photo_merged = result_df.merge(
            consolidated_cats['photo_z'][['galaxy_id', 'pointing', 'z_UNICORN', 'ra', 'dec']], 
            on=['galaxy_id', 'pointing'], 
            how='left'
        )
        result_df['photo_z_values'] = photo_merged['z_UNICORN'].fillna(0)
        result_df['ra_values'] = photo_merged['ra'].fillna(0)
        result_df['dec_values'] = photo_merged['dec'].fillna(0)
    
    # Calculate redshift mode using vectorized operations
    def calculate_redshift_mode(row):
        """Calculate mode of redshift values"""
        z_values = [
            round(row['redshift_at'], 2),
            round(row['redshift_bp'], 2),
            round(row['redshift_msaexp'], 2),
            round(row['redshift_marz'], 2),
            round(row['redshift_lime'], 2)
        ]
        # Remove zeros for mode calculation
        z_values = [z for z in z_values if z != 0]
        if z_values:
            return scipy_mode(z_values)[0]
        return 0
    
    result_df['redshifts_mode'] = result_df.apply(calculate_redshift_mode, axis=1)
    
    # Create final catalog DataFrame
    catalog = pd.DataFrame({
        "Galaxy": result_df['galaxy_id'].astype(int),
        "1d_Spectrum": result_df['file1d'],
        "2d_Spectrum": result_df['file2d'],
        "RA": result_df['ra_values'],
        "DEC": result_df['dec_values'],
        "photo_Redshift": result_df['photo_z_values'],
        "Redshift": np.zeros(len(result_df)),  # Initialize with zeros as in original
        "Field": [field] * len(result_df),
        "Pointing": result_df['pointing'],
        "Flag": [""] * len(result_df),
        "AT_Redshift": result_df['redshift_at'],
        "MSAEXP_Redshift": result_df['redshift_msaexp'],
        "LiMe_Redshift": result_df['redshift_lime'],
        "MARZ_Redshift": result_df['redshift_marz'],
        "Cigale_Redshift": result_df['redshift_cigale'],
        "BAGPIPES_Redshift": result_df['redshift_bp'],
        "Redshift_Mode": result_df['redshifts_mode'],
        "[OIII]+Hβ": [""] * len(result_df),
        "Hα": [""] * len(result_df),
        "Lyα": [""] * len(result_df),
        "Paα": [""] * len(result_df),
        "Lyman Break": [""] * len(result_df),
        "Balmer Break": [""] * len(result_df),
        "Comments": [""] * len(result_df)
    })

    # Clean up and filter the catalog
    catalog["Redshift"] = catalog["Redshift"].apply(lambda x: 0 if pd.isna(x) or x <= 0 else x)
    catalog_filtered = catalog[catalog["1d_Spectrum"].notnull() & (catalog["1d_Spectrum"] != "")]
    catalog_sorted = catalog_filtered.sort_values(by="Galaxy").reset_index(drop=True)
    
    return catalog_sorted


# Alternative version with even more aggressive optimization for very large datasets
def process_field_ultra_optimized(field, google_drive_url, output_folder, pointings):
    """
    Ultra-optimized version using pure vectorized operations and minimal memory usage.
    Use this for very large datasets (>100k galaxies).
    """
    print(f"\n{'='*50}")
    print(f"Processing {field} field (Ultra-Optimized)")
    print(f"{'='*50}")
    
    # Download logic (same as original)
    version_file = os.path.join(output_folder, "version.txt")
    expected_version = "v0.1"
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if os.path.exists(version_file):
        with open(version_file, 'r') as f:
            current_version = f.read().strip()
        if current_version == expected_version:
            print(f"Data for {field} is already downloaded and up-to-date.")
        else:
            print(f"Data version mismatch for {field}. Downloading the latest version.")
            download_and_extract_zip(google_drive_url, output_folder)
            with open(version_file, 'w') as f:
                f.write(expected_version)
    else:
        print(f"Downloading data for {field}...")
        download_and_extract_zip(google_drive_url, output_folder)
        with open(version_file, 'w') as f:
            f.write(expected_version)
    
    # Collect all data in single pass
    all_data = []
    
    for pointing in pointings:
        print(f"Processing pointing {pointing} for {field}")
        
        # Get source data
        source_ids, file1d, file2d = create_pointing_table(pointing, output_folder, field)
        
        # Read catalogs
        catalogs = read_redshift_catalogs(pointing, field)
        BAGPIPES, AT, ID_AT, msaexp, lime, ID_lime, MARZ, file_MARZ, CIGALE, CIGALE_ID, photo_z, photo_z_ID = catalogs
        
        # Create lookup dictionaries for O(1) access
        lookups = {}
        
        # Build lookup dictionaries
        if BAGPIPES is not None:
            lookups['bagpipes'] = dict(zip(BAGPIPES['ID'].astype(float), BAGPIPES['z_bp']))
        
        if AT is not None and ID_AT is not None:
            valid_ids = pd.to_numeric(ID_AT, errors='coerce').dropna()
            lookups['at'] = dict(zip(valid_ids, AT['z'].iloc[:len(valid_ids)]))
        
        if msaexp is not None:
            lookups['msaexp'] = dict(zip(msaexp['msa_id'].astype(float), msaexp['z_spec']))
        
        if lime is not None and ID_lime is not None:
            valid_ids = pd.to_numeric(ID_lime, errors='coerce').dropna()
            lookups['lime'] = dict(zip(valid_ids, lime['z_manual'].iloc[:len(valid_ids)]))
        
        if MARZ is not None and file_MARZ is not None:
            valid_ids = pd.to_numeric(file_MARZ, errors='coerce').dropna()
            lookups['marz'] = dict(zip(valid_ids, MARZ['AutoZ'].iloc[:len(valid_ids)]))
        
        if CIGALE is not None and CIGALE_ID is not None:
            lookups['cigale'] = dict(zip(CIGALE_ID.astype(float), CIGALE['best.universe.redshift']))
        
        if photo_z is not None and photo_z_ID is not None:
            photo_ids = photo_z_ID.astype(float)
            lookups['photo_z'] = dict(zip(photo_ids, zip(photo_z['z_UNICORN'], photo_z['ra'], photo_z['dec'])))
        
        # Process each source using lookups
        for sid, f1d, f2d in zip(source_ids, file1d, file2d):
            galaxy_id = float(sid)
            
            # Get redshift values using dictionary lookups (O(1))
            redshift_bp = lookups.get('bagpipes', {}).get(galaxy_id, 0)
            redshift_at = lookups.get('at', {}).get(galaxy_id, 0)
            redshift_msaexp = lookups.get('msaexp', {}).get(galaxy_id, 0)
            redshift_lime = lookups.get('lime', {}).get(galaxy_id, 0)
            redshift_marz = lookups.get('marz', {}).get(galaxy_id, 0)
            redshift_cigale = lookups.get('cigale', {}).get(galaxy_id, 0)
            
            # Photo-z data
            photo_data = lookups.get('photo_z', {}).get(galaxy_id, (0, 0, 0))
            photo_z_val, ra_val, dec_val = photo_data
            
            # Calculate mode
            z_values = [round(z, 2) for z in [redshift_at, redshift_bp, redshift_msaexp, redshift_marz, redshift_lime] if z != 0]
            redshift_mode = scipy_mode(z_values)[0] if z_values else 0
            
            all_data.append({
                'Galaxy': int(galaxy_id),
                '1d_Spectrum': f1d,
                '2d_Spectrum': f2d,
                'RA': ra_val,
                'DEC': dec_val,
                'photo_Redshift': photo_z_val,
                'Redshift': 0,
                'Field': field,
                'Pointing': pointing,
                'Flag': "",
                'AT_Redshift': redshift_at,
                'MSAEXP_Redshift': redshift_msaexp,
                'LiMe_Redshift': redshift_lime,
                'MARZ_Redshift': redshift_marz,
                'Cigale_Redshift': redshift_cigale,
                'BAGPIPES_Redshift': redshift_bp,
                'Redshift_Mode': redshift_mode,
                '[OIII]+Hβ': "",
                'Hα': "",
                'Lyα': "",
                'Paα': "",
                'Lyman Break': "",
                'Balmer Break': "",
                'Comments': ""
            })
    
    # Create DataFrame from collected data
    catalog = pd.DataFrame(all_data)
    
    # Filter and sort
    catalog_filtered = catalog[catalog["1d_Spectrum"].notnull() & (catalog["1d_Spectrum"] != "")]
    catalog_sorted = catalog_filtered.sort_values(by="Galaxy").reset_index(drop=True)
    
    return catalog_sorted

def main():
    """
    Main function to process data from multiple fields (UDS, EGS, COSMOS).
    """
    # Define data for each field
    fields_data = {
        "UDS": {
            "google_drive_url": "https://drive.google.com/uc?id=1nLwEX3edQI9ktHof9eAsH0z7lvUYFz7i",
            "output_folder": "data_UDS",
            "pointings": ['P1', 'P2', 'P3', 'P5']
        },
        "EGS": {
            "google_drive_url": "https://drive.google.com/uc?id=12ilzJDKhUdu00nbsU_4dxwdWIxuu61xA",
            "output_folder": "data_EGS",
            "pointings": ['P1', 'P2', 'P3', 'P4', 'P5', 'P7']  # Update with actual EGS pointings
        },
        "COSMOS": {
            "google_drive_url": "https://drive.google.com/uc?id=1lZri9sfRakE9lQwsWhaXbVTqIZNl0IUw",
            "output_folder": "data_COSMOS",
            "pointings": ['P1', 'P2', 'P3', 'P4', 'P5', 'P7']  # Update with actual COSMOS pointings
        }
    }
    
    # Ask user which field to process
    print("Available fields: UDS, EGS, COSMOS")
    selected_field = input("Enter field to process (or 'all' for all fields): ").strip().upper()
    
    all_catalogs = {}
    
    if selected_field == "ALL":
        fields_to_process = list(fields_data.keys())
    else:
        fields_to_process = [selected_field]
    
    # Process each selected field
    for field in fields_to_process:
        if field not in fields_data:
            print(f"Error: {field} is not a valid field. Skipping.")
            continue
        
        field_info = fields_data[field]
        catalog = process_field(
            field, 
            field_info["google_drive_url"], 
            field_info["output_folder"],
            field_info["pointings"]
        )
        
        all_catalogs[field] = catalog
    
    # Combine catalogs if multiple fields were processed
    if len(all_catalogs) > 1:
        combined_catalog = pd.concat(list(all_catalogs.values())).reset_index(drop=True)
        print(f"Combined catalog created with {len(combined_catalog)} total entries.")
        catalog_sorted = combined_catalog
        field_name = "combined"
    else:
        field_name = list(all_catalogs.keys())[0]
        catalog_sorted = all_catalogs[field_name]
    
    # Visualization options
    visualize = input("Do you want to visualize sources within a specific ID range, a specific pointing, a specific group or a specific field? (range/pointing/field/upload/all): ").strip().lower()
    
    if visualize == 'range':
        min_id = int(input("Enter the minimum ID: "))
        max_id = int(input("Enter the maximum ID: "))
        catalog_sorted = catalog_sorted[(catalog_sorted['Galaxy'].astype(int) >= min_id) & (catalog_sorted['Galaxy'].astype(int) <= max_id)]
    elif visualize == 'pointing':
        chosen_pointing = input("Enter the pointing (e.g., P1, P2, P3, P5): ").strip().upper()
        catalog_sorted = catalog_sorted[catalog_sorted['Pointing'] == chosen_pointing]
    elif visualize == 'upload':
        root = tk.Tk()
        root.withdraw()
        file_path = filedialog.askopenfilename(title="Select CSV file", filetypes=[("CSV files", "*.csv")])
        if file_path:
            uploaded_catalog = pd.read_csv(file_path)
            galaxy_ids = uploaded_catalog['Galaxy'].astype(int).tolist()
            catalog_sorted = catalog_sorted[catalog_sorted['Galaxy'].astype(int).isin(galaxy_ids)]
    elif visualize == 'field' and len(all_catalogs) > 1:
        chosen_field = input("Enter the field name (UDS, EGS, COSMOS): ").strip().upper()
        catalog_sorted = catalog_sorted[catalog_sorted['Field'] == chosen_field]
    
    # Add index column
    catalog_sorted.insert(0, 'Index', range(len(catalog_sorted)))
    
    # Save the catalog
    default_catalog_name = f'initial_catalog_{field_name}.csv'
    catalog_sorted.to_csv(default_catalog_name, index=False)
    print(f"Catalog saved as {default_catalog_name}")
    
    # Let user choose a custom name
    new_catalog_name = input("Choose a name for the output catalog (remember the extension!): (e.g., final_catalog_YourName.csv) ")
    catalog_sorted.to_csv(new_catalog_name, index=False)
    print(f"Catalog saved as {new_catalog_name}")
    
    return catalog_sorted, new_catalog_name

if __name__ == "__main__":
    catalog_sorted, new_catalog_name = main()
