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
    
    # Check if folder exists and is accessible
    if not os.path.exists(folder):
        print(f"Warning: Folder does not exist: {folder}")
        return spectra_dict
    
    try:
        filenames = os.listdir(folder)
    except PermissionError:
        print(f"Warning: Permission denied for folder: {folder}")
        return spectra_dict
    
    for filename in filenames:
        # Skip hidden files
        if filename.startswith('.'):
            continue
            
        # Define the filename pattern to match based on field
        init_name = f'CAPERS-{field}_{pointing}_s'
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
    
    if not spectra_dict:
        print(f"Warning: No spectra files found in {folder}")

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
    
    folder = os.path.join(folder_path, f'CAPERS-{field}_DR0.4/Spectra/CAPERS-{field}_{pointing}/')
    
    # Check if folder exists
    if not os.path.exists(folder):
        print(f"Warning: Spectra folder not found: {folder}")
        return [], [], []
    
    spectra_files = get_spectra_file_paths_and_ids(folder, pointing, field)

    source_ids = []
    file1d = []
    file2d = []

    for source_id, files in spectra_files.items():
        source_ids.append(source_id)
        file1d.append(files['x1d'])
        file2d.append(files['s2d'])

    return source_ids, file1d, file2d

def load_catalog_safe(filepath, file_type='csv', **kwargs):
    """
    Safely load a catalog file, return None if file doesn't exist.
    
    Args:
        filepath: Path to the file
        file_type: Type of file ('csv', 'fits', 'txt')
        **kwargs: Additional arguments for reading functions
        
    Returns:
        DataFrame or None if file doesn't exist
    """
    try:
        if not os.path.exists(filepath):
            print(f"Warning: File not found - {filepath}")
            return None
            
        if file_type == 'fits':
            with fits.open(filepath) as hdul:
                data = hdul[1].data
                return pd.DataFrame(data)
        elif file_type == 'txt':
            return pd.read_csv(filepath, **kwargs)
        else:  # csv
            return pd.read_csv(filepath, **kwargs)
            
    except (FileNotFoundError, IOError) as e:
        print(f"Warning: File not found - {filepath}")
        return None
    except Exception as e:
        print(f"Warning: Error loading {filepath}: {str(e)}")
        return None

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
    
    # BAGPIPES catalog
    if field == "UDS" or field == "COSMOS":
        bagpipes_file = os.path.join(catalog_path, f'capers_{field}_bagpipes_{pointing}.cat')
        BAGPIPES = load_catalog_safe(bagpipes_file, file_type='txt', sep=r'\s+')
    elif field == "EGS": 
        bagpipes_file = os.path.join(catalog_path, f'capers_{field}_bagpipes_{pointing}.cat')
        BAGPIPES = load_catalog_safe(bagpipes_file, file_type='txt', sep=r'\s+')
    else:
        BAGPIPES = None
    
    # AT catalog - using comprehensive catalog
    at_file = os.path.join(os.path.dirname(catalog_path), 'AT_all_zfits.csv')
    AT_all = load_catalog_safe(at_file, file_type='csv')
    
    # Filter for current field only
    if AT_all is not None:
        AT_all['Pointing'] = 'P' + AT_all['Pointing'].astype(str)
        # Filter by Field only (capitalized column name)
        AT = AT_all[AT_all['Field'] == field].copy()
        ID_AT = AT['ID'] if not AT.empty else None
    
    # MSAEXP catalog
    if field == "UDS":
        msaexp_file = os.path.join(catalog_path, f'capers_{pointing.lower()}_msaexp_zcat.csv')
    elif field == "EGS":
        msaexp_file = os.path.join(catalog_path, f'capers_{field.lower()}_msaexp_zcat.csv')
    else:
        msaexp_file = None
    
    msaexp = load_catalog_safe(msaexp_file, file_type='csv') if msaexp_file else None
    
    # LiMe catalog
    if field == "UDS":
        lime_file = os.path.join(catalog_path, f'CAPERS_{field}_V0.4_automatic_redshifts_inspected_v0.5.csv')
        lime = load_catalog_safe(lime_file, file_type='csv')
        ID_lime = lime['id'].astype(str).str.extract(r's(\d{1,9})_')[0] if lime is not None else None
    elif field == "EGS":
        lime_file = os.path.join(catalog_path, f'CAPERS_{field}_V0.2_lime_redshifts_v2.txt')
        lime = load_catalog_safe(lime_file, file_type='txt', sep=r'\s+')
        ID_lime = lime['optext'].astype(str).str.extract(r's(\d{1,9})_')[0] if lime is not None else None
    elif field == "COSMOS":
        lime_file = os.path.join(catalog_path, f'CAPERS_{field}_V0.2.1_redshifts_aspectv0.3_and_manual_inspection.csv')
        lime = load_catalog_safe(lime_file, file_type='csv')
        ID_lime = lime['file_name'].astype(str).str.extract(r's(\d{1,9})_')[0] if lime is not None else None
    else:
        lime = None
        ID_lime = None
    
    # MARZ catalog
    if field == "UDS":
        marz_file = os.path.join(catalog_path, f'capers_{field.lower()}_{pointing.lower()}_redshift.fits')
        MARZ = load_catalog_safe(marz_file, file_type='fits')
        file_MARZ = MARZ['Name'].str.extract(r's(\d+)')[0] if MARZ is not None else None
    elif field == "EGS":
        if pointing == "P6":
            marz_file = os.path.join(catalog_path, f'capers_{field.lower()}_{pointing.lower()}_vetted.csv')
            MARZ = load_catalog_safe(marz_file, file_type='csv')
            file_MARZ = MARZ['id'] if MARZ is not None else None
        else:
            marz_file = os.path.join(catalog_path, f'capers_{field.lower()}_{pointing.lower()}_marz_result.fits')
            MARZ = load_catalog_safe(marz_file, file_type='fits')
            file_MARZ = MARZ['Name'].str.extract(r's(\d+)')[0] if MARZ is not None else None
    elif field == "COSMOS":
        marz_file = os.path.join(catalog_path, f'capers_{field.lower()}_{pointing.lower()}_marz_result.fits')
        MARZ = load_catalog_safe(marz_file, file_type='fits')
        file_MARZ = MARZ['Name'].str.extract(r's(\d+)')[0] if MARZ is not None else None
    else:
        MARZ = None
        file_MARZ = None
		
    # CIGALE catalog
    if field == "UDS":
        cigale_file = os.path.join(catalog_path, f'redshifts_CIGALE.csv')
    elif field == "EGS":
        cigale_file = os.path.join(catalog_path, f'redshift-cigale-v2.csv')
    elif field == "COSMOS":
        cigale_file = os.path.join(catalog_path, f'CIGALE-redshift-COSMOS-V0.2.1.csv')
    else:
        cigale_file = None
    
    CIGALE = load_catalog_safe(cigale_file, file_type='csv') if cigale_file else None
    CIGALE_ID = CIGALE['id'] if CIGALE is not None else None

    # Photometric catalog
    photo_file = os.path.join(catalog_path, f'CAPERS-{field}_master_yield_v6.csv')
    
    if photo_file:
        photometric_catalog = load_catalog_safe(photo_file, file_type='csv')
        if photometric_catalog is not None:
            photo_z = photometric_catalog[['ID', 'z_UNICORN', 'ra', 'dec']]
            photo_z_ID = photo_z['ID']
        else:
            photo_z = None
            photo_z_ID = None
    else:
        photo_z = None
        photo_z_ID = None
    
    return BAGPIPES, AT, ID_AT, msaexp, lime, ID_lime, MARZ, file_MARZ, CIGALE, CIGALE_ID, photo_z, photo_z_ID

def process_field(field, google_drive_url, output_folder, pointings):
    """
    Process data for a specific field - IMPROVED VERSION with better data preservation.
    
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
    
    # Use photometric catalog as the base
    print(f"\nLoading base photometric catalog for {field}...")
    
    # Load the photometric catalog for this field
    photo_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "solutions", field, f'CAPERS-{field}_master_yield_v6.csv')
    
    if photo_file and os.path.exists(photo_file):
        base_catalog = pd.read_csv(photo_file)
        print(f"Loaded {len(base_catalog)} sources from photometric catalog")
        
        # Create base DataFrame with galaxy info
        result_df = pd.DataFrame({
            'galaxy_id': base_catalog['ID'].astype(float),
            'ra_values': base_catalog['ra'],
            'dec_values': base_catalog['dec'],
            'photo_z_values': base_catalog['z_UNICORN']
        })
    else:
        print(f"Warning: Photometric catalog not found at {photo_file}")
        print(f"Will create catalog from spectra and redshift data only")
        result_df = pd.DataFrame()
    
    # Collect all galaxy IDs from spectra and redshift catalogs
    all_galaxy_ids = set()
    
    # Add file paths by matching with spectra files
    print(f"\nMatching spectra files...")
    source_data = []
    
    # Process each pointing for this field to get spectra file paths
    for pointing in pointings:
        print(f"Processing pointing {pointing} for {field}")
        
        # Get source IDs and file paths
        source_ids, file1d, file2d = create_pointing_table(pointing, output_folder, field)
        
        print(f"  Found {len(source_ids)} spectra files in pointing {pointing}")
        
        # Store source data with pointing info
        for sid, f1d, f2d in zip(source_ids, file1d, file2d):
            galaxy_id = float(sid)
            all_galaxy_ids.add(galaxy_id)
            source_data.append({
                'galaxy_id': galaxy_id,
                'file1d': f1d,
                'file2d': f2d,
                'pointing': pointing
            })
    
    # Convert spectra data to DataFrame
    if source_data:
        spectra_df = pd.DataFrame(source_data)
        print(f"Total spectra files found: {len(spectra_df)}")
    else:
        print("Warning: No spectra files found for any pointing")
        spectra_df = pd.DataFrame(columns=['galaxy_id', 'file1d', 'file2d', 'pointing'])
    
    # Prepare catalog data for redshift merging
    catalog_data = []
    
    # Read redshift catalogs for each pointing
    for pointing in pointings:
        print(f"Reading redshift catalogs for {pointing}")
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
    
    # Create consolidated catalog DataFrames with proper ID conversion
    def create_consolidated_catalogs(catalog_data, all_galaxy_ids):
        """Create consolidated DataFrames for each catalog type and collect all galaxy IDs"""
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
                # Add galaxy IDs to the set
                all_galaxy_ids.update(bp_df['galaxy_id'].tolist())
                bagpipes_list.append(bp_df[['galaxy_id', 'z_bp', 'pointing']])
            
            # AT - updated to use direct ID column
            # Note: AT is filtered by field only, not by pointing
            if cat_info['AT'] is not None and not cat_info['AT'].empty:
                at_df = cat_info['AT'].copy()
                at_df['galaxy_id'] = pd.to_numeric(at_df['ID'], errors='coerce')
                # Pointing is already in the AT dataframe from filtering (capitalized)
                # Convert to string to match other catalogs
                at_df['pointing'] = at_df['Pointing']
                at_clean = at_df[['galaxy_id', 'zfit', 'pointing']].dropna(subset=['galaxy_id'])
                # Add galaxy IDs to the set
                all_galaxy_ids.update(at_clean['galaxy_id'].tolist())
                at_list.append(at_clean)
            
            # MSAEXP
            if cat_info['msaexp'] is not None:
                msaexp_df = cat_info['msaexp'].copy()
                msaexp_df['galaxy_id'] = msaexp_df['msa_id'].astype(float)
                msaexp_df['pointing'] = pointing
                # Add galaxy IDs to the set
                all_galaxy_ids.update(msaexp_df['galaxy_id'].tolist())
                msaexp_list.append(msaexp_df[['galaxy_id', 'z_spec', 'pointing']])
            
            # LIME
            if cat_info['lime'] is not None and cat_info['ID_lime'] is not None:
                lime_df = cat_info['lime'].copy()
                lime_df['galaxy_id'] = pd.to_numeric(cat_info['ID_lime'], errors='coerce')
                lime_df['pointing'] = pointing
                lime_clean = lime_df[['galaxy_id', 'z_manual', 'pointing']].dropna(subset=['galaxy_id'])
                # Add galaxy IDs to the set
                all_galaxy_ids.update(lime_clean['galaxy_id'].tolist())
                lime_list.append(lime_clean)
            
            # MARZ processing
            if cat_info['MARZ'] is not None and cat_info['file_MARZ'] is not None:
                marz_df = cat_info['MARZ'].copy()
                if field == "UDS":
                    marz_df['galaxy_id'] = pd.to_numeric(cat_info['file_MARZ'], errors='coerce')
                    marz_df['pointing'] = pointing
                    marz_df['redshift'] = marz_df['AutoZ']
                    # UDS uses 'redshift' column directly from the FITS file
                    if 'AutoZ' not in marz_df.columns:
                        print(f"Warning: 'AutoZ' column not found in MARZ data for {pointing}")
                        continue
                if field == "EGS" and pointing != "P6":
                    marz_df['galaxy_id'] = pd.to_numeric(cat_info['file_MARZ'], errors='coerce')
                    marz_df['pointing'] = pointing
                    marz_df['redshift'] = marz_df['AutoZ'] 
                if field == "EGS" and pointing == "P6":
                    marz_df['galaxy_id'] = cat_info['file_MARZ']
                    marz_df['pointing'] = "P6"
                    marz_df['redshift'] = marz_df['redshift'] 
                if field == "COSMOS":
                    marz_df['galaxy_id'] = pd.to_numeric(cat_info['file_MARZ'], errors='coerce')
                    marz_df['pointing'] = pointing
                    marz_df['redshift'] = marz_df['AutoZ'] 

                marz_clean = marz_df[['galaxy_id', 'redshift', 'pointing']].dropna(subset=['galaxy_id'])
                # Add galaxy IDs to the set
                all_galaxy_ids.update(marz_clean['galaxy_id'].tolist())
                marz_list.append(marz_clean)
            
            # CIGALE
            if cat_info['CIGALE'] is not None and cat_info['CIGALE_ID'] is not None:
                cigale_df = cat_info['CIGALE'].copy()
                cigale_df['galaxy_id'] = cat_info['CIGALE_ID'].astype(float)
                cigale_df['pointing'] = pointing
                # Add galaxy IDs to the set
                all_galaxy_ids.update(cigale_df['galaxy_id'].tolist())
                cigale_list.append(cigale_df[['galaxy_id', 'best.universe.redshift']])
            
            # Photo-z
            if cat_info['photo_z'] is not None and cat_info['photo_z_ID'] is not None:
                photo_df = cat_info['photo_z'].copy()
                photo_df['galaxy_id'] = cat_info['photo_z_ID'].astype(float)
                photo_df['pointing'] = pointing
                photo_z_list.append(photo_df[['galaxy_id', 'z_UNICORN', 'ra', 'dec', 'pointing']])
        
        # Concatenate all DataFrames with proper empty DataFrame fallbacks
        consolidated['bagpipes'] = (pd.concat(bagpipes_list, ignore_index=True) if bagpipes_list 
                                    else pd.DataFrame(columns=['galaxy_id', 'z_bp', 'pointing']))
        # AT: drop duplicates since we filter by field only (same data added multiple times)
        consolidated['at'] = (pd.concat(at_list, ignore_index=True).drop_duplicates(subset=['galaxy_id', 'pointing']) if at_list 
                              else pd.DataFrame(columns=['galaxy_id', 'zfit', 'pointing']))
        consolidated['msaexp'] = (pd.concat(msaexp_list, ignore_index=True) if msaexp_list 
                                  else pd.DataFrame(columns=['galaxy_id', 'z_spec', 'pointing']))
        consolidated['lime'] = (pd.concat(lime_list, ignore_index=True) if lime_list 
                                else pd.DataFrame(columns=['galaxy_id', 'z_manual', 'pointing']))
        consolidated['marz'] = (pd.concat(marz_list, ignore_index=True) if marz_list 
                                else pd.DataFrame(columns=['galaxy_id', 'redshift', 'pointing']))
        consolidated['cigale'] = (pd.concat(cigale_list, ignore_index=True) if cigale_list 
                                  else pd.DataFrame(columns=['galaxy_id', 'best.universe.redshift', 'pointing']))
        consolidated['photo_z'] = (pd.concat(photo_z_list, ignore_index=True) if photo_z_list 
                                   else pd.DataFrame(columns=['galaxy_id', 'z_UNICORN', 'ra', 'dec', 'pointing']))
        
        return consolidated
    
    # Create consolidated catalogs and collect all galaxy IDs
    consolidated_cats = create_consolidated_catalogs(catalog_data, all_galaxy_ids)
    
    print(f"\nTotal unique galaxy IDs found across all catalogs: {len(all_galaxy_ids)}")
    
    # Create comprehensive result_df from all sources
    if not result_df.empty:
        # We have a photometric catalog base
        print(f"Photometric catalog IDs: {len(result_df)}")
        # Add any IDs from redshift catalogs that aren't in photometric catalog
        photo_ids = set(result_df['galaxy_id'].tolist())
        missing_ids = all_galaxy_ids - photo_ids
        if missing_ids:
            print(f"Found {len(missing_ids)} additional IDs in redshift catalogs not in photometric catalog")
            missing_df = pd.DataFrame({
                'galaxy_id': list(missing_ids),
                'ra_values': 0.0,
                'dec_values': 0.0,
                'photo_z_values': 0.0
            })
            result_df = pd.concat([result_df, missing_df], ignore_index=True)
    else:
        # No photometric catalog, create from all collected IDs
        print(f"Creating catalog from {len(all_galaxy_ids)} IDs found in redshift catalogs and spectra")
        result_df = pd.DataFrame({
            'galaxy_id': list(all_galaxy_ids),
            'ra_values': 0.0,
            'dec_values': 0.0,
            'photo_z_values': 0.0
        })
    
    # Initialize file path columns BEFORE merging
    result_df['file1d'] = ""
    result_df['file2d'] = ""
    result_df['pointing'] = ""
    
    # Merge spectra file paths - IMPROVED: handle duplicates by keeping first match
    if not spectra_df.empty:
        #print(f"\nMerging {len(spectra_df)} spectra file paths...")
        # Remove duplicates in spectra_df, keeping first occurrence
        spectra_df_unique = spectra_df.drop_duplicates(subset=['galaxy_id'], keep='first')
        #print(f"After removing duplicates: {len(spectra_df_unique)} unique spectra entries")
        
        # Merge using a left join to preserve all rows in result_df
        result_df = result_df.merge(
            spectra_df_unique[['galaxy_id', 'file1d', 'file2d', 'pointing']], 
            on='galaxy_id', 
            how='left',
            suffixes=('_old', '')
        )
        
        # Drop old columns if they exist
        result_df = result_df.drop(columns=[col for col in result_df.columns if col.endswith('_old')])
        
        # Fill NaN values with empty strings
        result_df['file1d'] = result_df['file1d'].fillna("")
        result_df['file2d'] = result_df['file2d'].fillna("")
        result_df['pointing'] = result_df['pointing'].fillna("")
        # Verify merge
        spectra_count = (result_df['file1d'] != "").sum()
        #print(f"Successfully merged file paths for {spectra_count} sources")
    
    #print(f"Final merged catalog has {len(result_df)} total sources")
    
    # Initialize redshift columns with NaN (will be filled with 0 later for display)
    redshift_columns = ['redshift_bp', 'redshift_at', 'redshift_msaexp', 
                       'redshift_lime', 'redshift_marz', 'redshift_cigale']
    
    for col in redshift_columns:
        result_df[col] = np.nan
    
    # Efficient merging using pandas merge operations - IMPROVED
    def safe_merge(df, cat_df, cat_name, value_col, result_col):
        """Safely merge catalog data, handling duplicates and preserving data"""
        if not cat_df.empty:
            # Check if required columns exist
            required_cols = ['galaxy_id', value_col]
            
            # Determine merge strategy based on pointing availability
            has_pointing_in_df = 'pointing' in df.columns and (df['pointing'] != "").any()
            has_pointing_in_cat = 'pointing' in cat_df.columns
            
            if has_pointing_in_df and has_pointing_in_cat:
                required_cols.append('pointing')
                merge_on = ['galaxy_id', 'pointing']
            else:
                merge_on = ['galaxy_id']
            
            missing_cols = [col for col in required_cols if col not in cat_df.columns]
            if missing_cols:
                print(f"Warning: {cat_name} catalog missing columns {missing_cols}. Skipping merge.")
                return df
            
            # Remove duplicates from catalog, keeping first occurrence
            cat_df_unique = cat_df[required_cols].drop_duplicates(subset=merge_on, keep='first')
            
            #print(f"  Merging {cat_name}: {len(cat_df_unique)} unique entries")
            
            # Merge the data
            merged = df.merge(
                cat_df_unique, 
                on=merge_on, 
                how='left',
                suffixes=('', '_new')
            )
            
            # Update the result column, preserving existing non-NaN values
            if value_col + '_new' in merged.columns:
                # Use new values where old is NaN
                df[result_col] = merged[value_col + '_new'].combine_first(merged.get(result_col, pd.Series(dtype=float)))
            elif value_col in merged.columns:
                df[result_col] = merged[value_col].combine_first(df.get(result_col, pd.Series(dtype=float)))
            
            matches = df[result_col].notna().sum()
            #print(f"  {cat_name} matched {matches} sources")
        
        return df
    
    # Merge each catalog type - updated AT to use 'zfit' column
    print("\nMerging redshift catalogs...")
    result_df = safe_merge(result_df, consolidated_cats['bagpipes'], 'BAGPIPES', 'z_bp', 'redshift_bp')
    result_df = safe_merge(result_df, consolidated_cats['at'], 'AT', 'zfit', 'redshift_at')
    result_df = safe_merge(result_df, consolidated_cats['msaexp'], 'MSAEXP', 'z_spec', 'redshift_msaexp')
    result_df = safe_merge(result_df, consolidated_cats['lime'], 'LiMe', 'z_manual', 'redshift_lime')
    result_df = safe_merge(result_df, consolidated_cats['marz'], 'MARZ', 'redshift', 'redshift_marz')
    result_df = safe_merge(result_df, consolidated_cats['cigale'], 'CIGALE', 'best.universe.redshift', 'redshift_cigale')
    
    # Calculate redshift mode using vectorized operations - IMPROVED
    def calculate_redshift_mode(row):
        """Calculate mode of redshift values, ignoring zeros and NaNs"""
        z_values = [
            row['redshift_at'],
            row['redshift_bp'],
            row['redshift_msaexp'],
            row['redshift_marz'],
            row['redshift_lime']
        ]
        
        # Remove NaN and zeros, then round to 2 decimal places
        z_values = [round(z, 2) for z in z_values if pd.notna(z) and z > 0]
        
        if z_values:
            mode_result = scipy_mode(z_values, keepdims=True)
            return mode_result[0][0]
        return np.nan
    
    #print("\nCalculating redshift mode...")
    result_df['redshifts_mode'] = result_df.apply(calculate_redshift_mode, axis=1)
    
    # Calculate best redshift estimate - prefer mode, then individual catalog values
    def get_best_redshift(row):
        """Get the best available redshift estimate"""
        # First try mode
        if pd.notna(row['redshifts_mode']) and row['redshifts_mode'] > 0:
            return row['redshifts_mode']
        
        # Try each catalog in order of preference
        for z_col in ['redshift_lime', 'redshift_at', 'redshift_msaexp', 'redshift_marz', 'redshift_bp']:
            if pd.notna(row[z_col]) and row[z_col] > 0:
                return row[z_col]
        
        # Fall back to photometric redshift
        if pd.notna(row['photo_z_values']) and row['photo_z_values'] > 0:
            return row['photo_z_values']
        
        return 0.0
    
    result_df['best_redshift'] = result_df.apply(get_best_redshift, axis=1)
    
    # Replace NaN with 0 for display purposes
    for col in redshift_columns + ['redshifts_mode']:
        result_df[col] = result_df[col].fillna(0)
    
    # Create final catalog DataFrame - IMPROVED with all data preserved
    catalog = pd.DataFrame({
        "Galaxy": result_df['galaxy_id'].astype(int),
        "1d_Spectrum": result_df['file1d'],
        "2d_Spectrum": result_df['file2d'],
        "RA": result_df['ra_values'],
        "DEC": result_df['dec_values'],
        "photo_Redshift": result_df['photo_z_values'],
        "Redshift": result_df['best_redshift'],  # Use calculated best redshift
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

    # Clean up the catalog - only set invalid values to 0
    catalog["Redshift"] = catalog["Redshift"].apply(lambda x: 0 if pd.isna(x) or x < 0 else x)
    
    # Sort by Galaxy ID
    catalog_sorted = catalog.sort_values(by="Galaxy").reset_index(drop=True)
    
    # Print summary statistics
    print(f"\n{'='*50}")
    print(f"CATALOG SUMMARY FOR {field}")
    #print(f"{'='*50}")
    #print(f"Total sources: {len(catalog_sorted)}")
    spectra_count = (catalog_sorted["1d_Spectrum"] != "").sum()
    #print(f"Sources with 1D spectra: {spectra_count}")
    spectra_2d_count = (catalog_sorted["2d_Spectrum"] != "").sum()
    #print(f"Sources with 2D spectra: {spectra_2d_count}")
    #print(f"Sources without spectra: {len(catalog_sorted) - spectra_count}")
    
    # Redshift statistics
    redshift_with_data = (catalog_sorted["Redshift"] > 0).sum()
    print(f"\nRedshift coverage:")
    print(f"  Total with redshift > 0: {redshift_with_data}")
    print(f"  AT redshifts: {(catalog_sorted['AT_Redshift'] > 0).sum()}")
    print(f"  MSAEXP redshifts: {(catalog_sorted['MSAEXP_Redshift'] > 0).sum()}")
    print(f"  LiMe redshifts: {(catalog_sorted['LiMe_Redshift'] > 0).sum()}")
    print(f"  MARZ redshifts: {(catalog_sorted['MARZ_Redshift'] > 0).sum()}")
    print(f"  CIGALE redshifts: {(catalog_sorted['Cigale_Redshift'] > 0).sum()}")
    print(f"  BAGPIPES redshifts: {(catalog_sorted['BAGPIPES_Redshift'] > 0).sum()}")
    #print(f"  Redshift mode calculated: {(catalog_sorted['Redshift_Mode'] > 0).sum()}")
    print(f"{'='*50}\n")
    
    return catalog_sorted


def main():
    """
    Main function to process data from multiple fields (UDS, EGS, COSMOS).
    """
    # Define data for each field
    fields_data = {
        "UDS": {
            "google_drive_url": "https://drive.google.com/uc?id=1DW5MsjzxYKi4bsNPh8WF8VUIhvblRgha",
            "output_folder": "data_UDS",
            "pointings": ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7']
        },
        "EGS": {
            "google_drive_url": "https://drive.google.com/uc?id=1mU0b3Glco7xHs4S9f2-lYjJUrvFPx8s9",
            "output_folder": "data_EGS",
            "pointings": ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7']
        },
        "COSMOS": {
            "google_drive_url": "https://drive.google.com/uc?id=1qx7-l1ZrgAr4OzOuEcQT-_BGF7bTOtLT",
            "output_folder": "data_COSMOS", 
            "pointings": ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7']
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
    print(f"\nCatalog saved as {default_catalog_name}")
    
    # Let user choose a custom name
    new_catalog_name = input("Choose a name for the output catalog (e.g., final_catalog_YourName.csv): ")
    if new_catalog_name:
        catalog_sorted.to_csv(new_catalog_name, index=False)
        print(f"Catalog saved as {new_catalog_name}")
    
    return catalog_sorted, new_catalog_name if new_catalog_name else default_catalog_name

if __name__ == "__main__":
    catalog_sorted, new_catalog_name = main()