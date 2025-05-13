import os
import re
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.stats import mode as scipy_mode
import gdown
import zipfile

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
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Path to catalogs for this field (within the script directory)
    catalog_path = os.path.join(script_dir, "solutions", field)
    
    # BAGPIPES catalog
    if field == "UDS":
        BAGPIPES = pd.read_csv(os.path.join(catalog_path, f'capers_{field}_bagpipes_{pointing}.cat'), sep=r'\s+')
    elif field == "EGS":
        BAGPIPES = pd.read_csv(os.path.join(catalog_path, 'bagpipes_zfits.cat'), sep=r'\s+')
        print(BAGPIPES)
    # AT catalog
    AT = pd.read_csv(os.path.join(catalog_path, f'CAPERS_{pointing}_zspec.txt'), sep=',')
    ID_AT = AT['file']
    for ids in range(len(ID_AT)):
        match = re.search(r's(\d+)_', ID_AT[ids])
        if match:
            extracted_number = match.group(1)
            ID_AT[ids] = extracted_number

    # MSAEXP catalog
    if field == "UDS":
        msaexp = pd.read_csv(os.path.join(catalog_path, f'capers_{pointing.lower()}_msaexp_zcat.csv'))
    elif field == "EGS":
        msaexp = pd.read_csv(os.path.join(catalog_path, f'capers_{field.lower()}_msaexp_zcat.csv'))

    # LiMe catalog
    if field == "UDS":
        lime = pd.read_csv(os.path.join(catalog_path, f'CAPERS_{field}_V01_ASPECT_redshifts_v2.csv'))
        ID_lime = lime['file']
        for ids in range(len(ID_lime)):
            match = re.search(r's(\d{1,9})_', str(ID_lime[ids]))
            if match:
                extracted_number = match.group(1)
                ID_lime[ids] = extracted_number
    elif field == "EGS":
        lime = pd.read_csv(os.path.join(catalog_path, f'CAPERS_{field}_V0.2_lime_redshifts_v2.txt'), sep='\s+')
        ID_lime = lime['optext']
        for ids in range(len(ID_lime)):
            match = re.search(r's(\d{1,9})_', str(ID_lime[ids]))
            if match:
                extracted_number = match.group(1)
                ID_lime[ids] = extracted_number
    
    # MARZ catalog
    if field == "UDS":
        marz_file = os.path.join(catalog_path, f'capers_{field.lower()}_{pointing.lower()}_redshift.fits')
    elif field == "EGS":
        marz_file = os.path.join(catalog_path, f'capers_{field.lower()}_{pointing.lower()}_marz_result.fits')
    marz = fits.open(marz_file)
    marz_data = marz[1].data 
    MARZ = pd.DataFrame(marz_data)
    file_MARZ = MARZ['Name']
    for ids in range(len(file_MARZ)):
        match = re.search(r's(\d+)', file_MARZ[ids])
        if match:
            extracted_number = match.group(1)
            file_MARZ[ids] = extracted_number

    # CIGALE catalog
    if field == "UDS":
        CIGALE = pd.read_csv(os.path.join(catalog_path, f'CAPERS_v0.1_cigale_redshift_with_photometry.csv'))
    elif field == "EGS":
        CIGALE = pd.read_csv(os.path.join(catalog_path, f'redshift-cigale-v2.csv'))
    CIGALE_ID = CIGALE['id']

    # Photometric catalog
    if field == "UDS":
        photometric_catalog = pd.read_csv(os.path.join(catalog_path, f'CAPERS_{field}_master_yield_v2.1_actual.csv'))
    elif field == "EGS":
        photometric_catalog = pd.read_csv(os.path.join(catalog_path, f'CAPERS_{field}_master_yield_v2.csv'))
    photo_z = photometric_catalog[['ID', 'z_UNICORN', 'ra', 'dec']]
    photo_z_ID = photo_z['ID']

    return BAGPIPES, AT, ID_AT, msaexp, lime, ID_lime, MARZ, file_MARZ, CIGALE, CIGALE_ID, photo_z, photo_z_ID

def process_field(field, google_drive_url, output_folder, pointings):
    """
    Process data for a specific field.
    
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
    
    # Create empty lists to store data
    total_source_ids = []
    total_file1d = []
    total_file2d = []
    total_pointing = []

    all_bagpipes = []
    all_at = []
    all_id_at = []
    all_msaexp = []
    all_lime = []
    all_id_lime = []
    all_marz = []
    all_file_marz = []
    all_cigale = []
    all_cigale_id = []
    all_photo_z = []
    all_photo_z_id = []

    # Process each pointing for this field
    for pointing in pointings:
        print(f"Processing pointing {pointing} for {field}")
        
        # Get source IDs and file paths
        source_ids, file1d, file2d = create_pointing_table(pointing, output_folder, field)
        total_source_ids.extend(source_ids)
        total_file1d.extend(file1d)
        total_file2d.extend(file2d)
        total_pointing.extend([pointing] * len(source_ids))

        # Read redshift catalogs
        BAGPIPES, AT, ID_AT, msaexp, lime, ID_lime, MARZ, file_MARZ, CIGALE, CIGALE_ID, photo_z, photo_z_ID = read_redshift_catalogs(pointing, field)
        all_bagpipes.append(BAGPIPES)
        all_at.append(AT)
        all_id_at.append(ID_AT)
        all_msaexp.append(msaexp)
        all_lime.append(lime)
        all_id_lime.append(ID_lime)
        all_marz.append(MARZ)
        all_file_marz.append(file_MARZ)
        all_cigale.append(CIGALE)
        all_cigale_id.append(CIGALE_ID)
        all_photo_z.append(photo_z)
        all_photo_z_id.append(photo_z_ID)

    # Process the collected data
    galaxies = np.array(total_source_ids)
    redshift_bp = np.zeros(len(galaxies))
    redshift_at = np.zeros(len(galaxies))
    redshift_msaexp = np.zeros(len(galaxies))
    redshift_lime = np.zeros(len(galaxies))
    redshift_marz = np.zeros(len(galaxies))
    redshifts_mode = np.zeros(len(galaxies))
    redshifts = np.zeros(len(galaxies))
    redshift_cigale = np.zeros(len(galaxies))
    photo_z_values = np.zeros(len(galaxies))
    ra_values = np.zeros(len(galaxies))
    dec_values = np.zeros(len(galaxies))

    # Match galaxies with redshift values
    for i in range(len(galaxies)):
        # Process BAGPIPES data
        for BAGPIPES in all_bagpipes:
            for j in range(len(BAGPIPES['ID'])):
                if float(BAGPIPES['ID'][j]) == float(galaxies[i]):
                    redshift_bp[i] = BAGPIPES['z_bp'][j]

        # Process AT data
        for ID_AT, AT in zip(all_id_at, all_at):
            for m in range(len(ID_AT)):
                if float(ID_AT[m]) == float(galaxies[i]):
                    redshift_at[i] = AT['z'][m]

        # Process MSAEXP data
        for msaexp in all_msaexp:
            for n in range(len(msaexp['msa_id'])):
                if float(msaexp['msa_id'][n]) == float(galaxies[i]):
                    redshift_msaexp[i] = msaexp['z_spec'][n]

        # Process MARZ data
        for file_MARZ, MARZ in zip(all_file_marz, all_marz):
            for jj in range(len(file_MARZ)):
                if float(file_MARZ[jj]) == float(galaxies[i]):
                    redshift_marz[i] = MARZ['AutoZ'][jj]

        # Process LiMe data
        for ID_lime, lime in zip(all_id_lime, all_lime):
            for mm in range(len(ID_lime)):
                if float(ID_lime[mm]) == float(galaxies[i]):
                    redshift_lime[i] = lime['z_manual'][mm]

        # Process CIGALE data
        for CIGALE_ID, CIGALE in zip(all_cigale_id, all_cigale):
            for kk in range(len(CIGALE_ID)):
                if float(CIGALE_ID[kk]) == float(galaxies[i]):
                    redshift_cigale[i] = CIGALE['best.universe.redshift'][kk]

        # Process photometric redshifts
        for photo_z_id, photo_z in zip(all_photo_z_id, all_photo_z):
            for ll in range(len(photo_z_id)):
                if float(photo_z_id[ll]) == float(galaxies[i]):
                    photo_z_values[i] = photo_z['z_UNICORN'][ll]
                    ra_values[i] = photo_z['ra'][ll]
                    dec_values[i] = photo_z['dec'][ll]

    # Convert to numpy arrays
    spectra_files_1d = np.array(total_file1d)
    spectra_files_2d = np.array(total_file2d)

    # Calculate redshift mode
    for i in range(len(galaxies)):
        z1 = round(redshift_at[i], 2)
        z3 = round(redshift_bp[i], 2)
        z4 = round(redshift_msaexp[i], 2)
        z5 = round(redshift_marz[i], 2)
        z6 = round(redshift_lime[i], 2)
        z = [z1, z3, z4, z5, z6]
        
        redshifts_mode[i] = scipy_mode(z)[0]
    
    # Create catalog DataFrame
    catalog = pd.DataFrame({
        "Galaxy": galaxies,
        "1d_Spectrum": spectra_files_1d,
        "2d_Spectrum": spectra_files_2d,
        "RA": ra_values,
        "DEC": dec_values,
        "photo_Redshift": photo_z_values,
        "Redshift": redshifts,
        "Field": [field] * len(galaxies),  # Add field name to catalog
        "Pointing": total_pointing,
        "Flag": [""]*len(galaxies),
        "AT_Redshift": redshift_at,
        "MSAEXP_Redshift": redshift_msaexp,
        "LiMe_Redshift": redshift_lime,
        "MARZ_Redshift": redshift_marz,
        "Cigale_Redshift": redshift_cigale,
        "BAGPIPES_Redshift": redshift_bp,
        "Redshift_Mode": redshifts_mode,
        "[OIII]+Hβ": [""]*len(galaxies),
        "Hα": [""]*len(galaxies),
        "Lyα": [""]*len(galaxies),
        "Paα": [""]*len(galaxies),
        "Lyman Break": [""]*len(galaxies),
        "Balmer Break": [""]*len(galaxies),
        "Comments": [""]*len(galaxies)
    })

    # Clean up and filter the catalog
    catalog["Redshift"] = catalog["Redshift"].apply(lambda x: 0 if pd.isna(x) or x <= 0 else x)
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
            "google_drive_url": "https://drive.google.com/uc?id=1-lBbb9AT9HcBaKe067YHTjIQfOU6FazL",
            "output_folder": "data_EGS",
            "pointings": ['P1', 'P2', 'P3', 'P4', 'P5', 'P7']  # Update with actual EGS pointings
        },
        "COSMOS": {
            "google_drive_url": "YOUR_COSMOS_GOOGLE_DRIVE_LINK_HERE",
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
    visualize = input("Do you want to visualize sources within a specific ID range, a specific pointing, or a specific field? (range/pointing/field/all): ").strip().lower()
    
    if visualize == 'range':
        min_id = int(input("Enter the minimum ID: "))
        max_id = int(input("Enter the maximum ID: "))
        catalog_sorted = catalog_sorted[(catalog_sorted['Galaxy'].astype(int) >= min_id) & (catalog_sorted['Galaxy'].astype(int) <= max_id)]
    elif visualize == 'pointing':
        chosen_pointing = input("Enter the pointing (e.g., P1, P2, P3, P5): ").strip().upper()
        catalog_sorted = catalog_sorted[catalog_sorted['Pointing'] == chosen_pointing]
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