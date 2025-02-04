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

def get_spectra_file_paths_and_ids(folder, pointing):
    spectra_dict = {}
    for filename in os.listdir(folder):
        if filename.startswith('.'):
            continue

        init_name = f'CAPERS_UDS_{pointing}_s'
        match = re.search(init_name + r'(\d{9})_(x1d_optext|s2d).fits', filename)
        if match:
            source_id = match.group(1)
            file_type = match.group(2)

            if source_id not in spectra_dict:
                spectra_dict[source_id] = {'x1d': None, 's2d': None}

            if file_type == 'x1d_optext':
                spectra_dict[source_id]['x1d'] = os.path.join(folder, filename)
            elif file_type == 's2d':
                spectra_dict[source_id]['s2d'] = os.path.join(folder, filename)

    return spectra_dict

def create_pointing_table(pointing, folder_path):
    folder = folder_path + f'/CAPERS_UDS_V0.1/{pointing}'
    spectra_files = get_spectra_file_paths_and_ids(folder, pointing)

    source_ids = []
    file1d = []
    file2d = []

    for source_id, files in spectra_files.items():
        source_ids.append(source_id)
        file1d.append(files['x1d'])
        file2d.append(files['s2d'])

    return source_ids, file1d, file2d

def read_redshift_catalogs(pointing):
    BAGPIPES = pd.read_csv(f'solutions/capers_UDS_bagpipes_{pointing}.cat', sep=r'\s+')
    AT = pd.read_csv(f'solutions/CAPERS_{pointing}_zspec.txt', sep=',')
    ID_AT = AT['file']
    for ids in range(len(ID_AT)):
        match = re.search(r's(\d+)_', ID_AT[ids])
        if match:
            extracted_number = match.group(1)
            ID_AT[ids] = extracted_number

    msaexp = pd.read_csv(f'solutions/capers_{pointing.lower()}_msaexp_zcat.csv')

    lime = pd.read_csv('solutions/CAPERS_UDS_V01_ASPECT_redshifts_v2.csv')
    ID_lime = lime['file']

    for ids in range(len(ID_lime)):
        match = re.search(r's(\d+)_', ID_lime[ids])
        if match:
            extracted_number = match.group(1)
            ID_lime[ids] = extracted_number

    marz_file = f'solutions/capers_uds_{pointing}_redshift.fits'
    marz = fits.open(marz_file)
    marz_data = marz[1].data 
    MARZ = pd.DataFrame(marz_data)
    file_MARZ = MARZ['Name']
    
    for ids in range(len(file_MARZ)):
        match = re.search(r's(\d+)', file_MARZ[ids])
        if match:
            extracted_number = match.group(1)
            file_MARZ[ids] = extracted_number

    CIGALE = pd.read_csv('solutions/CAPERS_v0.1_cigale_redshift_with_photometry.csv')
    CIGALE_ID = CIGALE['id']

    return BAGPIPES, AT, ID_AT, msaexp, lime, ID_lime, MARZ, file_MARZ, CIGALE, CIGALE_ID

def main():
    google_drive_url = "https://drive.google.com/uc?id=1nLwEX3edQI9ktHof9eAsH0z7lvUYFz7i"
    output_folder = "data"
    version_file = os.path.join(output_folder, "version.txt")
    expected_version = "v0.1"

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if os.path.exists(version_file):
        with open(version_file, 'r') as f:
            current_version = f.read().strip()
        if current_version == expected_version:
            print("Data is already downloaded and up-to-date.")
        else:
            print("Data version mismatch. Downloading the latest version.")
            download_and_extract_zip(google_drive_url, output_folder)
            with open(version_file, 'w') as f:
                f.write(expected_version)
    else:
        download_and_extract_zip(google_drive_url, output_folder)
        with open(version_file, 'w') as f:
            f.write(expected_version)

    folder_path = output_folder

    pointings = ['P1', 'P2', 'P3', 'P5']
    
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

    for pointing in pointings:
        source_ids, file1d, file2d = create_pointing_table(pointing, folder_path)
        total_source_ids.extend(source_ids)
        total_file1d.extend(file1d)
        total_file2d.extend(file2d)
        total_pointing.extend([pointing] * len(source_ids))

        BAGPIPES, AT, ID_AT, msaexp, lime, ID_lime, MARZ, file_MARZ, CIGALE, CIGALE_ID = read_redshift_catalogs(pointing)
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

    galaxies = np.array(total_source_ids)
    redshift_bp = np.zeros(len(galaxies))
    redshift_at = np.zeros(len(galaxies))
    redshift_msaexp = np.zeros(len(galaxies))
    redshift_lime = np.zeros(len(galaxies))
    redshift_marz = np.zeros(len(galaxies))
    redshifts_mode = np.zeros(len(galaxies))
    redshifts = np.zeros(len(galaxies))
    redshift_cigale = np.zeros(len(galaxies))

    for i in range(len(galaxies)):
        for BAGPIPES in all_bagpipes:
            for j in range(len(BAGPIPES['ID'])):
                if float(BAGPIPES['ID'][j]) == float(galaxies[i]):
                    redshift_bp[i] = BAGPIPES['z_bp'][j]

        for ID_AT, AT in zip(all_id_at, all_at):
            for m in range(len(ID_AT)):
                if float(ID_AT[m]) == float(galaxies[i]):
                    redshift_at[i] = AT['z'][m]

        for msaexp in all_msaexp:
            for n in range(len(msaexp['msa_id'])):
                if float(msaexp['msa_id'][n]) == float(galaxies[i]):
                    redshift_msaexp[i] = msaexp['z_spec'][n]

        for file_MARZ, MARZ in zip(all_file_marz, all_marz):
            for jj in range(len(file_MARZ)):
                if float(file_MARZ[jj]) == float(galaxies[i]):
                    redshift_marz[i] = MARZ['AutoZ'][jj]

        for ID_lime, lime in zip(all_id_lime, all_lime):
            for mm in range(len(ID_lime)):
                if float(ID_lime[mm]) == float(galaxies[i]):
                    redshift_lime[i] = lime['z_manual'][mm]

        for CIGALE_ID, CIGALE in zip(all_cigale_id, all_cigale):
            for kk in range(len(CIGALE_ID)):
                if float(CIGALE_ID[kk]) == float(galaxies[i]):
                    redshift_cigale[i] = CIGALE['best.universe.redshift'][kk]

    spectra_files_1d = np.array(total_file1d)
    spectra_files_2d = np.array(total_file2d)

    for i in range(len(galaxies)):
        z1 = round(redshift_at[i], 2)
        z3 = round(redshift_bp[i], 2)
        z4 = round(redshift_msaexp[i], 2)
        z5 = round(redshift_marz[i], 2)
        z6 = round(redshift_lime[i], 2)
        z = [z1, z3, z4, z5, z6]
        
        redshifts_mode[i] = scipy_mode(z)[0]
        
    catalog = pd.DataFrame({
        "Galaxy": galaxies,
        "1d_Spectrum": spectra_files_1d,
        "2d_Spectrum": spectra_files_2d,
        "Redshift": redshifts,
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

    catalog["Redshift"] = catalog["Redshift"].apply(lambda x: 0 if pd.isna(x) or x <= 0 else x)
    catalog_filtered = catalog[catalog["1d_Spectrum"].notnull() & (catalog["1d_Spectrum"] != "")]
    
    catalog_sorted1 = catalog_filtered.sort_values(by="Galaxy").reset_index(drop=True)
    
    visualize = input("Do you want to visualize sources within a specific ID range or a specific pointing? (range/pointing/no): ").strip().lower()
    if visualize == 'range':
        min_id = int(input("Enter the minimum ID: "))
        max_id = int(input("Enter the maximum ID: "))
        catalog_sorted = catalog_sorted1[(catalog_sorted1['Galaxy'].astype(int) >= min_id) & (catalog_sorted1['Galaxy'].astype(int) <= max_id)]
    elif visualize == 'pointing':
        chosen_pointing = input("Enter the pointing (P1, P2, P3, P5): ").strip().upper()
        catalog_sorted = catalog_sorted1[catalog_sorted1['Pointing'] == chosen_pointing]

    catalog_name = 'initial_catalog_total.csv'
    catalog_sorted.to_csv(catalog_name, index=False)
    
    new_catalog_name = input("Choose a name for the output catalog (remember the extension!): (e.g., final_pointing_YourName.csv) ")
    catalog_sorted.to_csv(new_catalog_name, index=False)
  
    
    return catalog_sorted, new_catalog_name

if __name__ == "__main__":
    catalog_sorted, new_catalog_name = main()
