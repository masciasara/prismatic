# PRISMATIC Redshifting Program

<center><img width="537" alt="prismatic_logo" src="https://github.com/user-attachments/assets/be684848-9304-40fe-8420-851bfad3f844" /></center>
  
**PRISMATIC** is a Python-based tool designed to analyze and determine redshifts and redshift quality flags from **JWST NIRSpec prism spectra**. This program is optimized for handling low-resolution spectra obtained from the NIRSpec prism, in particular spectra from the **The CANDELS-Area Prism Epoch of Reionization Survey (CAPERS)** JWST program (GO 6368, PI Dickinson).


---

## Features  
- **Redshift Estimation**: Automatically includes best redshifts solution from several independent codes, obtainesd using spectral line identification and/or template fitting.  
- **JWST Optimization**: Tailored for NIRSpec prism spectra, ensuring compatibility with JWST data products.  
- **Ease of Use**: Intuitive user interface and simple configuration options for fast setup and analysis.   

---

## Installation  

To install PRISMATIC, ensure you have Python 3.8 or later installed. Then clone this repository and install the required dependencies:

```bash
git clone https://github.com/your_username/PRISMATIC.git  
cd PRISMATIC  
pip install -r requirements.txt
