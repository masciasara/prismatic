# PRISMATIC Redshifting Program

<center><img width="537" alt="PRISMATIC_logo" src="https://github.com/user-attachments/assets/be684848-9304-40fe-8420-851bfad3f844" /></center>
  
**PRISMATIC** is a Python-based tool designed to analyze and determine redshifts and quality flags from **JWST NIRSpec prism spectra**. This program is optimized for handling low-resolution spectra obtained from the NIRSpec prism, in particular spectra from the **The CANDELS-Area Prism Epoch of Reionization Survey (CAPERS)** JWST program (GO 6368, PI Dickinson).

---
<h2>Features</h2>
    <ol>
        <li>
            <strong>Redshift selection:</strong>
            <p>PRISMATIC allows users to visually inspect spectral data and choose the correct redshift by matching emission lines with observed spectral features.</p>
        </li>
      <li>
        <strong>Solution selection from independent codes:</strong>
          <p>The tool integrates outputs from various automatic codes (<strong>msaexp</strong>, <strong>bagpipes</strong>, <strong>Cigale</strong>, <strong>Marz</strong>, <strong>LiMe</strong>, ...), enabling users to compare suggested redshift solution from each code and select the most accurate one based on visual inspection and feature matching.</li>
        <li>
            <strong>Flag assignment:</strong>
            <p>Users can associate a confidence flag to each selected redshift, based on the following legend:</p>
            <table>
                <thead>
                    <tr>
                        <th>Flag</th>
                        <th>Description</th>
                    </tr>
                </thead>
                <tbody>
                    <tr>
                        <td>4</td>
                        <td>Very secure redshift. Several clear features matching.</td>
                    </tr>
                    <tr>
                        <td>3</td>
                        <td>Secure redshift. Various high and low confidence features matching.</td>
                    </tr>
                    <tr>
                        <td>2</td>
                        <td>Probable redshift. Various low confidence features.</td>
                    </tr>
                    <tr>
                        <td>1</td>
                        <td>Uncertain redshift. Single low confidence feature.</td>
                    </tr>
                    <tr>
                        <td>0</td>
                        <td>Unknown. No detection or absence of features.</td>
                    </tr>
                    <tr>
                        <td>9</td>
                        <td>Best guess for single high confidence feature.</td>
                    </tr>
                </tbody>
            </table>
        </li>
        <li>
            <strong>Review Needed:</strong>
            <p>If the source exhibits some features but the redshift is challenging to determine, users can click the "Review needed" button. This action will:</p>
            <ul>
                <li>Assign a redshift of <code>-2</code>.</li>
                <li>Add a comment: "Uncertain solution, flagged for review"</li>
            </ul>
        </li>
        <li>
            <strong>Common Spectral Features:</strong>
            <p>The tool includes an option to select and highlight the most common spectral features for quick reference.</p>
        </li>
        <li>
            <strong>Comments Section:</strong>
            <p>Users can leave additional comments or notes for each source.</p>
        </li>
      <li> <strong> Remember to save each source's information using the <code>Save to csv</code> button before moving to the next one.</li>
</strong></li>
    </ol>
    
---

# PRISMATIC Analysis Notebook

```bash
git clone https://github.com/masciasara/PRISMATIC.git  
cd PRISMATIC  
pip install -r requirements.txt
jupyter notebook PRISMATIC.ipynb
```

---
## Contributions

Contributions are welcome! If you have suggestions, bug reports, or requests, feel free to open an issue or contact me at sara.mascia@ista.ac.at.

---
## Acknowledgments

PRISMATIC was developed as part of the CAPERS project. 
