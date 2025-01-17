import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from astropy.io import fits
from ipywidgets import FloatSlider, Text, Dropdown, Textarea, Button, Checkbox, VBox, HBox, Label, Layout, Output, HTML
from scipy.stats import mode as scipy_mode

emission_lines = {"Lyα": 1215.67,
                  "HeII": 1640.42, 
                  "CIV": 1550.,
                  "CIII": 1909.,
                  "[OII](1)": 3727.104102,
                  "[OII](2)": 3729.887902,
                  "[NeIII]": 3869.873137,
                  "Hδ": 4102.922122,
                  "Hγ": 4341.719762,
                  "Hβ": 4862.73153,
                  "[OIII](1)": 4960.337588,
                  "[OIII](2)": 5008.283371,
                  "Hα": 6564.706817,
                  "[SII](1)": 6718.371995,
                  "[SII](2)": 6732.746127,
                  "[SIII]": 9533.841463,
                  "Paδ": 10052.25852,
                  "HeI": 10833.45512,
                  "Paγ": 10941.23211,
                  "[FeII]": 12570.38253,
                  "Paβ": 12821.7568,
                  "Paα": 18756.36916,
                  "Brβ": 26259.07141,
                  "Pfundβ": 37406.11053,
                  "Brα": 40534.349,
                  "Pfundα": 46551.12201
                 }

additional_emission_lines = {"OIII]": 1660.809,
                             "OIII]": 1666.15,
                             "MgII]": 2796.33262,
                             "MgII]": 2803.511683,
                             "[NeIII]": 3968.651338,
                             "[SIII]": 6313.875704,
                             "[NII]": 6549.933569,
                             "[NII]": 6585.353753,
                             "[SIII]": 9071.20845,
                             "Brδ": 19451.17493,
                             "Brγ": 21661.53048
                            }

def plot_spectrum(file, file_2d, redshift, galaxy_id, c_values, em_line_check):
    global catalog_filtered

    codes = ["AT", "MSAEXP", "LiMe", "MARZ", "Cigale", "BAGPIPES"]

    with fits.open(file) as hdul:
        data = hdul[1].data
        wavelength = data['WAVELENGTH'] * 10**4
        flux = data['FLUX']
        noise = data['FLUX_ERROR']

    with fits.open(file_2d) as hdul1:
        spectrum_2d = hdul1['SCI'].data

    wave_new = np.concatenate([wavelength, [wavelength[-1]+1]])
    y_pixels = np.arange(spectrum_2d.shape[0]+1)
    X, Y = np.meshgrid(wave_new, y_pixels)
    flux2d_not_nan = spectrum_2d[~np.isnan(spectrum_2d)]
    vmin = np.percentile(flux2d_not_nan, 15)
    vmax = np.percentile(flux2d_not_nan, 95)
    norm = Normalize(vmin=vmin, vmax=vmax)

    redshifts_dict = {code: {"redshift": catalog_filtered.loc[catalog_filtered['Galaxy'] == galaxy_id, f"{code}_Redshift"].iloc[0], "checked": c_values[idx]} for idx, code in enumerate(codes)}

    fig = plt.figure(figsize=(14, 7))
    gs = plt.GridSpec(2, 2, width_ratios=[2.5, 0.3], height_ratios=[0.5, 1])

    ax2d = fig.add_subplot(gs[0, 0])
    ax2d.set_title(f"Galaxy ID: {galaxy_id}", fontsize=12)
    pcm = ax2d.pcolormesh(X, Y, spectrum_2d, cmap='inferno', vmin=vmin, vmax=vmax, shading='auto')

    ax_hist = fig.add_subplot(gs[0, 1])
    pixel_values = spectrum_2d.flatten()
    pixel_values = pixel_values[~np.isnan(pixel_values)]
    cmap = plt.cm.inferno
    hist, bins = np.histogram(pixel_values, bins=20, range=(vmin, vmax))

    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    colors = cmap(norm(bin_centers))

    for left, right, h, color in zip(bins[:-1], bins[1:], hist, colors):
        ax_hist.barh(y=0.5 * (left + right), width=h, height=right - left, color=color, edgecolor='none')

    ax_hist.axis('off')
    ax_hist.set_ylim([vmin, vmax])
    ax_hist.set_xlim([0, 2])
    ax_hist.set_xlabel('Pixel Count')
    ax_hist.yaxis.tick_right()
    ax_hist.yaxis.set_label_position("right")

    ax1d = fig.add_subplot(gs[1, 0])
    ax1d.plot(wavelength, flux, c='#003049', label=f"z = {redshift:.3f}")
    ax1d.errorbar(wavelength, flux, yerr=noise, c='#669bbc', fmt='none', alpha=0.4)

    colors = ['#390099', '#006747', '#ff006e', '#8338ec', '#3a86ff', '#ffa600']
    offsets = [0.15, 0.3, 0.4, 0.55, 0.7, 0.85]

    for idx, (code, data) in enumerate(redshifts_dict.items()):
        if not data["checked"]:
            continue
        z = data["redshift"]
        color = colors[idx % len(colors)]
        shift = offsets[idx % len(offsets)]
        for line, rest_wavelength in emission_lines.items():
            observed_line = rest_wavelength * (1 + z)
            if np.nanmin(wavelength) < observed_line < np.nanmax(wavelength):
                ax1d.axvline(observed_line, color=color, ls=":", alpha=0.7, lw=0.8)
                ax1d.text(observed_line + 5, shift * np.nanmax(flux),
                          f"{line}", color=color, rotation=90,
                          verticalalignment="center", horizontalalignment="left", fontsize=6)

    for line, rest_wavelength in emission_lines.items():
        observed_line = rest_wavelength * (1 + redshift)
        if observed_line < np.nanmax(wavelength) and observed_line > np.nanmin(wavelength):
            ax1d.axvline(observed_line, color='k', ls=':', alpha=0.7, lw=0.8)
            ax1d.text(observed_line + 5, 0.6 * np.nanmax(flux), f"{line}", color='k',
                      rotation=90, verticalalignment='center', horizontalalignment='left', fontsize=8)
    
    if em_line_check:
        for line, rest_wavelength in additional_emission_lines.items():
            observed_line = rest_wavelength * (1 + redshift)
            if observed_line < np.nanmax(wavelength) and observed_line > np.nanmin(wavelength):
                ax1d.axvline(observed_line, color='grey', ls=':', alpha=0.7, lw=0.8)
                ax1d.text(observed_line + 5, 0.6 * np.nanmax(flux), f"{line}", color='grey',
                          rotation=90, verticalalignment='center', horizontalalignment='left', fontsize=8)        

    ax1d.set_xlabel("Wavelength [Å]")
    ax1d.set_ylabel("Flux [Jy]")
    ax1d.legend(loc="upper right", fontsize="small", ncol=2)
    ax1d.grid(c='lightgrey', ls=':')

    ax2d.sharex(ax1d)
    plt.tight_layout()
    plt.show()

    return ax_hist, pcm, fig

# Create an output widget
output = Output()

def interactive_spectrum_viewer(index=0):
    """
    View spectra interactively following the catalog order.
    """

    codes = ["AT", "MSAEXP", "LiMe", "MARZ", "Cigale", "BAGPIPES"]
    
    def get_mode_redshift(index):
        z_values = [catalog_filtered.iloc[index][f'{code}_Redshift'] for code in codes]
        rounded_values = np.array([round(val, 2) for val in z_values])
        return float(scipy_mode(rounded_values)[0])
    
    redshift_slider = FloatSlider(value=get_mode_redshift(index), min=0, max=20, step=0.001, description="Redshift:")
    redshift_slider.style.handle_color = 'lightblue'
    galaxy_id_input = Text(value=str(catalog_filtered.iloc[index]["Galaxy"]), description="ID:", placeholder="Enter ID")
    em_line_button = Button(description="+ lines", button_style='danger')
    em_line_check = False
    
    def toggle_em_line_state(button):
        nonlocal em_line_check
        em_line_check = not em_line_check
        button.button_style = "success" if em_line_check else "danger"
        update_viewer()

    flag_dropdown = Dropdown(options=["", "4", "3", "2", "1", "0", "9"], value=catalog_filtered.iloc[index]["Flag"], description="Flag:")
    comments_box = Textarea(value='', description='Comments:', placeholder='Enter any comments here...', layout=Layout(width='100%', height='50px'))

    next_button = Button(description="Next", button_style='primary')
    prev_button = Button(description="Previous", button_style='primary')
    save_button = Button(description="Save to CSV", button_style='warning')

    features = ["[OIII]+Hβ", "Hα", "Lyα", "Lyman Break", "Balmer Break"]
    checkboxes = [Checkbox(value=False, description=feature) for feature in features]
    left_column = VBox(checkboxes[:2])
    right_column = VBox(checkboxes[2:])
    checkbox_columns = HBox([left_column, right_column])
    features_list = VBox([Label("Select features:")] + [checkbox_columns])
    features_list.layout = Layout(border='1px solid lightgrey', padding='1px', width='100%')

    top_widgets = HBox([galaxy_id_input, redshift_slider, em_line_button], layout=Layout(align_items='center', justify_content='center'))
    left_widgets = VBox([top_widgets, HBox([flag_dropdown, comments_box])], layout=Layout(width='50%', align_items='center', justify_content='center'))

    right_widgets = VBox(
        [HBox([prev_button, next_button, save_button], layout=Layout(justify_content='flex-end')),
         features_list],
        layout=Layout(width='45%', align_items='flex-end')
    )

    ui = HBox([left_widgets, right_widgets], layout=Layout(width='100%', justify_content='space-between'))
    display(ui)

    def update_viewer():
        """
        Update the spectrum viewer for the current source.
        """
        galaxy = catalog_filtered.iloc[index]
        file = galaxy["1d_Spectrum"]
        file_2d = galaxy["2d_Spectrum"]
        redshift = redshift_slider.value
        flag = galaxy["Flag"]
        galaxy_id = galaxy["Galaxy"]
        comments = galaxy["Comments"]
        
        flag_dropdown.value = flag
        comments_box.value = comments
        galaxy_id_input.value = str(galaxy_id)
        
        for cb, feature in zip(checkboxes, features):
            cb.value = catalog_filtered.iloc[index][feature] == 1

        colorlist = ['#390099', '#006747', '#ff006e', '#8338ec', '#3a86ff', '#ffa600']
        checkboxes_redshift = [Checkbox(value=False, description='', style={'description_width': 'initial'}, layout={'width': 'auto'}) for _ in codes]

        global z_list
        try:
            z_list.close()
        except:
            pass

        # Now, create a list of HBox to pair each checkbox with its colored description
        styled_checkboxes = [
            HBox([HTML(f'<span style="color: {color};">{code} ({round(catalog_filtered[f"{code}_Redshift"].iloc[index], 2)})</span>'), checkbox])
            for code, color, checkbox in zip(codes, colorlist, checkboxes_redshift)
        ]

        # Create the grid layout for the checkboxes
        left_column = VBox(styled_checkboxes[:1])
        middle_1column = VBox(styled_checkboxes[1:2])
        middle_2column = VBox(styled_checkboxes[2:3])
        middle_3column = VBox(styled_checkboxes[3:4])
        middle_4column = VBox(styled_checkboxes[4:5])
        right_column = VBox(styled_checkboxes[5:])
        checkbox_columns_z = HBox([left_column, middle_1column, middle_2column, middle_3column, middle_4column, right_column])
        z_list = VBox([checkbox_columns_z])
        z_list.layout = Layout(border='1px solid lightgrey', padding='1px', width='100%')

        def on_solution_change(change):
            with output:
                output.clear_output(wait=True)
                ax_hist, pcm, fig = plot_spectrum(file, file_2d, redshift_slider.value, galaxy_id, [cb.value for cb in checkboxes_redshift], em_line_check)
                display(z_list)
                def on_hist_zoom(event):
                    new_vmin, new_vmax = ax_hist.get_ylim()
                    norm = pcm.norm
                    norm.vmin = new_vmin
                    norm.vmax = new_vmax
                    pcm.set_clim(vmin=new_vmin, vmax=new_vmax)
                    fig.canvas.draw_idle()

                ax_hist.callbacks.connect('ylim_changed', on_hist_zoom)

        for cb_z in checkboxes_redshift:
            cb_z.observe(on_solution_change, names='value')

        with output:
            output.clear_output(wait=True)
            ax_hist, pcm, fig = plot_spectrum(file, file_2d, redshift, galaxy_id, [cb.value for cb in checkboxes_redshift], em_line_check)
            display(z_list)
            def on_hist_zoom(event):
                new_vmin, new_vmax = ax_hist.get_ylim()
                norm = pcm.norm
                norm.vmin = new_vmin
                norm.vmax = new_vmax
                pcm.set_clim(vmin=new_vmin, vmax=new_vmax)
                fig.canvas.draw_idle()

            ax_hist.callbacks.connect('ylim_changed', on_hist_zoom)

    def on_next(_):
        nonlocal index
        if index < len(catalog_filtered) - 1:
            index += 1
            redshift_slider.value = get_mode_redshift(index)
            update_viewer()

    def on_prev(_):
        nonlocal index
        if index > 0:
            index -= 1
            redshift_slider.value = get_mode_redshift(index)
            update_viewer()

    def save_to_csv(_):
        catalog_filtered.loc[index, "Comments"] = comments_box.value
        for cb, feature in zip(checkboxes, features):
            catalog_filtered.loc[index, feature] = 1 if cb.value else 0
        catalog_filtered.to_csv(f"{new_catalog_name}", index=False)
        with output:
            print(f"Catalog saved to '{new_catalog_name}'.")

    def on_redshift_change(change):
        catalog_filtered.loc[index, "Redshift"] = change["new"]
        update_viewer()

    def on_flag_change(change):
        catalog_filtered.loc[index, "Flag"] = change["new"]

    redshift_slider.observe(on_redshift_change, names="value")
    flag_dropdown.observe(on_flag_change, names="value")
    next_button.on_click(on_next)
    prev_button.on_click(on_prev)
    save_button.on_click(save_to_csv)
    em_line_button.on_click(toggle_em_line_state)

    update_viewer()

    controls = HBox([prev_button, next_button, save_button])
    display(VBox([output]))

if __name__ == "__main__":
    catalog_filtered = pd.read_csv("path_to_your_catalog.csv")
    new_catalog_name = "output_catalog.csv"
    interactive_spectrum_viewer()