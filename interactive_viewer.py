def interactive_spectrum_viewer(index=0):
    """
    View spectra interactively following the catalog order.
    """

    codes = ["AT", "MSAEXP", "LiMe", "MARZ", "Cigale", "BAGPIPES"]
    
    def get_Redshift_Mode(index):
        z_values = [catalog_filtered.iloc[index][f'{code}_Redshift'] for code in codes]
        rounded_values = np.array([round(val, 2) for val in z_values])
        return float(scipy_mode(rounded_values)[0])
    
    redshift_slider = FloatSlider(value=get_Redshift_Mode(index), min=0, max=20, step=0.001, description="Redshift:")
    redshift_slider.style.handle_color = 'lightblue'
    galaxy_id_input = Text(value=str(catalog_filtered.iloc[index]["Galaxy"]), layout=Layout(width='250px'), description="ID:", placeholder="Enter ID", disabled=True)
    em_line_button = Button(description="+ lines", button_style='danger', layout=Layout(width='60px'))
    attention_button = Button(description="Review needed", button_style='danger')
    em_line_check = False
    
    def toggle_em_line_state(button):
        nonlocal em_line_check
        em_line_check = not em_line_check
        button.button_style = "success" if em_line_check else "danger"
        update_viewer()

    def mark_attention(button):
        catalog_filtered.loc[index, "Redshift"] = -2
        catalog_filtered.loc[index, "Comments"] = "Uncertain solution, flagged for review"
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

    top_widgets = HBox([galaxy_id_input, redshift_slider, em_line_button, attention_button], layout=Layout(align_items='center', justify_content='center'))
    left_widgets = VBox([top_widgets, HBox([flag_dropdown, comments_box])], layout=Layout(width='50%', align_items='center', justify_content='center'))

    right_widgets = VBox(
        [HBox([prev_button, next_button, save_button], layout=Layout(justify_content='flex-end')),
         features_list],
        layout=Layout(width='35%', align_items='flex-end')
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

        styled_checkboxes = [
            HBox([HTML(f'<span style="color: {color};">{code} ({round(catalog_filtered[f"{code}_Redshift"].iloc[index], 2)})</span>'), checkbox])
            for code, color, checkbox in zip(codes, colorlist, checkboxes_redshift)
        ]

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
            redshift_slider.value = get_Redshift_Mode(index)
            update_viewer()

    def on_prev(_):
        nonlocal index
        if index > 0:
            index -= 1
            redshift_slider.value = get_Redshift_Mode(index)
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
    attention_button.on_click(mark_attention)

    update_viewer()

    controls = HBox([prev_button, next_button, save_button])
    display(VBox([output]))

if __name__ == "__main__":
    catalog_filtered = pd.read_csv("path_to_your_catalog.csv")
    new_catalog_name = "output_catalog.csv"
    interactive_spectrum_viewer()
