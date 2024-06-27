# Check TED/HMM matches

Check domain boundaries between HMM matches and TED consensus predictions.

The incoming data is a subset of all HMM/TED domains corresponding to UniProt accessions from human.

Currently the processing is done in a Jupyter notebook which:
* filters the incoming data based on minimum domain length 
* creates bins according to overlap and evalue
* randomly selects 10 entries per bin  

### Visualizing the domain boundaries

   ```
   python -m venv myenv
   source myenv/bin/activate
   pip install -r requirements.txt
   ```

For domain visualization and analysis, please refer to the Jupyter notebook: `domain_visualization.ipynb` 

The notebook uses a custom function `visualize_domains_and_plot(df)` which provides an interactive visualization of domain boundaries and related metrics. This function does the following:

1. 3D Protein Structure Visualization:
   - Displays the 3D structure of a selected protein using PyMOL.
   - Colors different regions of the protein to indicate:
     - HMM-only domains (blue)
     - TED-only domains (red)
     - Overlapping HMM and TED domains (purple)
     - Regions without domain predictions (grey)

2. Interactive Scatter Plot:
   - Creates a scatter plot of HMM E-value vs. Overlap Percentage for all proteins in the dataset.
   - The x-axis (HMM E-value) uses a logarithmic scale.
   - Each point represents a protein, with the currently selected protein highlighted in red.

3. User Controls:
   - A slider and navigation buttons to move through the dataset.
   - A dropdown menu to sort the data by Index, HMM E-value, or Overlap Percentage.



