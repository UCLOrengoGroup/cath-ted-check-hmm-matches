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
For domain visualization, please refer to the Jupyter notebook: `domain_visualization.ipynb` 

Use `vis_domain(df, index)` to visualize the domain boundaries in the human subset dataframe (df) from a specified row (index) in the Jupyter notebook using a PyMOL viewer.

