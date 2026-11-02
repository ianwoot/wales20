# wales20

Scripts and data to replicate Walker's Wales 20 mph analysis

You should be able to run the whole analysis by running do_analysis.R. The data files are too large to upload
to Github (e.g., 2017-2025 collision data ~48 MB) but do_analysis.R has details of where they are obtained from, 
and build_panel.R will recreate the data files. 

That said, you also shouldn't need to run the analysis: I've uploaded the "output" and "figures"
folders from after I ran the analysis, so what would be produced if you ran the scripts should already 
be there for inspection.

## Sources (downloaded under data/raw)
- STATS19 road safety open data, collisions and casualties 2017-2025,
  DfT (record-level, includes DfT severity-adjustment probabilities).
- TRA8904: motor vehicle traffic (vehicle kilometres) by local authority, DfT.
- TRA0413: pedal cycle traffic (vehicle kilometres) by region and country, DfT.
- ONS mid-year population estimates (MYEB2 back-series, current boundaries).
- ONS Open Geography: LAD -> county lookup; LAD boundaries (Dec 2024, BUC).
