# wales20

Scripts and data to replicate Walker's Wales 20 mph analysis

You should be able to run the whole analysis by running do_analysis.R, which is a wrapper script to call
the various components of the analysis. The early stages of this wrapper script will extract and build 
the panel data from the raw data sources below. However, as some of the raw data files from DfT are very 
large, it is more likely you will want to just run the later stages of the analysis that are listed in 
do_analysis.R, such as analysis.R, perm_all.R, etc. These scripts run on the pre-prepared panel data, 
which are provided here in the file data/panel.rds

That said, you also shouldn't need to run the analysis: the "output" and "figures"
folders from running the analysis are provided, so what would be produced if you ran the scripts should 
already be there for inspection.

## Sources (downloaded under data/raw)
- STATS19 road safety open data, collisions and casualties 2017-2025,
  DfT (record-level, includes DfT severity-adjustment probabilities).
- TRA8904: motor vehicle traffic (vehicle kilometres) by local authority, DfT.
- TRA0413: pedal cycle traffic (vehicle kilometres) by region and country, DfT.
- ONS mid-year population estimates (MYEB2 back-series, current boundaries).
- ONS Open Geography: LAD -> county lookup; LAD boundaries (Dec 2024, BUC).
