# POMZ-UVP

Author: Jacob Cram
Email: jcram@umces.edu


Particle size data from the 2016-2017 POMZ Cruise to the Eastern Tropical North Pacific Oxygen Minimum Zone, from the RV Sikuliaq.

This project focuses specifically on one station, Station P2, which we surveyed over the course of a week Rachael Lekanoff ran the UVP, Gabrielle Rocap was the chief scientist in charge of the CTD ops. Jessica Pretty processed the UVP data. Andrew McDonnell lead the UVP processing operations. Data from the P16 cruise was provided by Andrew McDonnell and Jessica Pretty.
Flux data was provided by Clara Fuchsman.
CTD data was processed and provide by Al Devol
Almost all analysis herein was carried out by Jacob Cram

# To Do
  * Update particle math data to account for Weber and Burchfield's edits
  * Look through each data file to ensure proper commenting
  * Add description of files, or at least copy out of that other folder.
  * Move more stuff to archive
  * Possibly git filter-branch or equivalent out all of the raw files that aren't going in the analysis.
  * Share the code for processing the EK60 data, though not the actual raw EK60 data which are huge. Those, people can get from UNOLS directly.

# Description of Files
What follows is a description of the contents of this repository
Note: I'm using rstudio notebooks, which mean that bost of the .Rmd files, which are plain text and can be read by any text editor, also generate nb.html files, which are readable in a browser and which I am leaving.

## Files

## Directories
  * archive -- Analysis files that I used in early data exploration, model development, and so on that are not part of my primary analysis pipeline, as used in the manuscript. They are maintained for historical purposes
  * data -- raw data imported, includes UVP data from the ETNP-POMZ 2017 project, and staion 100 from the P16 project, as well as ctd data. There is also processed data from the EK60.
  * dataOut -- All of these files have been generated by running scripts in this repository, and then are re-used elsewhere
  * DisagDoc -- Math that describes how the Eularian version of the PRiSM model, against which the observed particle size spectra are compared, are calculated.
  * figures -- figures generated by scripts in this directory
  * manuscript -- Manuscript draft documents, some of which direclty import figures from the figures file.
  * renv -- package managment data
