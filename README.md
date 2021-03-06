# POMZ-ETNP-UVP-2017

Author: Jacob Cram
Email: jcram@umces.edu


Particle size data from the 2016-2017 POMZ Cruise to the Eastern Tropical North Pacific Oxygen Minimum Zone, from the RV Sikuliaq.

This fork is for sharing with the public and avoids sharing data files that are not actually part of the project.

This project focuses specifically on one station, Station P2, which we surveyed over the course of a week Rachael Lekanoff ran the UVP, Gabrielle Rocap was the chief scientist in charge of the CTD ops. Jessica Pretty processed the UVP data. Andrew McDonnell lead the UVP processing operations. Data from the P16 cruise was provided by Andrew McDonnell and Jessica Pretty.
Flux data was provided by Clara Fuchsman.
CTD data was processed and provide by Al Devol
Almost all analysis herein was carried out by Jacob Cram

# To do:
 * Update disaggregation math. Make sure I am using the correct math.
 * Then update the readme.


# Description of Files
What follows is a description of the contents of this repository
Note: I'm using rstudio notebooks, which mean that bost of the .Rmd files, which are plain text and can be read by any text editor, also generate nb.html files, which are readable in a browser and which I am leaving.

## Files

# Description of Files

## Notes
This document exists to describe what the other files in this project do.

(Archive) means that the file isn't part of the main analyiss and can be moved to some other folder.

## Files in Main Directory

### Files that generate figures

  * GenerateFigures.Rmd -- Most of the figures that go in the manuscript, and some data analysis reported in the main text are in this file.

  * WaterMassEvans2 -- The water mass analyisis, creates the watermass figure

  * P16IntegratedPlots -- Figures about the oxic control station from the P16 transect in 2016

  * ExamineCTDProfiles -- Generates figure 1 with all of the CTD data

### Files that process data

  * UVP-2017-proc.R -- Run the main data processing

  * SmoothsAndFluxRevisited.R -- Optimizes fit of UVP data to trap data and saves output. (Keep)

### Libraries of R functions

  * ModelStuff.R -- The code for running the disaggregation model

  * UVP_2017_library.r -- Most of the functions for processing my data live here, except for the main disaggregation model functions.

  * bring_in_ctd.R -- Functions for loading CTD data (keep)

### Helper Files

  * references.bib -- Bibiography for manuscript

  * .gitignore -- tell git not to save silly stuff

  * README.md -- This document


## Directories
  * archive -- Analysis files that I used in early data exploration, model development, and so on that are not part of my primary analysis pipeline, as used in the manuscript. They are maintained for historical purposes
  * data -- raw data imported, includes UVP data from the ETNP-POMZ 2017 project, and staion 100 from the P16 project, as well as ctd data. There is also processed data from the EK60.
  * dataOut -- All of these files have been generated by running scripts in this repository, and then are re-used elsewhere
  * DisagDoc -- Math that describes how the Eularian version of the PRiSM model, against which the observed particle size spectra are compared, are calculated.
  * figures -- figures generated by scripts in this directory
  * manuscript -- Manuscript draft documents, some of which direclty import figures from the figures file.
  * renv -- package managment data

## Sub Structure of Directories

  * Archive
    * figures/ -- Figures that are not in the final analysis
    * WaterMassAnalysis -- an old version of water mass analysis (Archive)
    * DisagAttemptBivariateSmooth -- Some experiments with disaggreagation (Archive)
    * DisagModelTests -- Some experiments with disaggreagation (Archive)
    * ClaraFig5- Makes one figure that clara needed for a grant proposal (Archive)
    * Explore_C_r -- Uses the disaggregation model to explore how remineralization rate varies, but is complicated by active transport (Archive)
    * TestModelPipeline.r -- When I have problems with my analyis and model pipeline, I used some of this code for debbugging. (Archive)
    * P16S100CTDExplore -- An initial look at the oxic P16 station 100 ctd data. (Archive)
    * UVPGraveyard -- Graveyards are a collection of code and functions I am no longer using (Archive)
    * FrequentistParticleSmoothing -- An example of the smoothing approach that I am using that I shared with Dong Liang so he could make a baysean example of this same process. (Archive)
    * Normalize_UVP_Flux.Rmd -- An earlier version of the UVP flux normalization than SmoothsAndFluxRevisited.R (Archive)
    * DebugPSDGam.R -- Some debugging (Archive)
    * StackOverflowNb.Rmd -- A call for help online (Archive)
    * Zooplankton Sequestration Flux -- An investegation into how active transport relates to transfer efficiency. Interesting, but didn't make the manuscript (Archive)
    * DisagAttemptAlldgedgeAG -- Some experiments with disaggreagation (Archive). Used Alldrege estemates of fractal dimension for the first time here. (Archive)
    * GamVsGlmConfedenceQuestion -- A call to stack overflow for help. The answer was that I can use negative binomial distributions in gam functions but not glm functions, even if there are no smooths (Archive)
    * DisagAttemptNonPowerLaw -- Early version of disaggregation Model. (Archive)
    * Notes 24Sep2020 -- Some very isolated notes (Archive)
    * StackOverflowQuestions -- More calls for help online (Archive)
    * Remin_Library.R -- Early version of the remineralization model. Current versions are all in ModelStuff.R (Archive)
    * Explore_Big_And_Small -- Early look into large vs small particles. Size cutoff was I think 53 microns. (Archive) 
    * scratch.R -- playing around with stuff, who even knows what I was doing (Archive)
    * Explore_Quantiles_And_Size_Cutoff -- an early look into whether there were ever so few large particles that I could show there were more zeros than we would expect from chance. There weren't. (Archive)
    * Explore_Processed_UVP -- My earlist look at the UVP data (Archive)
    * SmoothDevelopment -- For figuring out how to properly smooth the data for some analyses (Archive)
    * Notes -- Mostly about some presentation I was writing once (Archive)
    * Notes24Sep2020.Rmd -- Another single session notes file. Essentially equivalent to a piece of scratch paper on my desk.
    * ExplorePAR -- Looking at photosynthetically available radiation.
    * MoreGamSmooth -- Jacob and Klaus Hubert debug an odd problem with the smooths -- it turns out it is a bad idea to try to average over profiles that don't extend through the depth interval
    * StackOverFlowGolorGradientLables -- A queston I needed to ask online to make sure the 0 and 24 both showed up on the color bars
    * CTDUVPDates - Confirming which casts happened when
    * Identify_Files_To_Remove_From_History.R -- Before I got permission to share them from Andrew McDonnell, I was going to remove unused data files. This finds them
    * FilesToRemoveFromHistory.txt -- A list of those files.
    * GenerateFigures_files -- At one point, I was using latex and these got created
      * {Contents of directory not listed}
  * MinistryOfTruth -- I had to remove some files that were not actually used in this project. They are all available upon request. I used the files in this directory to do this removal using `git filter-repo`.
    * FilesToRemoveFromHistory.txt -- The files that I removed
    * Identify_Files_to_Remove_From_History.R -- Goes into the UVP metadata and identifies which files to remove. I added a file containing all of the P16 data to the FilesToRemoveFromHistory.txt file as well. 
  * data
    * uvpdata -- the particle abundance data for all stations on the P16 cruise
      * export_detailed_20190304_23_14_Export_metadata_summary.tsv -- Metadata file
      * export_detailed_20190304_23_14_PAR_stn_XXX.tsv -- Particle data for one cast
      * export_detailed_20190304_23_14_ZOO_stn_XXX.tsv -- Zooplankton data for one cast (not used in this project)
    * P162016 -- UVP data for P16
      * {Directory Contents not listed here}
    * p16ctd -- CTD data for P16
      * {Directory Contents not listed here}
    * p16uvp -- UVP data for the only station from which I sampled
      * {Directory Contents not listed here}
    * WODatabaseP16 -- Chlorophyl data corresponding to our station of interest
      * {Directory Contents not listed here}
    * A20021852020182.L3m_CU_CHL_chlor_a_9km.nc -- Don't remember
    * A20190012019365.L3m_YR_CHL_chlor_a_9km.nc -- Don't remember
    * backscatter_table_go7.csv -- Processed EK60 backscatter data for station P2 -- Files used to generate this file are available upon request/when I get around to cleaning them up. Original data files are available from UNOLS and would never fit on a github repository.
    * lno11412-sup-0001-supinfo02.csv - Data from Natalya Evans's 2020 paper
    * POMZ2017Flux_CF.xlsx - Fuchsman's flux calculations
    * skq201617s-castinfo.csv - Cast Metadat
    * skq201617s-combinedctd.csv -- CTD data
    * SKQ201617S_CTD_Profile_CF.xlsx -- Another version of the CTD data
    * SKQ201617S_CTD_Profile.csv -- Another version of the CTD data
    * SKQ_26_27.mat -- Water mass data directly from Natalya Evans
    * StP2_2017_out2.xlsx --Don't remember
    * StP2_2017_out_v2b.xlsx -- Don't remember
    * StP2_2017_out_v2.xlsx -- Don't remember
    * StP2_2017_out.xlsx -- Don't remember
  * dataOut
    * binned_DepthSummary.csv -- Processed data, binned by depth, one row per depth, generated by UVP-2017-proc.R
    * binned_EachSize.csv -- Processed data, binned by depth, one row per particle size * depth, generated by UVP-2017-proc.R
    * unbinned_DepthSummary.csv -- Processed data, not depth binned, one row per depth, generated by UVP-2017-proc.R
    * unbinned_EachSize.csv -- Processed data, not depth binned, one row per particle size * depth, generated by UVP-2017-proc.R
    * fluxMS_distilled.csv -- Flux extemates from trap data. Created by another project. Details and code available upon request.
    * archive -- Files that this project generated at one point but that aren't current.
      * {Contents of this directory are not listed}
  * DisagDoc -- Math supporting the eularin disaggregation model -- not current
    * {Contents of this directory are not listed}
  * figures -- svg and png versions of all main and supplemental figures
    * "*.svg" -- same as png files, but different extension
    * AllParticleSizes.png -- Figure S4
    * BigVsSmall.png -- Figure S8
    * CombinedP16S100Info.png -- Figure S1
    * CombinedP2Info.png -- Figure 1
    * DisagExample.png -- Figure S11
    * ETNP.png -- Part of Figure 1, from Clara Fuchsman, developed by Hilary Palevsky
    * evans16p5N107W.png -- Figure S2
    * ExamplePSD163m.png -- Figure S5
    * FittedFlux.png -- Figure 3
    * FluxDeepDive.png -- Figure 5
    * FluxGamPlot.png -- Figure S6
    * FluxSizeShift.png -- Figure 7
    * OSMSGamPlot.png -- Figure S12 (OSMS was later renamed DFM)
    * P16FluxRelate.png -- Figure S9
    * P16MapManualExport.png -- A part of figure S1
    * P16Map.png -- A part of figure S1 (or maybe its the other P16MapMnaualExport.png file, I'm not sure)
    * ParticlesAndPSD_ETNPVsP16.png -- Figure S7
    * ParticlesPSDMany.png-- Figure 4
    * stationP2_EK60_18kOnly.png -- Figure 2
    * stationP2_EK60_go7.png -- Figure S3
    * WBModelValidation_Oxic.png -- Figure S10
    * WBModelValidation.png -- Figure 7
  * manuscript -- Draft manuscripts
    * AGU Template Documents -- Some of which I used, all of which I looked at.
    * AttemptsAtLatexManuscript -- I had a really rough couple of days where I tried to do this in latex. The main problem was reference management. In the end I went with word.
    * OldVersions -- Earlier manuscript versions
    * ParticleSizePOMZ2017_11May2021.docs -- Main manuscript draft
    * Cram_Supplement_May2021.docs -- Supplemental figures and text caption.


