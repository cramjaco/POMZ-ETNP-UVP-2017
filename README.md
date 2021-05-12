# POMZ-UVP

Author: Jacob Cram
Email: jcram@umces.edu


Particle size data from the 2016-2017 POMZ Cruise to the Eastern Tropical North Pacific Oxygen Minimum Zone, from the RV Sikuliaq.

This project focuses specifically on one station, Station P2, which we surveyed over the course of a week Rachael Lekanoff ran the UVP, Gabrielle Rocap was the chief scientist in charge of the CTD ops. Jessica Pretty processed the UVP data. Andrew McDonnell lead the UVP processing operations. Data from the P16 cruise was provided by Andrew McDonnell and Jessica Pretty.
Flux data was provided by Clara Fuchsman.
CTD data was processed and provide by Al Devol
Almost all analysis herein was carried out by Jacob Cram


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


