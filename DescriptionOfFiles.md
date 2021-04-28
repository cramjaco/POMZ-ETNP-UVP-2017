# Description of Files

## Notes
This document exists to describe what the other files in this project do.

(Archive) means that the file isn't part of the main analyiss and can be moved to some other folder.

## Files in Main Directory

references.bib -- Bibiography for manuscript

GenerateFigures.Rmd -- Most of the figures that go in the manuscript, and some data analysis reported in the main text are in this file.

SmoothDevelopment -- For figuring out how to properly smooth the data for some analyses (Archive)

WaterMassEvans2 -- The water mass analyisis

P16IntegratedPlots -- Figures about the oxic control station from the P16 transect in 2016

ModelStuff.R -- The code for running the disaggregation model

UVP_2017_library.r -- Most of the functions for processing my data live here, except for the main disaggregation model functions.

SmoothsAndFluxRevisited.R -- Optimizes fit of UVP data to trap data and saves output. (Keep)

ExamineCTDProfiles -- Generates figure 1 with all of the CTD data

.gitignore -- tell git not to save silly stuff

bring_in_ctd.R -- Functions for loading CTD data (keep)

.Rprofile -- Functions that run every time I launch R no matter what

Notes -- Mostly about some presentation I was writing once

README.md -- I shoud make this look nicer and more informative.

UVP-2017-proc.R -- Run the main data processing

DescriptionOfFiles -- If I describe too thouroughly, I get infinate recursion. (Archive -- describe the keep files in the readme).


## Files in *archive* directory

WaterMassAnalysis -- an old version of water mass analysis (Archive)


DisagAttemptBivariateSmooth -- Some experiments with disaggreagation (Archive)

DisagModelTests -- Some experiments with disaggreagation (Archive)



ClaraFig5- Makes one figure that clara needed for a grant proposal (Archive)

Explore_C_r -- Uses the disaggregation model to explore how remineralization rate varies, but is complicated by active transport (Archive)

TestModelPipeline.r -- When I have problems with my analyis and model pipeline, I used some of this code for debbugging. (Archive)

P16S100CTDExplore -- An initial look at the oxic P16 station 100 ctd data. (Archive)

UVPGraveyard -- Graveyards are a collection of code and functions I am no longer using (Archive)

FrequentistParticleSmoothing -- An example of the smoothing approach that I am using that I shared with Dong Liang so he could make a baysean example of this same process. (Archive)

Normalize_UVP_Flux.Rmd -- An earlier version of the UVP flux normalization than SmoothsAndFluxRevisited.R (Archive)

DebugPSDGam.R -- Some debugging (Archive)

StackOverflowNb.Rmd -- A call for help online (Archive)

Zooplankton Sequestration Flux -- An investegation into how active transport relates to transfer efficiency. Interesting, but didn't make the manuscript (Archive)

DisagAttemptAlldgedgeAG -- Some experiments with disaggreagation (Archive). Used Alldrege estemates of fractal dimension for the first time here. (Archive)

GamVsGlmConfedenceQuestion -- A call to stack overflow for help. The answer was that I can use negative binomial distributions in gam functions but not glm functions, even if there are no smooths (Archive)

DisagAttemptNonPowerLaw -- Early version of disaggregation Model. (Archive)

Notes 24Sep2020 -- Some very isolated notes (Archive)

StackOverflowQuestions -- More calls for help online (Archive)

Remin_Library.R -- Early version of the remineralization model. Current versions are all in ModelStuff.R (Archive)

Explore_Big_And_Small -- Early look into large vs small particles. Size cutoff was I think 53 microns. (Archive) 

scratch.R -- playing around with stuff, who even knows what I was doing (Archive)

Explore_Quantiles_And_Size_Cutoff -- an early look into whether there were ever so few large particles that I could show there were more zeros than we would expect from chance. There weren't. (Archive)

Explore_Processed_UVP -- My earlist look at the UVP data (Archive)

