
Exploring the impact of trait number and type on functional diversity
metrics in real-world ecosystems

Analyses by Kaitlin Kimmel and Tim Ohlert
9/3/2023

<p>
This code corresponds to the analyses proposed in the stage 1 Registered
Report “Exploring the impact of trait number and type on functional
diversity metrics in real-world ecosystem” in PLoS One. Please contact
<kaitlinakimmel@gmail.com> with questions regarding the code and data.
Please contact <Timothy.Ohlert@colostate.edu> for all other questions
regarding the manuscript.
</p>
<p>
This code uses the here::here() function to call data and save files.
This ensures that working directories do not have to be changed for each
user contributing to the code. Create folders: ‘data’, ‘R’, and
‘Figures’ so code will run easily.
</p>

## <u>Code </u>

Code is separated currently by site and/or distance matrix used and/or
predictor variable. Run code in the order provided to obtain outputs
necessary in subsequent steps. Note that for 2 & 3, scripts without
’\_euc’ are gower distance and scripts with ’\_euc’ are euclidean
distance.

1.  <b> Data Cleaning Codes </b>
    1.  CDRCommunityDataCleaning.R
    2.  KNZCommunityDataCleaning.R
    3.  SEVCommunityDataCleaning.R
    4.  CDRTraitCleaning.R
    5.  KNZTraitCleaning.R
    6.  SEVTraitCleaning.R
2.  <b> Functional Diversity Metric Calculations </b> These codes are
    processing intensive and take a long time to run (mainly because of
    the KDE functions). All the output is saves in the data folder to
    avoid running these each time the code is used.
    1.  CDRFDmetrics.R
    2.  CDRFDmetrics_Euc.R
    3.  KNZFDmetrics.R –> Konza did not have any communities that met
        the trait coverage threshold.
    4.  SEVFDmetrics.R
    5.  SEVDFmetrics_euc.R
3.  <b> AIC model selection </b>A bulk of this code is finding out best
    fits for correlations between number of traits, mean correlation of
    traits, maximum correlation of trait, minimum correlation of traits,
    and the calculated FD values. There is currently a separate script
    for each correlate and distance matrix. I do not have a master
    script that runs all of these simultaneously yet because we are
    still in the QA/QC process of making sure we have the correct models
    chosen for each FD metric, community, and distance matrix. Upon
    publication, a script will be written to run through this all
    without running each script individually.
    1.  n_traitStats.R
    2.  n_traitStats_euc.R
    3.  max_corrStats.R
    4.  max_corrStats_euc.R
    5.  mean_corrStats.R
    6.  mean_corrStats_euc.R
    7.  min_corrStats.R
    8.  min_corrStats_euc.R
4.  <b> Graphs </b> Originally, graphs were separated out by correlate
    and distance matrix. These scripts are now not used and only four
    scripts are necessary to create the graphs that will be presented in
    the publication.
    1.  graphs.R - number of traits; this code should be updated to pull
        model fits directly from n_traitStats.R and n_traitStats_euc.R
    2.  max_cor_graphs.R
    3.  min_cor_graphs.R
    4.  mean_cor_graphs.R

## <u> Data </u>

1.  All Cedar Creek data is pulled directly from pasta.lternet.edu
2.  Konza trait data (located in data/Konza)
    1.  firegraze_trait.csv - Trait data
    2.  SANA_grass_spp.csv - species names
    3.  sp_list.xlsx - species names
3.  Konza community data is pulled directly from pasta.lternet.edu
4.  Sevilleta trait data is pulled directly from pasta.lternet.edu
    except for height
    1.  located in data/Sev - sev_npp_height_biomass.csv - should be
        available on pasta.lternet.edu
5.  Sevilleta community data (located in data/Sev)
    1.  sevcommunity_KOtrait.csv - blue and black gramma sites
=======
# FDiv
>>>>>>> 9047cef5d254f6d3ea8f1dd7e02136318b3178e3
=======
# FDiv
>>>>>>> 9047cef5d254f6d3ea8f1dd7e02136318b3178e3
