# AOS_2026_Workshop
Supporting materials for the AOS Workshop: Modeling bird movement with eBird and BirdFlowR

[Google doc page where you can type questions during the meeting.](https://docs.google.com/document/d/1ljsKB-lxMzf0NuCdcREeYlk_ajjlkixQvVzbOKEEmlE/edit?tab=t.0)

## Objective
To introduce users to the BirdFlowR R package and demonstrate how they can use it to 
visualize bird movement, generate synthetic routes, 
predict where birds inhabiting a specific location will be in the future, 
and quantify movement rates among arbitrary polygons such as states or counties.

## Preparing for Workshop
To prepare for the workshop or to work through these materials on your own 
please do the following:
1. Install R, RStudio, and BirdFlowR [using the "Standard Install" instructions](https://birdflow-science.github.io/BirdFlowR/articles/Installation.html#standard-install).
2. Download the workshop examples and data from this github page by clicking the green "<> Code" button near the top of the page and then selecting "Download Zip". If you don't see it make sure you are on the [main github page for the workshop](https://github.com/birdflow-science/AOS_2026_Workshop).
3. Run the following code in R or RStudio to install other packages used in the workshop but not required by BirdFlowR itself:
```r
cran_packages <- c(
  "remotes",
  "RColorBrewer",
  "patchwork",
  "cowplot",
  "rmarkdown",
  "knitr"
)

missing_cran <- cran_packages[!cran_packages %in% rownames(installed.packages())]

if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}
```

4. If you have institutional access to [Eduroam wifi](https://eduroam.org/)
   check that it is installed and working on your laptop. We also have guest wifi credentials that you should have received in an email and we can share at the workshop.

## Schedule
* 12:30–1:00 Optional installation help
* 1:00–1:10 Introduction to team and welcome
* 1:10–1:25 Introduction to BirdFlow  - Dan Sheldon
* 1:25–1:50 High-level walk through of R package with follow-along examples - Ethan Plunkett
* 1:50–2:15 Structured activity. Participants can pick a species of their choice and
interrogate the model.  - Ethan Plunkett
* 2:15–2:30 Break
* 2:30–3:00 Inferring population-specific phenology within the areas of interest - Yangkang Chen
* 3:00–3:30 Restructuring model movement data to correspond to specific regions (polygons) and time frames - Ethan Plunkett
* 3:30–4:00 Migratory connectivity comparing BirdFlow derived MC with tracking data derived MC  - Yuting Deng

## Files
* 



