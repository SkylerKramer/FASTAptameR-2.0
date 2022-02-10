# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

# source the sub-tabs
source("./uiTabs/motifSearchTab.R")
source("./uiTabs/motifTrackerTab.R")

# merge the sub-tabs
motifTab <- tabPanel(
  "Motif",
  tabsetPanel(
    motifSearchTab,
    motifTrackerTab
  )
)