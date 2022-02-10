# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

# source the sub-tabs
source("./uiTabs/seqEnrichTab.R")
source("./uiTabs/posEnrichTab.R")

# merge the sub-tabs
enrichTab <- tabPanel(
  "Enrichment",
  tabsetPanel(
    seqEnrichTab,
    posEnrichTab
  )
)