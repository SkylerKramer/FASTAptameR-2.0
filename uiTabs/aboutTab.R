# Copyright (C) 2022 Skyler T. Kramer
# Full GNU GPL version 3 license found in LICENSE.txt

aboutTab <- tabPanel(
  "About",
  # logo
  div(img(src = "FASTAptameR2-0_logo.png"), style = "text-align: center;"),
  
  # title
  h1("About FASTAptameR 2.0"),
  
  # general description
  p(
    "FASTAptameR 2.0 is an R-based update of FASTAptamer. Like its predecessor, FA2 is an open-source
         toolkit designed to analyze populations of sequences resulting from combinatorial selections. This
         updated version features a user interface (UI), interactive graphics, more modules, and a faster
         implementation of the original clustering algorithm."
  ),
  
  # user guide
  p(HTML(paste0(
    "The user guide is available ",
    a(href = "https://github.com/SkylerKramer/FASTAptameR-2.0/blob/main/UserGuide.pdf", "here", .noWS = "outside", target = "_blank"),
    "."
  ))),
  
  # inputs and outputs
  p(HTML(paste0(
    "To get started, supply a FASTA/Q file to the Count module. Sample FASTQ files can be downloaded ",
    a(href = "https://burkelab.missouri.edu/fastaptamer.html", "here", .noWS = "outside", target = "_blank"),
    "."
  ))),
  
  # availability section
  h3("Availability"),
  
  # list of accession points
  tags$ul(
    tags$li(HTML(paste0(
      "Website (direct access): ",
      a(href = "https://fastaptamer2.missouri.edu/", "fastaptamer2.missouri.edu", target = "blank")
    ))),
    
    tags$li(HTML(paste0(
      "Docker (local download of platform): ",
      a(href = "https://hub.docker.com/repository/docker/skylerkramer/fastaptamer2", "skylerkramer/fastaptamer2", target = "blank")
    ))),
    
    tags$li(HTML(paste0(
      "Github (code download): ",
      a(href = "https://github.com/SkylerKramer/FASTAptameR-2.0/", "SkylerKramer/FASTAptameR-2.0", target = "blank")
    )))
    
    # tags$li(HTML(paste0(
    #   "User guide: ",
    #   a(href = "https://github.com/SkylerKramer/FASTAptameR-2.0/blob/main/UserGuide.pdf", "UserGuide.pdf", target = "_blank")
    # )))
  ),
  
  # citation
  h3("Citation"),
  p("If you use or modify FASTAptameR 2.0, please cite: https://fastaptamer2.missouri.edu"),
  
  # contact section
  h3("Contact"),
  p(
    "FA2 is maintained by Skyler Kramer and Donald Burke. To report bugs or request features, please use one of the following: "
  ),
  
  # list of contacts
  tags$ul(
    tags$li(HTML(paste0(
      "Github: ",
      a(href = "https://github.com/SkylerKramer/FASTAptameR-2.0/issues", "SkylerKramer/FASTAptameR-2.0", target = "blank")
    ))),
    
    tags$li(HTML(paste0(
      "Twitter: ",
      a(href = "https://twitter.com/BurkeLabRNA", "@BurkeLabRNA", target = "_blank")
    ))),
    
    tags$li("Email: burkelab@missouri.edu")
  ),
  
  # update section
  h3("The most recent version of FA2 is from Feb. 11, 2022.")
)