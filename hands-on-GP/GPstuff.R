install.packages("RcppOctave") # "package ‘RcppOctave’ is not available for this version of R"
install_version("RcppOctave", "0.18.1") ### the latest version available in Cran

av <- available.packages(filters=list())
av[av[, "Package"] == 'RcppOctave', ]
View(av)
"RcppOctave" %in% rownames(av)

library(devtools)
# latest stable version
install_github('RcppOctave', 'renozao')
install_github("https://github.com/renozao/RcppOctave") #same error than with row 2

# development version
install_github('RcppOctave', 'renozao', 'develop')

