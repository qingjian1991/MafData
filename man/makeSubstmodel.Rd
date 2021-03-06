% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeSubstmodel.R
\name{makeSubstmodel}
\alias{makeSubstmodel}
\title{makeSubstmodel}
\usage{
makeSubstmodel(submod = submod, sm = 1, colnm = NULL)
}
\arguments{
\item{submod}{load from the data: data("submod", package = "MafData")}

\item{sm}{\strong{sm}, Substitution model: 1) sm =1:5; sm = c("Type1_7_CpG","Type2_13_CpG","Type3_17","Type4_96","Type5_192")}
}
\value{
substmodel for sm
}
\description{
create the new sm for dndscv sm models.
}
\details{
This is for make new substitution models for dndscv.

you can create a new sm based on submod.
}
