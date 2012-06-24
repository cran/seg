\name{nzcity}
\docType{data}
\alias{nzcity}
\title{Ethnic Groups in Sixteen New Zealand Cities}
\description{
This data set contains the digital boundaries of sixteen New Zealand cities and the population by ethnic groups (2006).
}
\usage{
data(nzcity)
}
\format{
A class of \code{list}. Each element contains the digital boundaries and the population data for one of the following New Zealand cities:

\tabular{rll}{
  [[1]] \tab AKC \tab Auckland \cr
  [[2]] \tab CHC \tab Christchurch \cr
  [[3]] \tab DUN \tab Dunedin \cr
  [[4]] \tab HMT \tab Hamilton \cr
  [[5]] \tab INV \tab Invercargill \cr
  [[6]] \tab LHT \tab Lower Hutt \cr
  [[7]] \tab MKC \tab Manukau \cr
  [[8]] \tab NAP \tab Napier \cr
  [[9]] \tab NEL \tab Nelson \cr
  [[10]] \tab NSC \tab North Shore \cr
  [[11]] \tab PAL \tab Palmerston North \cr
  [[12]] \tab POR \tab Porirua \cr
  [[13]] \tab TAU \tab Tauranga \cr
  [[14]] \tab UHT \tab Upper Hutt \cr
  [[15]] \tab WAI \tab Waitakere \cr
  [[16]] \tab WEL \tab Wellington
}

Each element is a class of \code{SpatialPolygonsDataFrame}, and the data slot contains population counts by census tracts (rows) and ethnic groups (columns). Each column indicates the following information:

\sQuote{AU} Area Unit names.

\sQuote{EUROPEAN} The number of Europeans who identified with one or more European ethnicities in the 2006 census. The largest ethnicity among European people is New Zealand European, followed by English, Dutch, Australian, South African, Scottish and Irish.

\sQuote{MAORI} The number of Maori people who identified with the Maori ethnic group in the 2006 census. 

\sQuote{PACIFIC} The number of Pacific peoples who identified with one or more Pacific ethnicities in the 2006 census. The seven largest ethnicities among Pacific peoples are Cook Island Maori, Fijian, Niuean, Samoan, Tokelauan, Tongan and Tuvaluan peoples.

\sQuote{ASIAN} The number of Asian people who identified with one or more Asian ethnicities in the 2006 census. The largest ethnicity among Asian people is Chinese, followed by Indian, Korean, Filipino, Japanese, Sri Lankan, Cambodian and Thai. 

\sQuote{MELAA} The number of people who were represented in MELAA in the 2006 census. The MELAA ethnic group includes people of Middle Eastern, African and Latin American ethnicities.

\sQuote{OTHER} The number of people who belonged to other groups that do not fit into the five other groupings.
}
\source{
\href{http://www.stats.govt.nz/Census/about-2006-census.aspx}{Statistics New Zealand} and licensed by Statistics New Zealand for re-use under the \href{http://creativecommons.org/licenses/by/3.0/nz/deed.en}{Creative Commons Attribution 3.0 New Zealand licence}.
}
\examples{
data(nzcity)
summary(nzcity$WEL[,2:7])
plot(nzcity$WEL, border = "Grey 50", col = "Grey 90")
plot(nzcity$WEL[nzcity$WEL$ASIAN > 500,], border = "Grey 50", col = "Grey 70", add = TRUE)
title(main = "Wellington City", xlab = "Areas with more than 500 Asians")
}
\keyword{datasets}
