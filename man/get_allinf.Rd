\name{get_allinf}
\alias{get_allinf}
\title{get_allinf}
\description{
This function allows users to map metabolite/chemical IDs
to additional fields such as KEGG BRITE IDs,
KEGG Pathway IDs, SMPDB IDs, HMDB Status, HMDB 
source, etc. depending upon the database.
}
\usage{
get_allinf(dataA,dbname)
}
\arguments{
  \item{dataA}{
data matrix with the first column corresponding to 
KEGG or HMDB IDs
}
  \item{dbname}{
Database name. e.g.: "KEGG", "HMDB", "T3DB"
}
}
\value{
Returns an object with chemical/metabolite IDs
merged with external IDs, BRITE categories,
pathway IDs, etc.
}
\author{
Karan Uppal
}

