#' get_kegg_map
#' 
#' This function allows users to color individiul objects in KEGG pathway maps
#' and generate the output as an PNG image.
#' 
#' 
#' @param keggids KEGG IDs to be colored. e.g.: c("C05345","C05378")
#' @param pathwayid KEGG pathway ID: e.g.: "map01230"
#' @param bgcolor Vector of the size of the keggids argument to indicate
#' background color for each object. e.g.: c("red", "green")
#' @param fgcolor Vector of the size of the keggids argument to indicate border
#' and text color for each object. e.g.: c("black", "blue")
#' @return PNG image for the KEGG pathway with colored objects
#' @author Karan Uppal
#' @keywords ~pathway
get_kegg_map <- function(keggids, pathwayid, bgcolor, fgcolor) {
    
    # example:
    # g1<-get_kegg_map(keggids=c('C05345','C05378'),pathwayid='map01230',bgcolor=rep('red',2),fgcolor=rep('black',2))
    # map001100
    url <- KEGGREST::color.pathway.by.objects(bg.color.list = bgcolor, 
        fg.color.list = fgcolor, object.id.list = keggids, 
        pathway.id = pathwayid)
    png_image <- readPNG(getURLContent(url))
    fname <- paste(pathwayid, ".png", sep = "")
    writePNG(png_image, target = fname)
    
}
