get_kegg_map <-
function(keggids,pathwayid,bgcolor,fgcolor){
    
    #example: g1<-get_kegg_map(keggids=c("C05345","C05378"),pathwayid="map01230",bgcolor=rep("red",2),fgcolor=rep("black",2))
    #map001100
    url<-KEGGREST::color.pathway.by.objects(bg.color.list=bgcolor,fg.color.list=fgcolor,object.id.list=keggids,pathway.id=pathwayid)
    png_image<-readPNG(getURLContent(url))
    fname<-paste(pathwayid,".png",sep="")
    writePNG(png_image,target=fname)
    
}
