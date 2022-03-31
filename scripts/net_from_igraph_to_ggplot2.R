library(igraph)
library(ggplot2)
library(data.table)
library(plotly)


# net est une data.frame qui contient ce qui est à afficher, 
# donc filtrée en amont en fonction des lignées sélectionnées par exemple
# avec au moins 2 colonnes (fragment 1, fragment 2) et 
# des colonnes d'annotations (ex: nb de lignées, distance, type)

# annots sera la data.frame décrivant les noeuds 
# (fragment.ID (pour les other: chr_start, bait: SYMBOL), 
# type.fragments (other/bait) 
# SYMBOL (NA if not a bait) etc) 

# Connstruction du graph
G <- graph_from_data_frame(net)

# récupérer le layout (en fait les positions des noeuds sur le plot)
# et les arrêtes
L <- as.data.frame(layout.kamada.kawai(G)) # exemple, d'autre sont dispo


## annotation des noeuds
# match noms de noeuds avec frangment.ID
m <- match(names(V(G)),annots$fragment.ID)
# ordonne (et duplique si besoin) les fragment.ID pour 
# correspondre aux noms des noeuds
ids <- annots$fragment.ID[na.omit(m)]
# Attribue les annotations 
V(G)$shape <- annots$type.fragments[match(ids,annots$fragment.ID)]
V(G)$label <- annots$SYMBOL[match(ids,annots$fragment.ID)]
# si on rajoute des données d'expression par exemple
V(G)$expression <- annots$expression[match(ids,annots$fragment.ID)]


## annotation du layout, récupérations des annotation du graph
# récupère les noeuds annotés
vs <- V(G)
# récupère les arrêtes annotées
es <- as.data.frame(get.edgelist(G))
# compte le nombre de noeuds
Nv <- length(vs)
# compte le nombre d'arrêtes'
Ne <- nrow(es)

# annotation du layout
L$name <- vs$label
L$shape <- vs$shape
# si on a des données continues type expression, on crée des bin/catégories (factor)
L$mode <- cut(vs$expression,breaks=seq(-15,15,0.1))
# créer un vecteur de couleurs correspondant aus niveaux du factor
color.logfc <- colorRampPalette(c(colours()[c(131,124)],"white",colours()[c(34,36)]))(length(levels(L$mode)))
names(color.logfc) <- levels(L$mode)

# création d'une data.frame avec les arrêtes annotées
# Pour chaque arrête:
edge_shapes <- rbindlist(lapply(1:Ne,function(i,es,edges.g,L){
  # récupère l'ID du fragment1
  v0 <- as.vector(es[i,]$V1)
  # récupère l'ID du fragment2
  v1 <- as.vector(es[i,]$V2)
  # Création de la data.frame 
  type <- as.vector(edges.g$type[i])
  distance <- as.vector(edges.g$distance[i])
  nb.lignees <- as.vector(edges.g$nb.lignees[i]) # peut-être à discretiser comme pour expression
  edge_shape <- data.table(frag1=v0,frag2=v1, # nom des noeuds
  							V1=L$V1[L$name==v0],V2=L$V2[L$name==v0], # position du fragment1
  							xend=L$V1[L$name==v1],yend=L$V2[L$name==v1], # position du fragment2
  							type, nb.lignees, distance,
  							name=NA # add a column name to be able to use name in ggplot()
                           )
  return(edge_shape)
},es=es,edges.g=E(G),L))


# ggplot 
g <- ggplot(L,aes(x=V1,y=V2,label=name)) + # layout
    geom_segment(data=edge_shapes,aes(xend=xend,yend=yend, color=nb.lignees,size=distance, linetype=type)) + # edges, if no name column=>error
    geom_point(aes(shape=shape,fill=mode),size = 2, stroke = 4) + # plot noeuds (mode si on a l'info)
    scale_fill_manual(values=color.logfc) + 
    scale_linetype_manual(values=c(bo="solid",bb="dashed")) +
    scale_shape_manual(values=c(other=21,bait=24)) + theme_void() 

g <- ggplotly(g,tooltip=c("x","y","label"))













