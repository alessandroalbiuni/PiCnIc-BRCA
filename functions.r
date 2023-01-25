
#funzione per selezionare un set di eventi di un ogetto TRONCO secondo delle condizioni
select<-function(x, min.freq, forced.genes) {
  #unisce multipli eventi per gene in un unico evento
  x.sel <- as.alterations(x)
  #selezione di eventi in base alla frequenza e i geni
  x.sel <- events.selection(x.sel, 
                            filter.freq = min.freq, 
                            filter.in.names = forced.genes)
  
  ## Subset input
  x <- events.selection(x, 
                        filter.in.names = as.genes(x.sel))
  return(x)
}


#fun<ione per creare il modello di ricostruzione per BRCA
model<-function(BRCA,label, label.short) {
  #selezione da BRCA con min freq e geni mutex
  BRCA.select <- select(BRCA,
                        min_freq,
                        unique(
                            unlist(BRCA.mutex)
                        )
                        )
  BRCA.select <- annotate.description(BRCA.select,
                                      paste('BRCA', label, 'selection'))
  
  ## consolidare il dataset
  del <- consolidate.data(BRCA.select,
                          print = TRUE)
  #rimozione degli eventi equivoci
  if (length(del[["indistinguishable"]]) > 0) {
    for (i in seq_along(del[["indistinguishable"]])) {
      for (j in seq_len(nrow(del[["indistinguishable"]][[i]]))) {
        gene <- del[["indistinguishable"]][[i]][j,][2][[1]]
        type <- del[["indistinguishable"]][[i]][j,][1][[1]]
        BRCA.select <- delete.event(BRCA.select,
                                    gene = as.character(gene),
                                    type = as.character(type))
      }
    }
  }
  
  ## aggiunta delle ipotes
  BRCA.hypo <- BRCA.select
  
  ## ipotesi da mutex
  if (!is.null(BRCA.mutex)) {
    for (group in BRCA.mutex) {
      group <- group[group %in% as.genes(BRCA.hypo)]
      if (length(group) >= 2) {
        print(group)
        BRCA.hypo <- hypothesis.add.group(BRCA.hypo,
                                          FUN = OR,
                                          group = group,
                                          dim.min = length(group)
        )
      }
    }
  }
  
  ## aggiunta delle ipotesi
  BRCA.hypo <- hypothesis.add.homologous(BRCA.hypo)
  
  BRCA.hypo <- annotate.description(BRCA.hypo,
                                    as.description(BRCA.select))

    # uso di CAPRI
  BRCA.model <- tronco.capri(BRCA.hypo,
                             boot.seed = 42,
                             nboot = num_boot_iter)
  
  ## TRONCO.PLOT
  ## DAG del modello
  
  tronco.plot(
    BRCA.model,
    pathways = pathway.list,
    edge.cex = 1.5,
    legend.cex = .35,
    scale.nodes = .3,
    confidence = c('tp', 'pr', 'hg'),
    pathways.color = pathways.color,
    disconnected = F,
    height.logic = .3,
  )


  return(BRCA.model)
  
}
