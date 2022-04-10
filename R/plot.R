#' plot_coExpress
#'
#' Plotting the dendro of metabolites modules.
#'

#' @return ploting
#' @importFrom WGCNA labels2colors plotDendroAndColors
#' @param net output of coExpress function
#' @examples
#' library(data.table)
#' path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')
#' load(path)
#' names(metabolites)[1] <- 'rn'
#' m <- mSet(m = metabolites)
#' net <- coExpress(m, minN = 2, power = 9, message = FALSE)
#'
plot_coExpress <-
        function(net){
                mergedColors <- labels2colors(net$colors)
                plotDendroAndColors(net$dendrograms[[1]],
                                    mergedColors[net$blockGenes[[1]]],
                                    "Module colors",
                                    dendroLabels = FALSE,
                                    hang = 0.03,
                                    addGuide = TRUE,
                                    guideHang = 0.05)
        }
#' plot_network
#'
#' plotting the network of metabolites and OTU. Orange nodes represent the OTU, while the other color represent the metabolite. Same color metabolites nodes are constructed in the same modules.
#'
#' @export
#' @param net result of function `coExpress`
#' @param corr result of function `corr.test`
#' @param return whether return the igraph
#' @param interaction plot method
#' @param seed set seed for layout in interaction ploting
#' @author Congcong Gong
#' @return  igraph or graph
#' @importFrom igraph V graph_from_data_frame plot.igraph degree "V<-" layout_with_kk
#' @importFrom visNetwork toVisNetworkData visNetwork visGroups visLegend visIgraphLayout
#' @examples
#' library(data.table)
#' path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')
#' load(path)
#' names(genera)[1] <- 'rn'
#' names(metabolites)[1] <- 'rn'
#' b <- bSet(b = genera)
#' m <- mSet(m = metabolites)
#' res <- corr(m, b, method = 'spearman')
#' net <- coExpress(m, minN = 2, power = 9, message = FALSE)
#' plot_network(net, res[abs(rho)>=0.85])
plot_network <-
        function(net, corr, seed = 123, return = FALSE, interaction = TRUE){
                b_n <- unique(corr$b)
                m_n <- unique(corr$m)
                nodes <- c(b_n, m_n)
                nodesCol <-
                        c(rep('orange', length(b_n)), net$colors)
                names(nodesCol)[seq_along(b_n)] <- b_n
                g <- graph_from_data_frame(corr[, c("m", "b"), with = FALSE], directed = FALSE)
                V(g)$shape <- ifelse(V(g)$name %in% corr$b,
                                     "sphere",
                                     "circle")
                V(g)$size <- degree(g, mode = 'all')
                V(g)$color <- nodesCol[V(g)$name]
                if (interaction) {
                        visData <- toVisNetworkData(g)
                        visData$nodes$group <-
                                ifelse(visData$nodes$id %in% corr$b,
                                       "OTU",
                                       "metabolites")
                        visData$nodes$shape <-
                                ifelse(visData$nodes$id %in% corr$b,
                                       "star",
                                       "dot")
                        visData$nodes$value <- visData$nodes$size*100
                        nodes <- data.frame(
                                label = c("OTU", paste0("Module: ", unique(net$colors))),
                                shape = c("star", rep("dot", length(unique(net$colors)))),
                                color = c("orange",unique(net$colors))
                        )
                        out <- visNetwork(visData$nodes, visData$edges) |>
                                visIgraphLayout(smooth = TRUE, type = 'full', layout = 'layout_with_kk') |>
                                visLegend(useGroups = FALSE,
                                          addNodes = nodes)
                        print(out)
                        if (return) {
                                out
                        }
                        } else {
                                plot.igraph(g, vertex.label="", layout = layout_with_kk)
                                if (return) {
                                        g
                                }
                }

        }


