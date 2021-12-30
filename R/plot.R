#' plot_coExpress
#'
#' Plotting the dendro of metabolites modules.
#'
#' @export
#' @return ploting
#' @importFrom WGCNA labels2colors plotDendroAndColors
#' @param net output of coExpress function
#' @examples
#' library(data.table)
#' path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')
#' load(path)
#' object <-
#'        mbSet(
#'              m = metabolites,
#'              b = genera
#'              )
#' res <- corr(object, method = 'spearman')
#' net <- coExpress(object, minN = 2, power = 9, message = FALSE)
#' plot_network(net, res[abs(rho)>=0.85])
#'
plot_coExpress <-
        function(net){
                mergedColors = labels2colors(net$colors)
                plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
        }
#' plot_network
#'
#' plotting the network of metabolites and OTU. Orange nodes represent the OTU, while the other color represent the metabolite. Same color metabolites nodes are constructed in the same modules.
#'
#' @importFrom tidygraph mutate as_tbl_graph centrality_degree
#' @import ggraph
#' @export
#' @param net result of function `wgcna`
#' @param corr result of function `corr.test`
#' @param centrality_degree_mode the mode for centrality_degress function
#' @param threshold, nodes whose centrality degree greater than the theshold will show labels
#' @param show_text logical whether add the nodes labels
#' @author Congcong Gong
#' @usage
#' plot_network(
#'     net,
#'     corr,
#'     centrality_degree_mode = 'out',
#'     threshold = 0.9,
#'     show_text = TRUE
#' )
#' @return ggraph
#' @importFrom ggplot2 aes coord_fixed scale_color_manual labs
#' @importFrom stats quantile
#' @importFrom tidygraph mutate
#' @examples
#' library(data.table)
#' path <- system.file('extdata', 'metabolites_and_genera.rda', package = 'mbOmic')
#' load(path)
#' mb <-
#'      mbSet(
#'            m = metabolites, b = genera
#'            )
#'            mb
#' ## Samples( 12 ):  BS1 BS2 BS3 BS4 BS5 BS6 SS1 SS2 SS3 SS4 SS5 SS6
#' ## number of OTU:  18
#' ## number of metabolites:  247
#' net <- coExpress(mb,message = TRUE,threshold.d = 0.02, threshold = 0.8, power = 9)
#' spearm <- corr(mb, method = 'spearman')
#' plot_network(net, spearm[abs(rho) >= 0.75 & p <= 0.05], threshold = 0.75, show_text = FALSE)
plot_network <-
        function(net, corr, centrality_degree_mode = 'out',threshold=0.9, show_text = TRUE){

                bN <- as.character(unique(corr$b))
                anno_c <- c(rep('orange',length(bN)),
                            as.character(interaction(net$colors, 'module', sep =' ')))
                names(anno_c) <- c(bN, names(net$colors))
                g <- as_tbl_graph(corr)
                g <- mutate(g, `centrality degree`=centrality_degree(mode = centrality_degree_mode),
                               key=anno_c[name],
                               label = ifelse(`centrality degree` >= quantile(`centrality degree`,threshold), name, ''))

                aes_col <- gsub(' module', '', unique(anno_c))
                names(aes_col) <- unique(anno_c)
                g <- ggraph(g, 'kk') +
                        geom_edge_link() +
                        geom_node_point(aes(colour = factor(key),
                                            size = `centrality degree`)) +
                        coord_fixed() +
                        labs(color = '') +
                        theme_graph(foreground = 'steelblue', fg_text_colour = 'white', base_family = 'Helvetica')+
                        scale_color_manual(values = aes_col)
                if (show_text) {
                        g <- g + geom_node_text(aes(label=label), size = 2)
                }
                g
        }


