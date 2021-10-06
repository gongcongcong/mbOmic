#'@title plot_coExpress
#'@export
#'@examples
#'\dontrun{
#'plot_coExpress(net)
#'}
#'@return ploting
#'@importFrom WGCNA labels2colors plotDendroAndColors
plot_coExpress <-
        function(net, ...){
                mergedColors = labels2colors(net$colors)
                plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                                           "Module colors",
                                           dendroLabels = FALSE, hang = 0.03,
                                           addGuide = TRUE, guideHang = 0.05,...)
        }
#' plot network
#' @title plot_network
#' @importFrom tidygraph mutate as_tbl_graph centrality_degree
#' @import ggraph
#' @export
#' @param net, the result of function `wgcna`
#' @param corr, the result of function `corr.test`
#' @param centrality_degree_mode, the mode for centrality_degress function
#' @param threshold, the nodes whose centrality degree greater than the theshold will show labels
#' @examples
#' \dontrun{
#' library(data.table)
#' path <- system.file('data',package = 'mbOmic')
#' load(file.path(path,'metabolites_and_genera.rda'))
#' object <-
#'        mbSet(
#'              m = metabolites,
#'              b = genera
#'              )
#' res <- corr(object, method = 'spearman')
#' net <- wgcna(object, minN = 2, power = 9,message = FALSE)
#' plot_network(net, res[abs(rho)>=0.85])
#' }
#' @author Congcong Gong
#' @usage
#' plot_network(
#'     net = net,
#'     corr = corr，
#'     centrality_degree_mode = 'out',
#'     threshold = 0.9
#' )
#' @return ggraph
#' @importFrom ggplot2 aes coord_fixed scale_color_manual labs
plot_network <-
        function(net, corr, centrality_degree_mode = 'out',threshold=0.9, show_text = TRUE){
                bN <- as.character(unique(corr$b))
                anno_c <- c(rep('orange',length(bN)),
                            as.character(interaction(net$colors, 'module', sep =' ')))
                names(anno_c) <- c(bN, names(net$colors))
                g <- as_tbl_graph(corr) %>%
                        tidygraph::mutate(`centrality degree`=centrality_degree(mode = centrality_degree_mode),
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


