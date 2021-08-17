#'@title plot_wgcna
#'@export
#'@examples
#'plot_wgcna(net)
#'@return ploting
#'@importFrom WGCNA labels2colors plotDendroAndColors
plot_wgcna <-
        function(net, ...){
                mergedColors = labels2colors(net$colors)
                plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                                           "Module colors",
                                           dendroLabels = FALSE, hang = 0.03,
                                           addGuide = TRUE, guideHang = 0.05,...)
        }
#' plot network
#' @importFrom tidygraph mutate as_tbl_graph centrality_degree
#' @import ggraph
#' @export
#' @param net, the result of function `wgcna`
#' @param corr, the result of function `corr.test`
#' @param centrality_degree_mode, the mode for centrality_degress function
#' @examples
#' library(data.table)
#' path <- system.file('data',package = 'mbOmic')
#' b <- readxl::read_xlsx(file.path(path,
#'                     '41598_2021_85433_MOESM2_ESM.xlsx'),
#'                     sheet = 'Table S2',skip = 1)
#' m <- readxl::read_xlsx(file.path(path,
#'                                  '41598_2021_85433_MOESM2_ESM.xlsx'),
#'                     sheet = 'Table S4',skip = 1)
#' setDT(b)
#' setDT(m)
#' cn <- c('Name_des',grep(x = names(m), 'BS|SS.*', value = TRUE))
#' m <- m[,..cn]
#' object <-
#'        mbSet(
#'              m = m,
#'              b = b
#'              )
#'
#' res <- cor.test(object, method = 'spearman')
#' net <- wgcna(object, minN = 2, power = 9,message = FALSE)
#' plot_network(net, res[abs(rho)>=0.85])
#' @author Congcong Gong
#' @usage
#' plot_network(
#'     net = net,
#'     corr = corr
#' )
#' @return ggraph
#' @importFrom ggplot2 aes coord_fixed scale_color_manual labs
plot_network <-
        function(net, corr, centrality_degree_mode = 'out',threshold=0.9){
                anno_c <- c(rep('orange',length(net$bs)),
                            as.character(interaction(net$colors, 'module', sep =' ')))
                names(anno_c) <- c(as.character(net$bs), names(net$colors))
                g <- as_tbl_graph(corr) %>%
                        tidygraph::mutate(`centrality degree`=centrality_degree(mode = centrality_degree_mode),
                               key=anno_c[name],
                               label = ifelse(`centrality degree` >= quantile(`centrality degree`,threshold), name, ''))

                aes_col <- gsub(' module', '',unique(anno_c))
                names(aes_col) <- unique(anno_c)
                ggraph(g, 'kk') +
                        geom_edge_link() +
                        geom_node_point(aes(colour = factor(key),
                                            size = `centrality degree`)) +
                        coord_fixed() +
                        labs(color = '') +
                        theme_graph(foreground = 'steelblue', fg_text_colour = 'white')+
                        geom_node_text(aes(label=label), size = 2)+
                        scale_color_manual(values = aes_col)
        }


