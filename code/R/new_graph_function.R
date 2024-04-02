## transmission tree
library(epicontacts)
library(igraph)
library(tidyverse)
library(readxl)


new_fun <- function (x, x_axis, edge_alpha = NULL, network_shape = "branching", 
                     root_order = "subtree_size", node_order = "subtree_size", 
                     reverse_root_order = FALSE, reverse_node_order = FALSE, rank_contact = x_axis, 
                     reverse_rank_contact = FALSE, lineend = c("butt", "round", 
                                                               "square"), unlinked_pos = c("bottom", "top", "middle"), 
                     position_dodge = FALSE, parent_pos = c("middle", "top", "bottom"), 
                     custom_parent_pos = NULL, y_label = NULL, node_label = NULL, shapes=NULL,
                     y_coor = NULL, igraph_type = NULL, col_pal = cases_pal, edge_col_pal = edges_pal, 
                     ...) 
{
  def <- as.list(args(vis_epicontacts))
  def$node_size <- 5
  def$edge_width <- 1
  def$size_range <- c(3, 10)
  args <- list(...)
  node_shape <- get_val("node_shape", def, args)
  node_color <- get_val("node_color", def, args)
  node_size <- get_val("node_size", def, args)
  edge_color <- get_val("edge_color", def, args)
  edge_width <- get_val("edge_width", def, args)
  edge_linetype <- get_val("edge_linetype", def, args)
  NA_col <- get_val("NA_col", def, args)
  size_range <- get_val("size_range", def, args)
  width_range <- get_val("width_range", def, args)
  thin <- get_val("thin", def, args)
  legend_max <- get_val("legend_max", def, args)
  parent_pos <- match.arg(parent_pos)
  unlinked_pos <- match.arg(unlinked_pos)
  lineend <- match.arg(lineend)
  x <- x[i = !is.na(x$linelist$id), j = !is.na(x$contacts$from) & 
           !is.na(x$contacts$to)]
  x_axis <- assert_x_axis(x, x_axis)
  not_in_ll <- sum(!get_id(x, "contacts") %in% get_id(x, "linelist"))
  contacts_rm <- sum(is.na(get_pairwise(x, x_axis)))
  na_x_axis <- is.na(x$linelist[[x_axis]])
  msg <- "%s nodes and %s edges removed as x_axis data is unavailable"
  sm <- not_in_ll + sum(na_x_axis) + contacts_rm
  if (sm > 0) {
    warning(sprintf(msg, not_in_ll + sum(na_x_axis), contacts_rm))
  }
  x <- x[!na_x_axis]
  x <- thin(x, what = "contacts")
  if (thin) {
    x <- thin(x)
  }
  if (nrow(x$contacts) == 0L) {
    stop("No contacts found between cases with available x_axis data")
  }
  #node_shape <- assert_node_shape(x$linelist, node_shape, "node_shape", shapes = c(`Travel associated` = "plane", `Locally acquired` = "circle"))
  node_color <- assert_node_color(x$linelist, node_color, "node_color")
  node_size <- assert_node_size(x$linelist, node_size, "node_size")
  edge_color <- assert_edge_color(x$contacts, edge_color, "edge_color")
  edge_width <- assert_edge_width(x$contacts, edge_width, "edge_width")
  edge_linetype <- assert_edge_linetype(x$contacts, edge_linetype, 
                                        "edge_linetype")
  edge_alpha <- assert_edge_alpha(x, edge_alpha)
  node_order <- assert_node_order(x, node_order)
  root_order <- assert_root_order(x, root_order)
  rank_contact <- assert_rank_contact(x, rank_contact)
  custom_parent_pos <- assert_custom_parent_pos(custom_parent_pos)
  if (length(edge_width) > 1 | !inherits(edge_width, c("numeric", 
                                                       "integer"))) {
    msg <- paste("edge width must be a single number if method = 'ggplot' (cannot be", 
                 "mapped to a variable because scale_size is reserved for node_size)")
    stop(msg)
  }
  if (!is.null(y_label)) {
    if (!y_label %in% names(x$linelist)) {
      stop("y_label does not exist in linelist")
    }
    if (!position_dodge) {
      stop("position_dodge must be TRUE if y-axis labels are specifed")
    }
    if (!is.null(igraph_type)) {
      stop("igraph_type cannot be specified with y-axis labels")
    }
  }
  if ("R_i" %in% c(node_color, node_size, node_order, root_order)) {
    x$linelist$R_i <- vapply(x$linelist$id, function(i) sum(x$contacts$from == 
                                                              i, na.rm = TRUE), numeric(1))
  }
  nodes <- x$linelist
  edges <- x$contacts
  if (is.null(y_coor)) {
    coor <- get_coor(x, x_axis = x_axis, position_dodge = position_dodge, 
                     root_order = root_order, reverse_root_order = reverse_root_order, 
                     node_order = node_order, reverse_node_order = reverse_node_order, 
                     unlinked_pos = unlinked_pos, axis_type = "none", 
                     parent_pos = parent_pos, custom_parent_pos = custom_parent_pos, 
                     method = "ggplot", igraph_type = igraph_type)
    nodes$y <- coor$y
  }
  else {
    nodes$y <- y_coor
  }
  nodes$x <- x$linelist[[x_axis]]
  nodes$subtree_size <- coor$subtree_size
  if (network_shape == "rectangle") {
    df <- get_g_rect(nodes, edges)
    i_ind <- match(edges$to, nodes$id)
    inf_ind <- match(edges$from, nodes$id)
    df1 <- data.frame(y = nodes$y[i_ind], yend = nodes$y[i_ind], 
                      x = nodes$x[inf_ind], xend = nodes$x[i_ind])
    df1 <- cbind(df1, edges[!names(edges) %in% c("from", 
                                                 "to")])
    if (!is.null(df)) {
      df <- df[apply(df[, 1:4], 1, function(xx) !any(is.na(xx))), 
      ]
    }
    df <- rbind(df, df1)
  }
  else if (network_shape == "branching") {
    i_ind <- match(edges$to, nodes$id)
    inf_ind <- match(edges$from, nodes$id)
    df <- data.frame(y = nodes$y[inf_ind], yend = nodes$y[i_ind], 
                     x = nodes$x[inf_ind], xend = nodes$x[i_ind])
    df <- cbind(df, edges[!names(edges) %in% c("from", "to")])
    df <- df[apply(df[, 1:4], 1, function(xx) !any(is.na(xx))), 
    ]
    to_node <- rep(TRUE, nrow(df))
  }
  if (!is.null(node_color)) {
    if (inherits(nodes[[node_color]], c("factor", "character", 
                                        "logical"))) {
      cols <- fac2col(factor(nodes[, node_color]), col_pal, 
                      NA_col, TRUE)
      vals <- cols$leg_col
      names(vals) <- cols$leg_lab
      col_pal <- scale_fill_manual(values = vals, na.value = NA_col)
    }
    else if (inherits(nodes[[node_color]], "Date")) {
      dates <- pretty(nodes[[node_color]])
      numeric_node_color <- as.numeric(nodes[[node_color]])
      node_color <- paste0(node_color, "_")
      nodes[[node_color]] <- numeric_node_color
      if (missing(col_pal)) {
        col_pal <- scale_fill_continuous(breaks = as.numeric(dates), 
                                         labels = dates)
      }
      else {
        cols <- col_pal(10)
        col_pal <- scale_fill_gradientn(colors = cols, 
                                        breaks = as.numeric(dates), labels = dates)
      }
    }
    else if (inherits(nodes[[node_color]], c("numeric", "integer"))) {
      if (missing(col_pal)) {
        col_pal <- scale_fill_continuous()
      }
      else {
        cols <- col_pal(10)
        col_pal <- scale_fill_gradientn(colors = cols)
      }
    }
  }
  else {
    col_pal <- NULL
  }
  if (!is.null(edge_color)) {
    if (inherits(edges[[edge_color]], c("factor", "character"))) {
      cols <- fac2col(factor(edges[, edge_color]), edge_col_pal, 
                      NA_col, TRUE)
      vals <- cols$leg_col
      names(vals) <- cols$leg_lab
      edge_col_pal <- scale_color_manual(values = vals, 
                                         na.value = NA_col)
    }
    else if (inherits(edges[[edge_color]], "Date")) {
      dates <- pretty(edges[[edge_color]])
      numeric_edge_color <- as.numeric(edges[[edge_color]])
      edge_color <- paste0(edge_color, "_")
      edges[[edge_color]] <- numeric_edge_color
      if (missing(edge_col_pal)) {
        edge_col_pal <- scale_color_continuous(breaks = as.numeric(dates), 
                                               labels = dates, na.value = NA_col)
      }
      else {
        cols <- edge_col_pal(10)
        edge_col_pal <- scale_color_gradientn(colors = cols, 
                                              breaks = as.numeric(dates), labels = dates, 
                                              na.value = NA_col)
      }
    }
    else if (inherits(edges[[edge_color]], c("numeric", "integer"))) {
      if (missing(edge_col_pal)) {
        edge_col_pal <- scale_color_continuous(na.value = NA_col)
      }
      else {
        cols <- edge_col_pal(10)
        edge_col_pal <- scale_color_gradientn(colors = cols, 
                                              na.value = NA_col)
      }
    }
  }
  else {
    edge_col_pal <- NULL
  }
  if (inherits(node_size, c("numeric", "integer"))) {
    point <- geom_point(data = nodes, aes_string(x = "x", shape=node_shape,
                                                 y = "y", fill = node_color), size = node_size)#, shape = 21)
    size_pal <- NULL
  }
  else {
    if (inherits(nodes[[node_size]], "character")) {
      stop("node_size cannot be mapped to character variable")
    }
    else if (inherits(nodes[[node_size]], "Date")) {
      dates <- pretty(nodes[[node_size]])
      numeric_node_size <- as.numeric(nodes[[node_size]])
      node_size <- paste0(node_size, "_")
      nodes[[node_size]] <- numeric_node_size
      size_pal <- scale_size(range = c(size_range[1], size_range[2]), 
                             breaks = scales::pretty_breaks(as.numeric(dates)), 
                             labels = dates)
    }
    else if (inherits(nodes[[node_size]], "factor")) {
      warning("Mapping factor to size; converting factors to integers.")
      lev <- levels(nodes[[node_size]])
      ind <- as.integer(nodes[[node_size]])
      node_size <- paste0(node_size, "_")
      nodes[[node_size]] <- ind
      size_pal <- scale_size(range = c(size_range[1], size_range[2]), 
                             breaks = scales::pretty_breaks(sort(unique(ind))), 
                             labels = lev)
    }
    else {
      size_pal <- scale_size(range = c(size_range[1], size_range[2]), 
                             breaks = scales::pretty_breaks())
    }
    point <- geom_point(data = nodes, aes_string(x = "x", shape=node_shape,
                                                 y = "y", size = node_size, fill = node_color))#, shape = 21)
  }
  if (x$directed) {
    arrow <- arrow(length = unit(0.015, "npc"), type = "closed", 
                   ends = "last")
  }
  else {
    arrow <- NULL
  }
  if (inherits(edge_alpha, c("numeric", "integer"))) {
    segment <- geom_segment(aes_string(x = "x", xend = "xend", 
                                       y = "y", yend = "yend", color = edge_color, linetype = edge_linetype), 
                            arrow = arrow, alpha = edge_alpha, lineend = lineend, 
                            size = edge_width)
  }
  else {
    segment1 <- geom_segment(aes_string(x = "xmid", xend = "xend", 
                                        y = "ymid", yend = "yend", color = edge_color, linetype = edge_linetype, 
                                        alpha = edge_alpha), lineend = lineend, size = edge_width)
    segment2 <- geom_segment(aes_string(x = "x", xend = "xend-2", 
                                        y = "y", yend = "yend", color = edge_color, 
                                        alpha = edge_alpha), lineend = lineend, arrow = arrow, linetype = "longdash", 
                             size = edge_width)
  }
  if (!is.null(node_label)) {
    nodes$y <- nodes$y
    lab <- geom_text(data = nodes, aes_string("x", "y", label = node_label), 
                     color = "black", size = 3)
  }
  else {
    lab <- NULL
  }
  if (!is.null(y_label)) {
    y_scale <- scale_y_continuous(name = NULL, breaks = sort(coor$y), 
                                  minor_breaks = NULL, labels = x$linelist[[y_label]][order(coor$y)], 
                                  expand = c(0.01, 0.01))
    gg_theme <- theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(), 
                      panel.grid.minor.y = element_blank())
  }
  else {
    y_scale <- NULL
    gg_theme <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
                      axis.title.y = element_blank(), panel.grid.major.y = element_blank(), 
                      panel.grid.minor.y = element_blank())
  }
  if (inherits(df$x, "Date")) {
    df$xmid <- as.Date((as.numeric(df$x) + as.numeric(df$xend))/2, 
                       origin = "1970-01-01")
  }
  else {
    df$xmid <- (df$x + df$xend)/2
  }
  df$ymid <- (df$y + df$yend)/2
  out <- ggplot(df)+  #segment1 + 
    segment2  + point +lab + col_pal + 
    edge_col_pal + size_pal + y_scale + theme_minimal() + 
    gg_theme + labs(x = x_axis)
  return(out)
}

environment(new_fun) <- asNamespace('epicontacts')

## current tree


make_tree_graph <- function(tree_out){
  
  cases <- tree_out%>%
    mutate(miami=ifelse(County=="Miami-Dade","Miami-Dade","Other")) %>%
    mutate(travel.status=factor(travel.status)) %>%
    mutate(infector=ifelse(travel.status=="Travel associated",
                           0,infector
    ))%>%
    filter(infector!=0 | ID %in% infector) %>%
    mutate(travel.status=str_replace(travel.status," ","\n"))
  
  linelist <- read_excel("data/Dengue 3 positives for phylogeographic analysis.xlsx") %>%
    dplyr::select(GenBank.acc.=`GenBank acc#`,zipcode=`ZIP Code1 (Miami only)`) %>%
    filter(GenBank.acc. %in% cases$GenBank.acc.) %>%
    #mutate(zipcode=factor(zipcode,labels=c(letters[1:22],LETTERS)))
    mutate(zipcode=factor(as.numeric(factor(zipcode))))
  
  ## Compute Re for each cluster
  epic <- make_epicontacts(
    linelist = cases %>%
      left_join(linelist)%>%
      mutate(Onset.Date=as.Date(Onset.Date,format="%Y-%m-%d")),
    contacts = cases %>%
      left_join(linelist)%>%
      filter(infector!=0) %>%
      dplyr::select(ID,infector,prob_source,travel.status,zipcode),
    id = "ID",
    from = "infector",
    to = "ID",
    directed = TRUE
  )
  
  ### fig attempt 1
  
  p <- new_fun(
    epic,
    edge_alpha="prob_source",
    x_axis = "Onset.Date",
    node_color = "travel.status",
    node_shape="travel.status",
    arrow_size = 0.1,
    node_size = 5,
    label = FALSE,
    height = 700,
    width = 700,
    #NA_col="black",
    #node_label="zipcode"
  )
  
  
  p + cowplot::theme_cowplot()+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(angle = 45,hjust=1)
    )+
    scale_x_date("Date of symptom onset",date_breaks = "1 months",date_labels = "%b %Y",
                 limits = as.Date(c("2022-05-01","2023-01-01"))
    )+
    scale_shape_manual("Case\nType",values=c(21,22)
                       )+
    scale_alpha_continuous("Compatibility\nScore",limits=c(0,1))+
    guides(fill="none")
  
}

p1 <- make_tree_graph(read.csv(file="output/20230929_istyruns/cases_assigned_PrunePi01_2023-09-29_maxPall.csv"))+
  ggtitle("1% reporting probability")

p5 <- make_tree_graph(read.csv(file="output/20230929_istyruns/cases_assigned_PrunePi05_2023-09-29_maxPall.csv"))+
  ggtitle("5% reporting probability")

p10 <- make_tree_graph(read.csv(file="output/20230929_istyruns/cases_assigned_PrunePi10_2023-09-29_maxPall.csv"))+
  ggtitle("10% reporting probability")

p15 <- make_tree_graph(read.csv(file="output/20230929_istyruns/cases_assigned_PrunePi20_2023-09-29_maxPall.csv"))+
  ggtitle("15% reporting probability")



cowplot::plot_grid(p1+theme(legend.position = "none"),
                   p5+theme(legend.position = "none"),
                   p10+theme(legend.position = "none"),
                   p15+theme(legend.position = "none"),
                   cowplot::get_legend(p10),
                   rel_widths = c(2,2,2,2,1),
                   nrow=1
                   )




