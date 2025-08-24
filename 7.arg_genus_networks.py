library(readxl)
library(dplyr)
library(stringr)
library(igraph)
library(ggraph)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

if (!file.exists("arg_data.xlsx") || !file.exists("metadata.xlsx")) {
  message("Data files not found")
  quit(save = "no", status = 1)
}



arg_data <- read_excel("arg_data.xlsx")
metadata <- read_excel("metadata.xlsx")

full_data <- arg_data %>%
  left_join(metadata, by = "sample_id")

#important args
important_args <- c("ccrA", "mupA", "mprF")
arg_map <- c(
  "CCRA" = "ccrA",
  "STAPHYLOCOCCUS MUPA CONFERRING RESISTANCE TO MUPIROCIN" = "mupA",
  "BRUCELLA SUIS MPRF" = "mprF"
)

full_data$ARG_simple <- arg_map[toupper(full_data$ARG)]
full_data$ARG_simple[is.na(full_data$ARG_simple)] <- full_data$ARG[is.na(full_data$ARG_simple)]

make_network_edges <- function(data, group_name) {
  data %>%
    filter(ARG_simple %in% important_args, group == group_name) %>%
    filter(!is.na(Taxa)) %>%
    mutate(
      Genus_level = sapply(strsplit(Taxa, ";"), function(x) ifelse(length(x) >= 6, x[6], NA))
    ) %>%
    filter(!is.na(Genus_level), Genus_level != "Unassigned") %>%
    count(ARG_simple, Genus_level, name = "weight") %>%
    rename(from = ARG_simple, to = Genus_level)
}

edges_influent <- make_network_edges(full_data, "influent")
edges_effluent <- make_network_edges(full_data, "final effluent")

#Function to plot network 
plot_network <- function(edges_df) {
  if (nrow(edges_df) == 0) stop("No edges to plot!")
  
  all_nodes <- data.frame(name = unique(c(edges_df$from, edges_df$to)), stringsAsFactors = FALSE)
  net_graph <- graph_from_data_frame(d = edges_df, vertices = all_nodes, directed = FALSE)
  
  V(net_graph)$NodeType <- ifelse(V(net_graph)$name %in% important_args, V(net_graph)$name, "Genus")
  V(net_graph)$degree <- degree(net_graph)
  
  node_colors <- c(
    "ccrA" = brewer.pal(8, "Set1")[1],
    "mupA" = brewer.pal(8, "Set1")[2],
    "mprF" = brewer.pal(8, "Set1")[3],
    "Genus" = "gray70"
  )
  
  ggraph(net_graph, layout = "fr") +
    geom_edge_link(aes(width = weight, alpha = weight), color = "gray50") +
    scale_edge_width(range = c(0.5, 3), guide = "none") +
    scale_edge_alpha(range = c(0.4, 0.8), guide = "none") +
    geom_node_point(aes(fill = NodeType, size = degree), shape = 21, color = "black", stroke = 0.3) +
    geom_node_text(aes(label = name), repel = TRUE, size = 4, family = "Arial") +
    scale_fill_manual(values = node_colors, name = "Node Type") +
    scale_size(range = c(4, 10), name = "Node Degree") +
    theme_minimal(base_family = "Arial") +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(15, 15, 15, 15),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
    )
}
plot_influent <- plot_network(edges_influent)
plot_effluent <- plot_network(edges_effluent)

#vertical
combined_net_plot <- (plot_influent / plot_effluent) +
  plot_annotation(tag_levels = 'a', tag_prefix = "", tag_suffix = ".", tag_sep = "\n") &
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0, 1.02))


ggsave(
  filename = "ARG_genus_network.png",
  plot = combined_net_plot,
  width = 10,
  height = 16,
  units = "in",
  dpi = 300
)

print(combined_net_plot)
