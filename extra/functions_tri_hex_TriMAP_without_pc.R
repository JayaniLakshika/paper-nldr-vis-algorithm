Fit_TriMAP_data <- function(df_2_without_class, tem_dir){
  write.csv(df_2_without_class, file.path(tem_dir, "df_2_without_class.csv"), row.names = FALSE,
            quote = TRUE)
}


### Plot tSNE in 2D

plot_TriMAP_2D <- function(tSNE_df){
  tSNE_df_plot <- tSNE_df %>%
    ggplot(aes(x = TriMAP1,
               y = TriMAP2))+
    geom_point() +
    coord_equal()
  return(tSNE_df_plot)
}

## To find optimal bin value

find_opt_bin_val <- function (x)
{
  #h <- stats::IQR(x)
  #if (h == 0)
  #  h <- stats::mad(x, constant = 2)
  #if (h > 0)
  #  ceiling(diff(range(x))/(2 * h * length(x)^(-1/3)))
  #else 1L
  h <- 2 * IQR(x) / length(x)^(1/3) # bin width
  return(h)
}


## To create hexbin object
create_hex_bin <- function(tSNE_df, num_bins, shape){
  hb <- hexbin(tSNE_df$TriMAP1, tSNE_df$TriMAP2, num_bins, IDs = TRUE, shape = shape)
  return(hb)
}

### To Draw a hex bins with the points
draw_hex_bin <- function(tSNE_df, num_bins){
  hex_bin_plot <- ggplot(tSNE_df, aes(TriMAP1, TriMAP2)) +
    geom_hex(bins = num_bins) + geom_point()+
    scale_fill_viridis_c()

  return(hex_bin_plot)
}

## To create the dataframe which contain coordinate values with the hexbin ID
create_hexbin_ID_for_all <- function(df_2, tSNE_df, num_bins, shape){
  ### Merge tSNE dataset in 2D with original dataset which contains 4D coordinates
  df_2 <- df_2 %>%
    mutate(ID=row_number())

  df_new <- df_2 %>% inner_join(tSNE_df, by = "ID")

  ### Fit hexbins and store hexbin IDs
  hb <- create_hex_bin(tSNE_df, num_bins, shape)

  ### Add hexbin Ids as a column to the original dataset
  df_new$hb_id <- hb@cID

  ## To change column names to lower case
  names(df_new) <- tolower(names(df_new))

  ### Select specific columns to get the average number of 4D points in each bin
  ## Column names starts with x
  df <- df_new %>% select(starts_with("x"), hb_id)
  return(list(df = df, hb = hb, df_new = df_new))
}
## To create the dataframe which contain average 4D points
create_hexbin_avg_data <- function(df){
  df_b <- df %>%
    group_by(hb_id) %>%
    summarise_all(mean)

  ## To change column names to lower case
  names(df_b) <- tolower(names(df_b))

  ## Column names starts with x
  df_b <- df_b %>% select(starts_with("x"), hb_id)
  return(df_b)
}


## To extract the dataframe with contain hexbin cetroids data and mean hexbin coordinate data(using hcell2xy())
create_df_with_hex_bin_centroids <- function(df_b, hb){

  ### Fit hexbins and store hexbin IDs
  #hb <- create_hex_bin(tSNE_df, num_bins)

  ### Computes x and y coordinates from hexagon cell id's
  xy <- hcell2xy(hb)

  d_cell <- data.frame(x_val_center = xy$x, y_val_center = xy$y)

  ### Data of each cell (bin) which contain ID as hex_bin ID
  df_cell_data1 <- data.frame(ID = hb@cell, Cell_count = hb@count)

  df_cell_data <- cbind(df_cell_data1, d_cell)

  df_cell_data <- df_cell_data %>% rename( hb_id = "ID")

  ### Merge hexbin data with original mean dataset
  df_b_with_center_data <- df_b %>% inner_join(df_cell_data, by = "hb_id")
  df_b_with_center_data <- df_b_with_center_data %>%
    mutate(ID=row_number())
  return(df_b_with_center_data)
}

## To extract the dataframe with contain hexbin cetroids data and mean hexbin coordinate data(using hcell2xy())
#create_df_with_hex_bin_centroids_new <- function(df_b, tSNE_df){

### Fit hexbins and store hexbin IDs
#hb <- create_hex_bin(tSNE_df)

### Computes x and y coordinates from hexagon cell id's
#xy <- hcell2xy(hb)

#d_cell <- data.frame(x_val_center = xy$x, y_val_center = xy$y)

### Data of each cell (bin) which contain ID as hex_bin ID
#df_cell_data1 <- data.frame(ID = hb@cell, Cell_count = hb@count)

#df_cell_data <- cbind(df_cell_data1, d_cell)

#df_cell_data <- df_cell_data %>% rename( hb_id = "ID")

### Merge hexbin data with original mean dataset
#df_b_with_center_data <- df_b %>% inner_join(df_cell_data, by = "hb_id")
#df_b_with_center_data <- df_b_with_center_data %>%
#  mutate(ID=row_number())
#return(df_b_with_center_data)
#}

## To triangulate bin cetroids
trangulate_bin_centroids <- function(df_b_with_center_data){
  tr1 <- tri.mesh(df_b_with_center_data$x_val_center, df_b_with_center_data$y_val_center)
  return(tr1)
}

### Function to draw the trimesh in ggplot2
plot_trimesh_in_ggplot <- function(tr1) { ## Parameter should be the tri.mesh object
  tr_df <- data.frame(x = tr1$x, y = tr1$y) ## Create a dataframe with tri.mesh x and y coordinate values
  tr_df <- tr_df %>%
    mutate(ID=row_number()) ## To add ID numbers, beacuse to join with from and to points in tri$arcs

  trang <- triangles(tr1) ## To extract a triangulation data structure from an triangulation object
  trang <- as.data.frame(trang)

  tr_arcs_df1 <- data.frame(from = trang$node1, to = trang$node2) ## Create dataframe with from and to edges
  tr_arcs_df2 <- data.frame(from = trang$node1, to = trang$node3)
  tr_arcs_df3 <- data.frame(from = trang$node2, to = trang$node3)
  tr_arcs_df <- rbind(tr_arcs_df1, tr_arcs_df2, tr_arcs_df3) ## Create dataframe with from and to edges

  ## To obtain x and values of from to in a dataframe

  tr_from_to_df_coord <- data.frame(matrix(ncol = 4, nrow = 0)) ## Initialize an empty dataframe to store data in a     specific format
  colnames(tr_from_to_df_coord) <- c("x_from", "y_from", "x_to", "y_to") ## Define column names

  for (i in 1:dim(tr_arcs_df)[1]) {
    from <- tr_df[tr_arcs_df$from[i],1:2]
    to <- tr_df[tr_arcs_df$to[i],1:2]
    tr_from_to_df_coord[i,] <- c(x_from = from$x, y_from = from$y, x_to = to$x, y_to = to$y) ## Add vector as an       appending row to the dataframe
  }

  ## To draw the tri.mesh plot using ggplot
  tri_mesh_plot <- ggplot(tr_df, aes(x = x, y = y)) +
    geom_segment(aes(x = x_from, y = y_from, xend = x_to, yend = y_to), data = tr_from_to_df_coord) +
    geom_point() +
    coord_equal()

  return(tri_mesh_plot)
}

## To plot trimesh in ggplotly
plot_trimesh_in_ggplotly <- function(tr1) {
  tri_plotly <- plot_trimesh_in_ggplot(tr1) %>% ggplotly()
  return(tri_plotly)
}

### Distance function
#calculate_distance <- function(tSNE_df, from_to_data_tri) { ## Parameter1: Dataframe which contain ID, x, and y values of the data points

## Parameter 2: Dataframe which contain from and to edge points (output of trimesh)

#result <- rep(NA, dim(from_to_data_tri)[1]) ## Initialize an empty vector to store distance values
#distance_df <- data.frame(matrix(ncol = 3, nrow = 0)) ## Initialize an empty dataframe to store data in a specific format
#colnames(distance_df) <- c("from", "to", "distance") ## Define column names

#for (i in 1:dim(from_to_data_tri)[1]) {
#  result[i] <- sqrt(((tSNE_df[tSNE_df$ID  == from_to_data_tri$from[i],1] - tSNE_df[tSNE_df$ID == from_to_data_tri$to[i],1])^2) +
#                     ((tSNE_df[tSNE_df$ID == from_to_data_tri$from[i],2] - tSNE_df[tSNE_df$ID == from_to_data_tri$to[i],2])^2)) ## Calculate distance in 2D
#  distance_df[i,] <- c(from = from_to_data_tri$from[i], to = from_to_data_tri$to[i], distance = result[i]) ## Add vector as an appending row to the dataframe
#}
#return(distance_df) ## To return the data set with distances
#}


get_tri_data_frame <- function(tr1){
  tr_df <- data.frame(x = tr1$x, y = tr1$y) ## Create a dataframe with tri.mesh x and y coordinate values
  tr_df <- tr_df %>%
    mutate(ID=row_number()) ## To add ID numbers, beacuse to join with from and to points in tri$arcs

  trang <- triangles(tr1)
  trang <- as.data.frame(trang)

  tr_arcs_df1 <- data.frame(from = trang$node1, to = trang$node2) ## Create dataframe with from and to edges
  tr_arcs_df2 <- data.frame(from = trang$node1, to = trang$node3)
  tr_arcs_df3 <- data.frame(from = trang$node2, to = trang$node3)
  tr_arcs_df <- rbind(tr_arcs_df1, tr_arcs_df2, tr_arcs_df3) ## Create dataframe with from and to edges

  ## To obtain x and values of from to in a dataframe

  tr_from_to_df_coord <- data.frame(matrix(ncol = 6, nrow = 0)) ## Initialize an empty dataframe to store data in a     specific format
  colnames(tr_from_to_df_coord) <- c("from", "to", "x_from", "y_from", "x_to", "y_to") ## Define column names

  for (i in 1:dim(tr_arcs_df)[1]) {
    from <- tr_df[tr_arcs_df$from[i],1:3]
    to <- tr_df[tr_arcs_df$to[i],1:3]
    tr_from_to_df_coord[i,] <- c(from = from$ID, to = to$ID, x_from = from$x, y_from = from$y, x_to = to$x, y_to = to$y) ## Add vector as an       appending row to the dataframe
  }

  return(tr_from_to_df_coord)
}



## To get langevitour
get_langevitour <- function(df, df_b, df_b_with_center_data){

  tr1 <- trangulate_bin_centroids(df_b_with_center_data)
  tr_from_to_df <- get_tri_data_frame(tr1)

  ### Define type column
  df <- df %>% mutate(type = "data") ## original dataset

  df_b <- df_b %>% mutate(type = "model") ## Data with summarized mean

  df_exe <- bind_rows(df_b, df)

  langevitour(df_exe[1:(length(df_exe)-2)], lineFrom = tr_from_to_df$from , lineTo = tr_from_to_df$to, group = df_exe$type)
}

## Distance function
cal_dist <- function(tr_from_to_df_coord){
  tr_from_to_df_coord$distance <- lapply(seq(nrow(tr_from_to_df_coord)), function(x) {
    start <- unlist(tr_from_to_df_coord[x, c("x_from","y_from")])
    end <- unlist(tr_from_to_df_coord[x, c("x_to","y_to")])
    sqrt(sum((start - end)^2))})
  distance_df <- tr_from_to_df_coord %>% select("from", "to", "distance")
  distance_df$distance <- unlist(distance_df$distance)
  return(distance_df)
}

## To plot the distribution of distance
plot_dist <- function(distance_df){
  distance_df$group <- "1"
  dist_plot <- ggplot(distance_df, aes(x = group, y = distance)) +
    geom_quasirandom()+
    ylim(0, max(unlist(distance_df$distance))+ 0.5) + coord_flip()
  return(dist_plot)
}

## To plot the distribution of distance (plotly)
plot_dist2 <- function(distance_df){
  distance_df$group <- "1"
  dist_plot <- ggplot(distance_df, aes(x = group, y = distance)) +
    geom_quasirandom()+
    ylim(0, max(unlist(distance_df$distance))+ 0.5) + coord_flip()

  dist_plot_plty <- ggplotly(dist_plot)%>%
    layout(dragmode = "select", hovermode = "x unified") %>%
    highlight(on = "plotly_selected", dynamic = TRUE, color = toRGB("red"), off = "plotly_doubleclick", persistent = TRUE) %>%
    event_register(event = "plotly_brushed")
  return(dist_plot_plty)
}

get_langevitour_with_dist_criteria <- function(df, df_b, df_b_with_center_data, distance_df){

  ### Define type column
  df <- df %>% mutate(type = "data") ## original dataset

  df_b <- df_b %>% mutate(type = "model") ## Data with summarized mean

  df_exe <- bind_rows(df_b, df)

  tr1 <- trangulate_bin_centroids(df_b_with_center_data)
  tr_from_to_df <- get_tri_data_frame(tr1)

  sorted_distance_df <- distance_df[order(distance_df$distance),] ## Sort the distances

  # diff calculates the difference between consecutive pairs of
  #  vector elements
  dist_difference <- cbind(sorted_distance_df, rbind(NA, apply(sorted_distance_df, 2, diff)))[,c(1,2,3,6)] # apply diff to each column of sorted_distance_df, bind an NA row to the beginning,
  #  and bind the resulting columns to the original sorted_distance_df
  names(dist_difference) <- c('from', 'to', 'distance', 'difference')

  ## Set the maximum difference as the criteria
  distance_df_small_edges <- distance_df %>% filter(distance < ceiling(dist_difference[which(dist_difference$difference == max(dist_difference$difference, na.rm=TRUE)),3]))

  langevitour(df_exe[1:(length(df_exe)-2)], lineFrom = distance_df_small_edges$from, lineTo = distance_df_small_edges$to, group = df_exe$type)
}

get_langevitour_with_dist_criteria2 <- function(df, df_b, df_bin_centroids, benchmark_value, distance_df){

  ### Define type column
  df <- df %>% mutate(type = "data") ## original dataset

  df_b <- df_b %>% mutate(type = "model") ## Data with summarized mean

  df_exe <- bind_rows(df_b, df)

  #tr1 <- trangulate_bin_centroids(df_b_with_center_data)

  ## Set the maximum difference as the criteria
  distance_df_small_edges <- distance_df %>% filter(distance < benchmark_value)
  ## Since erase brushing is considerd.



  langevitour(df_exe[1:(length(df_exe)-2)], lineFrom = distance_df_small_edges$from, lineTo = distance_df_small_edges$to, group = df_exe$type)
}

## To color according to hover selection (long edges in red)

color_long_edges <- function(distance_df, benchmark_value, tr1){
  tr_df <- data.frame(x = tr1$x, y = tr1$y)

  tr_from_to_df_coord <- get_tri_data_frame(tr1)
  distance_df_small_edges <- distance_df %>% filter(distance < benchmark_value)
  distance_df_long_edges <- distance_df %>% filter(distance >= benchmark_value)

  distance_df_small_edges$type <- "small_edges"
  distance_df_long_edges$type <- "long_edges"

  distance_edges <- bind_rows(distance_df_small_edges, distance_df_long_edges)

  tr_from_to_df_coord_with_group <- merge(tr_from_to_df_coord, distance_edges,  by=c("from","to"))


  ## To draw the tri.mesh plot using ggplot
  tri_mesh_plot <- ggplot(tr_df, aes(x = x, y = y)) +
    geom_segment(aes(x = x_from, y = y_from, xend = x_to, yend = y_to, color = type), data = tr_from_to_df_coord_with_group) +
    geom_point() +
    coord_equal() +
    scale_colour_manual(values = c("#de2d26", "#636363"))
  return(tri_mesh_plot)

}
