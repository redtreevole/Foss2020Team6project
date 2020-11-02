# ---
# title: "Visually_determine_k-mer_optimum"
# author: "JaeJin Choi"
# date: "9/23/2020"
# ---

require(dplyr)
#require(tidyverse)
#require(readr) 

require(ape)

## below functions are for plotting
#require(ggplot2) #http://docs.ggplot2.org/
#require(reshape)
#require(grid) #necessary for ggplot function, ex 'unit'
require(farver)

#' @name visual_optimum
#' @Description Input distance matrices, calculate, organize, and output a dataframe for plotting purpose
#' @param load_path path_to_distance_matrices (string)
#' @param name_regex regular expression to specify files (string; NULL)
#' @return dataframe for plotting FFP optimum
#' @examples  ffP_visual_optimum("./distance_matrices", "matrix"). Pick up only the files have 'matrix' string in names

#ffp_visual_optimum <- function(load_path="./", name_regex=NULL) #current path, pattern=NULL
ffp_visual_optimum <- function(file_list=list(), name_regex=NULL) #current path, pattern=NULL
{
  read_to_dist_df <- function(load_path)
  {
    tmp_colnames <-  c(seq(1:(length(readLines(load_path))))) #temporarly assign to fit a number of columns when filled in read.table
    
    input_table <- read.table(file = load_path
                              , sep = "" # all white spaces
                              , header= F
                              , skip = 1 #skip the first row (a number of OTUS)
                              , as.is = T #recycling rule to guess variable class
                              , col.names = tmp_colnames
                              , fill = T #able to read triangular matrix
                              
    )
    
    input_df <- as.data.frame(input_table[,-1]) #exclude the first column
    
    colnames(input_df) <- input_table$X1
    row.names(input_df) <- input_table$X1
    
    dist_df <- as.dist(input_df) #convert to low triangular distance matrix that is feasible to phylo
    
    return (dist_df)
    
    
  }
  
  get_tree_dist <- function(mphylo_trees)
  {
    m_tree_dist <- as.numeric(c())
    #compare topology of adjacent trees (according to index)
    for (cy1 in seq(1:(length(mphylo_trees)-1))) #dist.topo(cy1, cy1+1)
    {
      #print(dist.topo(mphylo_trees[[cy1]], mphylo_trees[[cy1+1]])[1])
      m_tree_dist <- append(m_tree_dist, dist.topo(mphylo_trees[[cy1]], mphylo_trees[[cy1+1]])[1]) #just one element return
      
    }
    
    m_tree_dist <- append(m_tree_dist, NA) #add last element as NA to fit size for m_combined_df
    
    #print(m_tree_dist)
    return (m_tree_dist)
  }
  
  
  
  #read distance matrix file, either symmetric or triangular matrix
  # file_list <- list.files(path=load_path
  #                         , full.names=T #full file path, F will get only file name
  #                         , pattern = name_regex
  # ) 
  
  #print(paste0("A number of input files: ", length(file_list)))
  
  
  #for sorting purpose; avoid lexigraphical order
  m_feature_length <- unlist(lapply(file_list, 
                                    function(fn) as.numeric(tail(strsplit(fn, split="\\.")[[1]], n=1))
  ))
  
  file_list <- file_list[sort.list(m_feature_length)] #sort file_list by ascending feature_length order
  m_feature_length <- sort(m_feature_length) #also sort
  
  
  ###the rest lapply process in the element order of file_list
  m_dist_mat <- lapply(file_list, function(fn) read_to_dist_df(fn)) #can only accept symmetric matrix
  
  
  ### calculate standard deviation of distance matrix
  m_jsd_sd <- lapply(m_dist_mat, function(fn) sd(fn))
  
  
  ##obtain total count of values >= 0.99999
  # m_jsd_upperlimit_count <- lapply(m_dist_mat, function(fn) get_jsd_upperlimit_count(fn)) #using custom function
  m_jsd_upperlimit_count <- lapply(m_dist_mat, function(fn) sum(fn>=0.99999)) #0.99999 is almost equivalent to 1.0
  
  
  # Get BIONJ, a Neighbor-Joining variant, tree from given dist.matrix 
  m_phylo_tree <- lapply(m_dist_mat, function(fn) bionj(as.dist(fn)))
  
  
  # Calculate Robinson-Foulds (RF) distance (topological distance between two trees)
  m_tree_dist <- get_tree_dist(m_phylo_tree)
  
  
  #check if a tree distance measured between adjacent trees (k-mer and k-mer+1), else replace to NA
  for (cy1 in seq(1:(length(m_tree_dist)-1))) #dist.topo(cy1, cy1+1)
  {
    if (m_feature_length[[cy1+1]] - m_feature_length[[cy1]]!=1)
    {
      m_tree_dist[cy1]=NA
      
    }
    
  }
  
  
  #unlist all lists and combine into one data.frame for plotting
  m_combined_df <- data.frame(feature_length = unlist(m_feature_length)
                              , file_name = unlist(file_list)
                              #, dist_mat = unlist(m_dist_mat)
                              , jsd_sd = unlist(m_jsd_sd)
                              , jsd_upperlimit = unlist(m_jsd_upperlimit_count)
                              , tree_dist = unlist(m_tree_dist)
  )
  
  return (m_combined_df) #output to dataframe
}
