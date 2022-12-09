#######################################################################
##     get_entrez_ids
##  
##
#######################################################################

get_entrez_ids = function (sig_names, proper_names) {
  
  # Join the data and the proper names.
  df = left_join(sig_names, proper_names, by = c("Original_name"))
  
  # remove NAs caused by results with no proper name (e.g. from multi-protein antibodies)
  df = na.omit(df)
  
  # Run query on the unique names to find the entrez ID for each gene. 
  ent_sym = queryMany(unique(df$Fixed_name), scopes="symbol", fields="entrezgene", species="human", entrezonly = TRUE)
  
  # join the enterez IDs with the original data.
  ent_sym$Fixed_name = ent_sym$query
  df = inner_join(df, ent_sym, by = "Fixed_name", copy = TRUE)
  
  # remove NAs.
  if (any(is.na(df$entrezgene))) {warning(paste("The genes/proteins", paste(df$Original_name[is.na(df$entrezgene)], collapse = ","),
                                                "could not be assigned entrez IDs and will be ignored" ))}
  df = df[is.na(df$entrezgene) == FALSE,]
  
  # Extract the useful data.
  df = subset(df, select = c(Fixed_name, entrezgene))
  df = dplyr::rename(df, Name = Fixed_name)
  
  return(df)
}