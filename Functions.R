##Folders
folder_data <- "../WormSynteny/Data/"
folder_ubuntu <- "../WormSynteny/Ubuntu/"

## Functions 

## subset_by_overlap : find overlapping genes 
subset_by_overlap <- function(species, aligned_coord_species){
  
  #Use Granges package : create Granges object first 
  aligned_coord_species_gr <- makeGRangesFromDataFrame(aligned_coord_species, keep.extra.columns = TRUE)
  species_gr <-  makeGRangesFromDataFrame(species, keep.extra.columns = TRUE)
  
  #Apply sybsetOverlaps on Granges objects 
  overlap_species <- subsetByOverlaps(species_gr, aligned_coord_species_gr, type = "any", ignore.strand = TRUE)
  
  #Convert grange objects to data frames 
  overlap_species <- overlap_species %>% as.data.frame()
  
  return(overlap_species)
}

## subset_by_intersect : find the overlapping genes of the C.elegans region in taking only the intersection part
subset_by_intersect <- function(species, aligned_coord_species){
  
  #Use Granges package : create Granges object first 
  aligned_coord_species_gr <- makeGRangesFromDataFrame(aligned_coord_species, keep.extra.columns = TRUE)
  species_gr <-  makeGRangesFromDataFrame(species, keep.extra.columns = TRUE)
  
  #Apply intersect on Granges objects 
  overlapping_genes <- GenomicRanges::intersect(species_gr, aligned_coord_species_gr, ignore.strand = TRUE)
  
  #Convert grange objects to data frames 
  overlapping_genes <- overlapping_genes %>% as.data.frame()
  
  return(overlapping_genes)
}

## create_region : create the region of C.elegans in taking the parameters that the user chose (seq_id, seqnames, start, end, width), check if there is genes in the considered region 
create_region <- function(seqnames, start, end){
  
  ## Step 1 : Create a data frame with the provided region 
  
  ## Region
  seq_id <- paste0("elegans_", seqnames)
  region <- data.frame(seq_id = seq_id, seqnames = seqnames, start = start, end = end, width = (end-start)+1)
  
  #Create a GTF file of the selected region
  region_gr <- makeGRangesFromDataFrame(region, keep.extra.columns = TRUE)
  export(region_gr, paste0(folder_ubuntu, "region.gtf"))
  
  return(region)
  
}

## system_linux : run halLiftover to find aligned region in the 10 other species 
system_linux <- function(species) {
  list_df <- list()
  system(paste0("rm ", folder_ubuntu, "region_*"))
  #Convert singletons_elegans GTF file to BED files
  system(paste0("gtf2bed < ", folder_ubuntu, "region.gtf > ", folder_ubuntu, "region.bed"))
  for(i in species){
    tryCatch({
      system(paste0("/home/lilly/hal/bin/halLiftover ",folder_ubuntu, "evolver_Caenorhabditis.hal --outPSLWithName caenorhabditis_elegans ", folder_ubuntu, "region.bed caenorhabditis_",i, " ",folder_ubuntu,"region_", i, ".psl --bedType 4"))
      list_df[[i]] <- read.table(paste0(folder_ubuntu, "region_", i,".psl"))
    }, error=function(e){print(i)})
  }
  return(list_df)
}

## pslToDataFrame : take a psl file in return a data frame (elegans_species version)
pslToDataFrame <- function(elegans_species, species_name){
  
  if(is.null(elegans_species)){
    return(elegans_species)
  }else{
    #Select the data frame
    colnames(elegans_species) <- c("gene_id", "matches", "misMatches" ,"repMatches","nCount","qNumInsert","qBaseInsert","tNumInsert","tBaseInsert" ,"strand","scaff","qSize","start","end","seq_id2","tSize","start2","end2","blockCount","blockSizes","qStarts","tStarts") 
    
    #Create seq_id and seq_id2 
    elegans_species$seq_id <- paste("elegans",elegans_species$scaff, sep = "_")
    elegans_species <- elegans_species %>% rename(seqnames = "seq_id2")
    elegans_species$seq_id2 <- paste(species_name,elegans_species$seqnames, sep = "_")
    
    #Create strand and strand2 (separate in two columns the column containing both strands)
    elegans_species <- elegans_species %>% mutate(new_strand = gsub("([[:punct:]])([[:punct:]])", "\\1 \\2", elegans_species$strand)) %>% select(-strand) %>% separate(new_strand, into = c("strand", "strand2"), sep = " ")
    return(elegans_species)
  }
}

## system_linux_2 : run halLiftover to find aligned region in adjacent species of the phylogenetic tree 
system_linux_2 <- function(species) {
  list_df <- list()
  for(i in 1 : length(species)){
    tryCatch({
      system(paste0("gtf2bed < ",folder_ubuntu, "reconfigured_singletons_",species[i],".gtf > ", folder_ubuntu,"reconfigured_singletons_",species[i],".bed"))
    }, error=function(e){print(i)})
  }
  for(i in 2 : length(species)-1){
    tryCatch({
      system(paste0("/home/lilly/hal/bin/halLiftover ",folder_ubuntu,"evolver_Caenorhabditis.hal --outPSLWithName caenorhabditis_",species[i]," ", folder_ubuntu,"reconfigured_singletons_",species[i],".bed caenorhabditis_",species[i+1]," ",folder_ubuntu, species[i],"_",species[i+1],".psl --bedType 4"))
    }, error=function(e){print(i)})
    myfile <- paste0(folder_ubuntu,species[i],"_",species[i+1],".psl")
    species1_species2 <- paste0(species[i],"_",species[i+1])
    if (file.exists(myfile)) {
      file_info <- file.info(myfile)
      if(file_info$size != 0){
        list_df[[species1_species2]] <- read.table(paste0(folder_ubuntu,species[i],"_",species[i+1],".psl"))
      }else{
        list_df[[species1_species2]] <- list()
      }
    }else{
      list_df[[species1_species2]] <- list()
    }
  }
  return(list_df)
}

## pslToDataFrame2 : take a psl file in return a data frame (species1_species2 version)
pslToDataFrame2 <- function(list_df, species_names){
  
  species1_species2 <- list_df[[species_names]]
  if(length(species1_species2)==0){
    species1_species2 <- NULL
    return(species1_species2)
  }else{
    #Make it clean
    species <- strsplit(species_names, "_")
    species1_species2 <- list_df[[species_names]]
    colnames(species1_species2) <- c("gene_id","matches", "misMatches" ,"repMatches","nCount","qNumInsert","qBaseInsert","tNumInsert","tBaseInsert" ,"strand","scaff","qSize","start","end","seq_id2","tSize","start2","end2","blockCount","blockSizes","qStarts","tStarts")
    
    #Create seq_id and seq_id2 to match gggenomes requirements 
    species1_species2$seq_id <- paste(  species[[1]][1],species1_species2$scaff, sep = "_")
    species1_species2$seq_id2 <- paste(  species[[1]][2],species1_species2$seq_id2, sep = "_")
    
    return(species1_species2)
  }
}

## raw : do not break the pipeline 
raw <- function(elegans_species){
  if(is.null(elegans_species)){
    return(elegans_species)
  }else{
    raw_aligned_coord <- elegans_species %>%  select(scaff, seq_id2, seqnames, start2, end2, strand2) %>% rename(start = "start2", end = "end2", seq_id = "seq_id2", strand = "strand2") %>% mutate(length = (end - start)+1)
    return(raw_aligned_coord)
  }
}

## stich : join reads close together (<20 b) 
stich <- function(df, stich_df, gap_max){
  i <- 1
  n <- nrow(df) 
  gap <- 0
  gap_cumul <- 0
  l <- list()
  while(n > 0){
    gap <- 0
    stich_df$end <- df$end[n]
    stich_df$start <- df$start[n]
    length <- df$length[n] - 1 
    stich_df$length <- length
    stich_df$gap <- gap
    if(n == 1){
      break
    }else{
      gap <- df$start[n] - df$end[n-1]
      gap_cumul <- gap 
    }
    while(gap < gap_max & gap >= 0){
      stich_df$start <- df$start[n-1]
      length <- (df$length[n-1] -1)  + length
      stich_df$length <- length
      stich_df$gap <- gap_cumul
      n <- n - 1
      if(n < 2){
        break
      }else{
        gap <- df$start[n] - df$end[n-1]
        gap_cumul <- gap + gap_cumul
      } 
    }
    l[[i]] <- stich_df 
    i <- i + 1
    n <- n - 1
  }
  l[[i]] <- stich_df 
  df_stich <- as.data.frame(do.call(rbind, l))
  return(df_stich)
}

## stich_envelop : select a sequence to after apply the stich function so only reads from the same sequence are stich
stich_envelop <- function(raw_coord, gap_max){
  
  if(is.null(raw_coord)){
    return(raw_coord)
  }else{
    l <- list()
    df_distinct <- raw_coord %>% distinct(seqnames)
    for(k in 1 : nrow(df_distinct)){
      df <- raw_coord %>% filter(seqnames %in% df_distinct$seqnames[k]) %>% arrange(seqnames, start)
      stich_df <- df[nrow(df),]
      df_stich <- stich(df, stich_df, gap_max)
      l[[k]] <- df_stich
    }
    df <- as.data.frame(do.call(rbind,l))
    return(df)
  }
}


## filter_by_5 : remove the aligned regions that are less than 5% of length of the selected region 
filter_by_5 <- function(stich_species,region, filter_percentage){
  
  if(is.null(stich_species)){
    return(stich_species)
  }else{
    if(nrow(stich_species)>0){
      filtered_df <- stich_species %>% filter(length >= filter_percentage/100 * region$width)
      return(filtered_df)
    }else{
      return(stich_species)
    }
  }
}


## trie function : to check if the aligned region have a gene (so it won't break the pipeline), keep only raw that have genes 
trie <- function(species, aligned_coord){
  if(is.null(aligned_coord)){
    return(aligned_coord)
  }
  else if(nrow(aligned_coord)==0){
    return(aligned_coord)
  }else{
    genes <- subset_by_overlap(species, aligned_coord)
    if(nrow(genes) == 0){
      return(aligned_coord)
    }else{
      j <- 1
      l <- list()
      for(i in 1:nrow(aligned_coord)){
        genes <- subset_by_overlap(species, aligned_coord[i,])
        if(nrow(genes)!=0){
          l[[j]] <- aligned_coord[i,]
          j <- j + 1 
        }
      }
      aligned_coord <- as.data.frame(do.call(rbind, l))
      return(aligned_coord)
    }
    return(trie_species)
  }
}


## final : take the max length of the trie seqs
final <- function(trie_species){
  if(is.null(trie_species)){
    return(trie_species)
  }
  else if(nrow(trie_species)==0){
    return(trie_species)
  }else{
    final_species <- trie_species %>% filter(length == max(trie_species$length))
    if(nrow(final_species) > 1){
      final_species <- final_species %>% arrange(seqnames) %>% dplyr::slice(1)
      return(final_species)
    }else{
      return(final_species)
    }
  }
}

## other_region : do not break the pipeline
other_region <- function(trie_species, final_species){
  if(is.null(final_species)){
    return(final_species)
  }
  else if(nrow(final_species)==0){
    return(final_species)
  }
  else{
    other_region_species <- anti_join(trie_species, final_species)
    return(other_region_species)
  }
}

## reconfigured : being able to plot the sequence even if there is no genes intersected
reconfigured_fct <- function(species_name, species, raw_species, final_species, region){
  if(is.null(final_species)){
    return(final_species)
  }else{
    genes <- subset_by_intersect(species, final_species)
    if(nrow(genes)!=0){
      reconfigured <- genes %>% group_by(seqnames) %>% mutate(start = min(start), end = max(end)) %>% dplyr::slice(1) %>% mutate(seq_id = paste0(species_name,"_",seqnames)) %>% select(seq_id, seqnames, start, end) %>% mutate(length = (end - start)+1)
    }else{
      raw_species <- raw_species %>% filter(length %in% max(length))
      reconfigured <- raw_species %>% mutate(start = min(start), end = start + region$width) %>% dplyr::slice(1) %>% mutate(length = (end - start)+1) %>% select(seq_id, seqnames, start, end, length)
    }
    return(reconfigured)
  }
}

## export : do not break the pipeline
export_reconfig <- function(reconfigured_species, species_name){
  if(!(is.null(reconfigured_species))){
    reconfigured_gr <- makeGRangesFromDataFrame(reconfigured_species, keep.extra.columns = TRUE)
    folder <- paste0(folder_ubuntu, "reconfigured_singletons_",species_name ,".gtf")
    export(reconfigured_gr, folder) 
  }
}

## links : do not break the pipeline
links_fct <- function(species1_species2){
  if(is.null(species1_species2)){
    return(species1_species2)
  }else{
    links_species <- species1_species2 %>%  select(c( "start", "end","seq_id", "start2", "end2", "seq_id2"))
    return(links_species)
  }
}

## genes : do not break the pipeline
genes_fct <- function(species, final_species, species_name){
  if(is.null(final_species)){
    return(final_species)
  }else{
    genes_species <- subset_by_intersect(species, final_species) %>% mutate(seq_id = paste0(species_name,"_",seqnames)) %>% select(seq_id,  seqnames, start, end) %>% mutate(length = (end - start)+1) %>% arrange(seqnames, start)
    return(genes_species)
  }
}

## control_of_noise_genes: control of repetition of seqs in top and bottom (if same seq find in both, remove the bottom one)
control_of_noise_genes <- function(genes_file, genes_noise){
  if(is.null(genes_file)){
    return(genes_file)
  }else{
    if(nrow(genes_noise)== 0){
      return(genes_noise)
    }
    else{
      genes_noise <- anti_join(genes_noise, genes_file, by = "Orthogroup")
      return(genes_noise)
    }
  }
}

## Find orthogroups : after intersecting genes, the orthogroups info is lost so apply overlap find it back and join it
find_orthogroups <- function(species, genes_species){
  if(is.null(genes_species)){
    return(genes_species)
  }else{
    if(nrow(genes_species) != 0){
      l <- list()
      for(i in 1:nrow(genes_species)){
        genes <- genes_species[i,]
        metadata <- subset_by_overlap(species, genes)
        if(nrow(metadata) > 1){
          OG_max_length <- metadata %>% filter(width == max(width)) %>% arrange(start) %>% dplyr::slice(1)
          genes <- genes %>% mutate(Orthogroup = OG_max_length$Orthogroup)
        }else{
          genes <- genes %>% mutate(Orthogroup = metadata$Orthogroup)
        }
        l[[i]] <- genes
      }
      genes_species <-  data.frame(do.call(rbind, l))
      return(genes_species)
      
    }else{
      genes_species <- genes_species %>% mutate(Orthogroup = NA)
      return(genes_species)
    }
  }
}

## Switch_filtered_genes : switch top and bottom genes when selected 
switch_filtered_genes <-  function(list_species, seqs, genes){
  
  if(nrow(seqs) == 1){
    return(genes)
  }else{
  
    #Select genes and genes_noise 
    genes_noise <- genes %>% filter(str_detect(seq_id, "\\(filtered\\)"))
    
    if(nrow(genes_noise)== 0){
      return(genes)
    }else{
      genes <- genes %>% filter(!(str_detect(seq_id, "\\(filtered\\)")))
      
      #Go through the species name list and switch genes and genes_noise for the corresponding species
      for (k in 1 : length(list_species)) {
        
        species <- list_species[[k]]
        
        #Selected the sequences 
        selected_genes <- genes %>% filter(str_detect(seq_id,species))
        selected_genes_noise <- genes_noise %>% filter(str_detect(seq_id,species)) %>% filter(length == max(length)) %>% dplyr::slice(1) #in case there's two genes with the same length...
        
        #Remove the old ones 
        genes <- anti_join(genes, selected_genes)
        genes_noise <- anti_join(genes_noise, selected_genes_noise)
        
        #Put them back in the opposite files
        genes <- rbind(genes, selected_genes_noise)
        genes_noise <- rbind(genes_noise, selected_genes)
        
      }
      # Order genes and genes_noise
      # Create a vector with the desired order of the species names
      species_order <-  c("bovis", "becei", "panamensis", "inopinata","elegans", "tropicalis", "remanei", "latens", "tribulationis", "briggsae", "nigoni")
      
      # Extract the species name from the seq_id column using regular expressions
      genes$species <- str_extract(genes$seq_id, "[:alpha:]+")
      genes_noise$species <- str_extract(genes_noise$seq_id, "[:alpha:]+")
      
      # Convert the species column to a factor with the desired order
      genes$species <- factor(genes$species, levels = species_order)
      genes_noise$species <- factor(genes_noise$species, levels = species_order)
      
      # Order the rows of the data frame by the species column
      genes <- genes[order(genes$species), ] %>%  select(-species) %>% mutate(seq_id = str_remove(seq_id, " \\(filtered\\)"))
      genes_noise <- genes_noise[order(genes_noise$species), ] %>%  select(-species) %>% mutate(seq_id = str_remove(seq_id, " \\(filtered\\)")) %>%  mutate(filtered =  "(filtered)") %>% unite(seq_id, c("seq_id", "filtered"), sep = " ")
      genes <- rbind(genes, genes_noise)
      
      return(genes)
    }
  }
}

## update_seqs : update sequences when switching genes 
update_seqs <- function(list_species, seqs, genes){
  
  #Select seqs and seqs_noise 
  seqs_noise <- seqs %>% filter(str_detect(seq_id, "\\(filtered\\)"))
  
  if(nrow(seqs_noise)== 0){
    return(seqs)
  }else{
    
    seqs <- seqs %>% filter(!(str_detect(seq_id, "\\(filtered\\)")))
    
    #Select genes and genes_noise 
    genes_noise <- genes %>% filter(str_detect(seq_id, "\\(filtered\\)"))
    genes <- genes %>% filter(!(str_detect(seq_id, "\\(filtered\\)")))
    
    #Go through the species name list 
    for (k in 1 : length(list_species)) {
      
      species <- list_species[[k]]
      
      #Remove species the sequences 
      seqs <- seqs %>% filter(!(str_detect(seq_id,species)))
      seqs_noise <- seqs_noise %>% filter(!(str_detect(seq_id,species)))
      
      #Select the sequences in genes 
      select_genes <- genes %>% filter(str_detect(seq_id,species))
      select_genes_noise <- genes_noise %>% filter(str_detect(seq_id,species))
      
      
      #create the reconfigured sequences 
      select_seqs <- select_genes %>% group_by(seq_id) %>% mutate(start = min(start), end = max(end)) %>% dplyr::slice(1) %>% select(seq_id, seqnames, start, end) %>% mutate(length = (end - start)+1)
      select_seqs_noise <- select_genes_noise  %>% select(seq_id, seqnames,start, end, length) 
      
      #Join everything 
      seqs <- rbind(seqs, select_seqs)
      seqs_noise <- rbind(seqs_noise, select_seqs_noise)
      
    }
    
    #Order seqs
    # Create a vector with the desired order of the species names
    species_order <-  c("bovis", "becei", "panamensis", "inopinata","elegans", "tropicalis", "remanei", "latens", "tribulationis", "briggsae", "nigoni")
    
    # Extract the species name from the seq_id column using regular expressions
    seqs$species <- str_extract(seqs$seq_id, "[:alpha:]+")
    seqs_noise$species <- str_extract(seqs_noise$seq_id, "[:alpha:]+")
    
    # Convert the species column to a factor with the desired order
    seqs$species <- factor(seqs$species, levels = species_order)
    seqs_noise$species <- factor(seqs_noise$species, levels = species_order)
    
    # Order the rows of the data frame by the species column
    seqs <- seqs[order(seqs$species), ] %>%  select(-species) %>% mutate(seq_id = str_remove(seq_id, " \\(filtered\\)"))
    seqs_noise <- seqs_noise[order(seqs_noise$species), ] %>%  select(-species) %>% mutate(seq_id = str_remove(seq_id, " \\(filtered\\)")) %>%  mutate(filtered =  "(filtered)") %>% unite(seq_id, c("seq_id", "filtered"), sep = " ")
    
    #Seqs total
    seqs <- rbind(seqs, seqs_noise)
    
    return(seqs)
  }
}

messages <- function(list_df){
  species <- c("bovis", "becei", "panamensis", "inopinata","tropicalis", "remanei", "latens", "tribulationis", "briggsae", "nigoni")
  l <- list()
  j <- 1
  for (i in species) {
    elegans_species <- list_df[[i]]
    if(is.null(elegans_species)){
      l[[j]] <- i
      j <- j + 1
    }
  }
  if(length(l)!=0){
    df <- as.data.frame(do.call(rbind, l))
  }else{
    df <- data.frame(V1 = NULL)
  }
  return(df)
}

#### FILES 

#Load GTF files with Orthogroups 
bovis <- read.csv( paste0(folder_data, "bovis.csv"))
becei <- read.csv( paste0(folder_data, "becei.csv"))
panamensis <- read.csv( paste0(folder_data, "panamensis.csv"))
inopinata <- read.csv(paste0(folder_data, "inopinata.csv"))
tropicalis <- read.csv( paste0(folder_data, "tropicalis.csv"))
remanei <- read.csv( paste0(folder_data, "remanei.csv"))
latens <- read.csv( paste0(folder_data, "latens.csv"))
tribulationis <- read.csv(paste0(folder_data, "tribulationis.csv"))
briggsae <- read.csv( paste0(folder_data, "briggsae.csv"))
nigoni <- read.csv( paste0(folder_data, "nigoni.csv"))

#Clean it 
bovis <-  bovis  %>% select(-X)
becei <- becei  %>% select(-X)
panamensis <- panamensis  %>% select(-X)
inopinata <- inopinata  %>% select(-X)
tropicalis <- tropicalis  %>% select(-X)
remanei <- remanei  %>% select(-X)
latens <- latens  %>% select(-X)
tribulationis <- tribulationis  %>% select(-X)
briggsae <- briggsae  %>% select(-X)
nigoni <- nigoni  %>% select(-X)

## C.elegans genes 
elegans <- read.csv(paste0(folder_data, "elegans.csv"))
elegans <- elegans %>% select(-X)

#### PIPELINE 
general_function <- function(region, gap_max, filter_percentage){
  
  
  ## Step 2 : Get the coordinates of regions aligned with singletons elegans in the 10 other species (HalLiftover)
  #Species names
  species <- c("bovis", "becei", "panamensis", "inopinata","tropicalis", "remanei", "latens", "tribulationis", "briggsae", "nigoni")
  list_df <- system_linux(species)
  species_no_alignment <- messages(list_df)
  write.csv(species_no_alignment, paste0(folder_ubuntu, "species_no_alignment.csv"))
  
  #becei
  species_name <- "becei"
  elegans_becei <- list_df$becei
  elegans_becei <- pslToDataFrame(elegans_becei,species_name)
  
  #becei
  species_name <- "bovis"
  elegans_bovis <- list_df$bovis
  elegans_bovis <- pslToDataFrame(elegans_bovis,species_name)
  
  #panamensis
  species_name <- "panamensis"
  elegans_panamensis <- list_df$panamensis
  elegans_panamensis <- pslToDataFrame(elegans_panamensis,species_name)
  
  #inopinata
  species_name <- "inopinata"
  elegans_inopinata <- list_df$inopinata
  elegans_inopinata <- pslToDataFrame(elegans_inopinata,species_name)
  
  #tropicalis
  species_name <- "tropicalis"
  elegans_tropicalis <- list_df$tropicalis
  elegans_tropicalis <- pslToDataFrame(elegans_tropicalis,species_name)
  
  #remanei
  species_name <- "remanei"
  elegans_remanei <- list_df$remanei
  elegans_remanei <- pslToDataFrame(elegans_remanei,species_name)
  
  #latens
  species_name <- "latens"
  elegans_latens <- list_df$latens
  elegans_latens <- pslToDataFrame(elegans_latens,species_name)
  
  #tribulationis
  species_name <- "tribulationis"
  elegans_tribulationis <- list_df$tribulationis
  elegans_tribulationis <- pslToDataFrame(elegans_tribulationis,species_name)
  
  #briggsae
  species_name <- "briggsae"
  elegans_briggsae <- list_df$briggsae
  elegans_briggsae <- pslToDataFrame(elegans_briggsae,species_name)
  
  #nigoni
  species_name <- "nigoni"
  elegans_nigoni <- list_df$nigoni
  elegans_nigoni <- pslToDataFrame(elegans_nigoni,species_name)
  
  ## Aligned coordinates of other species with the selected region of C.elegans
  raw_aligned_coord_bovis <- raw(elegans_bovis)
  raw_aligned_coord_becei <- raw(elegans_becei)
  raw_aligned_coord_panamensis <- raw(elegans_panamensis)
  raw_aligned_coord_inopinata <- raw(elegans_inopinata)
  raw_aligned_coord_tropicalis <-  raw(elegans_tropicalis)
  raw_aligned_coord_remanei <-  raw(elegans_remanei)
  raw_aligned_coord_latens <-  raw(elegans_latens)
  raw_aligned_coord_tribulationis <- raw(elegans_tribulationis)
  raw_aligned_coord_briggsae <-  raw(elegans_briggsae)
  raw_aligned_coord_nigoni <-  raw(elegans_nigoni)
  
  ## Step 3 : Stich the small reads together and filter by 5% of the length 
  stich_bovis <- stich_envelop(raw_aligned_coord_bovis, gap_max)
  stich_becei <- stich_envelop(raw_aligned_coord_becei, gap_max)
  stich_panamensis <- stich_envelop(raw_aligned_coord_panamensis, gap_max) 
  stich_inopinata <- stich_envelop(raw_aligned_coord_inopinata, gap_max ) 
  stich_tropicalis <- stich_envelop(raw_aligned_coord_tropicalis, gap_max) 
  stich_remanei <- stich_envelop(raw_aligned_coord_remanei, gap_max) 
  stich_latens <- stich_envelop(raw_aligned_coord_latens, gap_max) 
  stich_tribulationis <- stich_envelop(raw_aligned_coord_tribulationis, gap_max) 
  stich_briggsae <- stich_envelop(raw_aligned_coord_briggsae, gap_max) 
  stich_nigoni <- stich_envelop(raw_aligned_coord_nigoni, gap_max)
  
  ## Filter by 5% of the length of C.elegans
  aligned_coord_bovis <- filter_by_5(stich_bovis,region, filter_percentage) 
  aligned_coord_becei <- filter_by_5(stich_becei,region, filter_percentage)
  aligned_coord_panamensis <- filter_by_5(stich_panamensis,region, filter_percentage)
  aligned_coord_inopinata <- filter_by_5(stich_inopinata,region, filter_percentage)
  aligned_coord_tropicalis <- filter_by_5(stich_tropicalis,region, filter_percentage)
  aligned_coord_remanei <- filter_by_5(stich_remanei,region, filter_percentage)
  aligned_coord_latens <- filter_by_5(stich_latens,region, filter_percentage)
  aligned_coord_tribulationis <- filter_by_5(stich_tribulationis,region, filter_percentage)
  aligned_coord_briggsae <- filter_by_5(stich_briggsae,region, filter_percentage)
  aligned_coord_nigoni <- filter_by_5(stich_nigoni,region, filter_percentage)
  
  # Apply function trie 
  trie_bovis <- trie(bovis, aligned_coord_bovis) 
  trie_becei <- trie(becei, aligned_coord_becei) 
  trie_panamensis <- trie(panamensis, aligned_coord_panamensis) 
  trie_inopinata <- trie(inopinata, aligned_coord_inopinata) 
  trie_tropicalis <- trie(tropicalis, aligned_coord_tropicalis) 
  trie_remanei <- trie(remanei, aligned_coord_remanei) 
  trie_latens <- trie(latens, aligned_coord_latens) 
  trie_tribulationis <- trie(tribulationis, aligned_coord_tribulationis) 
  trie_briggsae <- trie(briggsae, aligned_coord_briggsae) 
  trie_nigoni <- trie(nigoni, aligned_coord_nigoni)
  
  
  ## Step 4 : Take ONLY the max length for displaying 
  final_bovis <- final(trie_bovis)
  final_becei <- final(trie_becei)
  final_panamensis <- final(trie_panamensis)
  final_inopinata <- final(trie_inopinata)
  final_tropicalis <- final(trie_tropicalis)
  final_remanei <- final(trie_remanei)
  final_latens <- final(trie_latens)
  final_tribulationis <- final(trie_tribulationis)
  final_briggsae <- final(trie_briggsae)
  final_nigoni <- final(trie_nigoni)
  
  ## Step 5 : Take the OTHER regions and put them in the bottom files 
  other_regions_bovis <- other_region(trie_bovis, final_bovis) 
  other_regions_becei <- other_region(trie_becei, final_becei) 
  other_regions_panamensis <- other_region(trie_panamensis, final_panamensis) 
  other_regions_inopinata <- other_region(trie_inopinata, final_inopinata) 
  other_regions_tropicalis <- other_region(trie_tropicalis, final_tropicalis) 
  other_regions_remanei <- other_region(trie_remanei, final_remanei) 
  other_regions_latens <- other_region(trie_latens, final_latens) 
  other_regions_tribulationis <- other_region(trie_tribulationis, final_tribulationis) 
  other_regions_briggsae <- other_region(trie_briggsae, final_briggsae) 
  other_regions_nigoni <- other_region(trie_nigoni, final_nigoni)
  
  
  ## Step 6 : Find the sequences length in the aligned regions in the 11 species
  species_name <- "bovis"
  reconfigured_singletons_bovis <- reconfigured_fct(species_name, bovis, raw_aligned_coord_bovis, final_bovis, region)
  
  species_name <- "becei"
  reconfigured_singletons_becei <- reconfigured_fct(species_name, becei, raw_aligned_coord_becei, final_becei, region)
  
  species_name <- "panamensis"
  reconfigured_singletons_panamensis <- reconfigured_fct(species_name, panamensis, raw_aligned_coord_panamensis, final_panamensis, region) 
  
  species_name <- "inopinata"
  reconfigured_singletons_inopinata <- reconfigured_fct(species_name, inopinata, raw_aligned_coord_inopinata, final_inopinata, region)
  
  species_name <- "tropicalis"
  reconfigured_singletons_tropicalis <- reconfigured_fct(species_name, tropicalis,raw_aligned_coord_tropicalis, final_tropicalis, region)
  
  species_name <- "remanei"
  reconfigured_singletons_remanei <- reconfigured_fct(species_name, remanei, raw_aligned_coord_remanei,final_remanei, region) 
  
  species_name <- "latens"
  reconfigured_singletons_latens <- reconfigured_fct(species_name, latens,raw_aligned_coord_latens, final_latens, region)
  
  species_name <- "tribulationis"
  reconfigured_singletons_tribulationis <- reconfigured_fct(species_name, tribulationis,raw_aligned_coord_tribulationis, final_tribulationis, region)
  
  species_name <- "briggsae"
  reconfigured_singletons_briggsae <- reconfigured_fct(species_name, briggsae, raw_aligned_coord_briggsae,final_briggsae, region)
  
  species_name <- "nigoni"
  reconfigured_singletons_nigoni <- reconfigured_fct(species_name, nigoni, raw_aligned_coord_nigoni, final_nigoni, region)
  
  reconfigured_singletons_elegans  <- region %>% mutate(seq_id = paste0("elegans_", seqnames)) %>% select(seq_id, seqnames, start, end, width) %>% dplyr::rename(length = "width")
  
  system(paste0("rm ", folder_ubuntu, "reconfigured_singletons_*"))
  system(paste0("rm ", folder_ubuntu, "bovis_becei.psl"))
  system(paste0("rm ", folder_ubuntu, "becei_panamensis.psl"))
  system(paste0("rm ", folder_ubuntu, "panamensis_inopinata.psl"))
  system(paste0("rm ", folder_ubuntu, "inopinata_elegans.psl"))
  system(paste0("rm ", folder_ubuntu, "elegans_tropicalis.psl"))
  system(paste0("rm ", folder_ubuntu, "tropicalis_remanei.psl"))
  system(paste0("rm ", folder_ubuntu, "remanei_latens.psl"))
  system(paste0("rm ", folder_ubuntu, "latens_tribulationis.psl"))
  system(paste0("rm ", folder_ubuntu, "tribulationis_briggsae.psl"))
  system(paste0("rm ", folder_ubuntu, "briggsae_nigoni.psl"))
  
  species_name <-  "bovis"
  export_reconfig(reconfigured_singletons_bovis, species_name)
  
  species_name <-  "becei"
  export_reconfig(reconfigured_singletons_becei, species_name)
  
  species_name <-  "panamensis"
  export_reconfig(reconfigured_singletons_panamensis, species_name)
  
  species_name <-  "inopinata"
  export_reconfig(reconfigured_singletons_inopinata, species_name)
  
  species_name <-  "tropicalis"
  export_reconfig(reconfigured_singletons_tropicalis, species_name)
  
  species_name <-  "remanei"
  export_reconfig(reconfigured_singletons_remanei, species_name)
  
  species_name <-  "latens"
  export_reconfig(reconfigured_singletons_latens, species_name)
  
  species_name <-  "tribulationis"
  export_reconfig(reconfigured_singletons_tribulationis, species_name)
  
  species_name <-  "briggsae"
  export_reconfig(reconfigured_singletons_briggsae, species_name)
  
  species_name <-  "nigoni"
  export_reconfig(reconfigured_singletons_nigoni, species_name)
  
  species_name <-  "elegans"
  export_reconfig(reconfigured_singletons_elegans, species_name)
  
  ## Step 7 : Use the alignment file to find the alignment between two neighboring species on the phylogenetic tree based on the coordinates above
  species <- c("bovis", "becei", "panamensis", "inopinata","elegans","tropicalis", "remanei", "latens", "tribulationis", "briggsae", "nigoni")
  list_df <- system_linux_2(species)
  
  #Bovis_Becei
  species_name <- "bovis_becei"
  bovis_becei <- pslToDataFrame2(list_df,species_name)
  
  #Becei_Panamensis
  species_name <- "becei_panamensis"
  becei_panamensis <- pslToDataFrame2(list_df,species_name)
  
  #Panamensis_Inopinata
  species_name <- "panamensis_inopinata"
  panamensis_inopinata <- pslToDataFrame2(list_df,species_name)
  
  #Inopinata_Elegans
  species_name <- "inopinata_elegans"
  inopinata_elegans <- pslToDataFrame2(list_df,species_name)
  
  #Elegans_Tropicalis
  species_name <- "elegans_tropicalis"
  elegans_tropicalis <- pslToDataFrame2(list_df,species_name)
  
  #Tropicalis_Remanei
  species_name <- "tropicalis_remanei"
  tropicalis_remanei <- pslToDataFrame2(list_df,species_name)
  
  #Remanei_Latens
  species_name <- "remanei_latens"
  remanei_latens <- pslToDataFrame2(list_df,species_name)
  
  #Latens_Tribulationis
  species_name <- "latens_tribulationis"
  latens_tribulationis <- pslToDataFrame2(list_df,species_name)
  
  #Tribulationis_Briggsae
  species_name <- "tribulationis_briggsae"
  tribulationis_briggsae <- pslToDataFrame2(list_df,species_name)
  
  #Briggsae_Nigoni
  species_name <- "briggsae_nigoni"
  briggsae_nigoni <- pslToDataFrame2(list_df,species_name)
  
  #Select only interesting columns : gene_id, start, end, seq_id, start2, end2,  seq_id2 
  links_bovis_becei <- links_fct(bovis_becei)
  links_becei_panamensis <- links_fct(becei_panamensis)
  links_panamensis_inopinata <- links_fct(panamensis_inopinata)
  links_inopinata_elegans <- links_fct(inopinata_elegans)
  links_elegans_tropicalis <- links_fct(elegans_tropicalis)
  links_tropicalis_remanei <- links_fct(tropicalis_remanei)
  links_remanei_latens <- links_fct(remanei_latens)
  links_latens_tribulationis <- links_fct(latens_tribulationis)
  links_tribulationis_briggsae <- links_fct(tribulationis_briggsae)
  links_briggsae_nigoni <- links_fct(briggsae_nigoni)
  
  ## Step 8 : Find genes overlapping region -> apply again the function subset_by_intersect but without reconfiguring the coordinates
  # Genes
  species_name <- "bovis"
  genes_bovis <- genes_fct(bovis, final_bovis, species_name)
  
  species_name <- "becei"
  genes_becei <- genes_fct(becei, final_becei, species_name)
  
  species_name <- "panamensis"
  genes_panamensis <- genes_fct(panamensis, final_panamensis, species_name)
  
  species_name <- "inopinata"
  genes_inopinata <- genes_fct(inopinata, final_inopinata, species_name)
  
  species_name <- "tropicalis"
  genes_tropicalis <- genes_fct(tropicalis, final_tropicalis, species_name)
  
  species_name <- "remanei"
  genes_remanei <- genes_fct(remanei, final_remanei, species_name)
  
  species_name <- "latens"
  genes_latens <- genes_fct(latens, final_latens, species_name)
  
  species_name <- "tribulationis"
  genes_tribulationis <- genes_fct(tribulationis, final_tribulationis, species_name)
  
  species_name <- "briggsae"
  genes_briggsae <- genes_fct(briggsae, final_briggsae, species_name)
  
  species_name <- "nigoni"
  genes_nigoni <- genes_fct(nigoni, final_nigoni, species_name)
  
  genes_elegans <- subset_by_intersect(elegans, region) %>% mutate(seq_id = paste0("elegans_",seqnames)) %>% select(seq_id,  seqnames, start, end) %>% mutate(length = (end - start)+1) %>% arrange(seqnames, start)
  
  #Apply function
  genes_bovis <- find_orthogroups(bovis, genes_bovis)
  genes_becei <- find_orthogroups(becei, genes_becei)
  genes_panamensis <- find_orthogroups(panamensis, genes_panamensis)
  genes_inopinata <- find_orthogroups(inopinata, genes_inopinata)
  genes_tropicalis <- find_orthogroups(tropicalis, genes_tropicalis)
  genes_remanei <- find_orthogroups(remanei, genes_remanei)
  genes_latens <- find_orthogroups(latens, genes_latens)
  genes_tribulationis <- find_orthogroups(tribulationis, genes_tribulationis)
  genes_briggsae <- find_orthogroups(briggsae, genes_briggsae)
  genes_nigoni <- find_orthogroups(nigoni, genes_nigoni)
  genes_elegans <- find_orthogroups(elegans, genes_elegans)
  
  #Genes noise
  species_name <- "bovis"
  genes_bovis_noise <- genes_fct(bovis, other_regions_bovis, species_name)
  
  species_name <- "becei"
  genes_becei_noise <- genes_fct(becei, other_regions_becei, species_name)
  
  species_name <- "panamensis"
  genes_panamensis_noise <- genes_fct(panamensis, other_regions_panamensis, species_name)
  
  species_name <- "inopinata"
  genes_inopinata_noise <- genes_fct(inopinata, other_regions_inopinata, species_name)
  
  species_name <- "tropicalis"
  genes_tropicalis_noise <- genes_fct(tropicalis, other_regions_tropicalis, species_name)
  
  species_name <- "remanei"
  genes_remanei_noise <- genes_fct(remanei, other_regions_remanei, species_name)
  
  species_name <- "latens"
  genes_latens_noise <- genes_fct(latens, other_regions_latens, species_name)
  
  species_name <- "tribulationis"
  genes_tribulationis_noise <- genes_fct(tribulationis, other_regions_tribulationis, species_name)
  
  species_name <- "briggsae"
  genes_briggsae_noise <- genes_fct(briggsae, other_regions_briggsae, species_name)
  
  species_name <- "nigoni"
  genes_nigoni_noise <- genes_fct(nigoni, other_regions_nigoni, species_name)
  
  ## Find orthogroups
  genes_bovis_noise <- find_orthogroups(bovis, genes_bovis_noise)
  genes_becei_noise <- find_orthogroups(becei, genes_becei_noise)
  genes_panamensis_noise <- find_orthogroups(panamensis, genes_panamensis_noise)
  genes_inopinata_noise <- find_orthogroups(inopinata, genes_inopinata_noise)
  genes_tropicalis_noise <- find_orthogroups(tropicalis, genes_tropicalis_noise)
  genes_remanei_noise <- find_orthogroups(remanei, genes_remanei_noise)
  genes_latens_noise <- find_orthogroups(latens, genes_latens_noise)
  genes_tribulationis_noise <- find_orthogroups(tribulationis, genes_tribulationis_noise)
  genes_briggsae_noise <- find_orthogroups(briggsae, genes_briggsae_noise)
  genes_nigoni_noise <- find_orthogroups(nigoni, genes_nigoni_noise)
  
  ## Step 9 : control of repetition of seqs in top and bottom
  genes_bovis_noise <- control_of_noise_genes(genes_bovis, genes_bovis_noise)
  genes_becei_noise <- control_of_noise_genes(genes_becei, genes_becei_noise)
  genes_panamensis_noise <- control_of_noise_genes(genes_panamensis, genes_panamensis_noise)
  genes_inopinata_noise <- control_of_noise_genes(genes_inopinata, genes_inopinata_noise)
  genes_tropicalis_noise <- control_of_noise_genes(genes_tropicalis, genes_tropicalis_noise)
  genes_remanei_noise <- control_of_noise_genes(genes_remanei, genes_remanei_noise)
  genes_latens_noise <- control_of_noise_genes(genes_latens, genes_latens_noise)
  genes_tribulationis_noise <- control_of_noise_genes(genes_tribulationis, genes_tribulationis_noise)
  genes_briggsae_noise <- control_of_noise_genes(genes_briggsae, genes_briggsae_noise)
  genes_nigoni_noise <- control_of_noise_genes(genes_nigoni, genes_nigoni_noise)
  
  print("before 10")
  
  ## Step 10 : Create seqs, genes and links files for gggenomes
  #Links
  links <- rbind(links_bovis_becei, links_becei_panamensis, links_panamensis_inopinata, links_inopinata_elegans,links_elegans_tropicalis, links_tropicalis_remanei,links_remanei_latens,  links_latens_tribulationis,links_tribulationis_briggsae, links_briggsae_nigoni)
  if(!(is.null(links))){
    links <- links %>% mutate(seq_id = str_replace(seq_id, "_"," "), seq_id2 = str_replace(seq_id2, "_"," "))
  }else{
    links <- NULL
    
  }
  
  #Genes
  genes_noise <-  rbind(genes_bovis_noise, genes_becei_noise, genes_panamensis_noise, genes_inopinata_noise, genes_tropicalis_noise, genes_remanei_noise, genes_latens_noise, genes_tribulationis_noise, genes_briggsae_noise, genes_nigoni_noise)
  
  if(!(is.null(genes_noise))){
    if(nrow(genes_noise) > 0){
      
      #Hundle genes noise 
      list_noise <- list(genes_bovis_noise, genes_becei_noise, genes_panamensis_noise, genes_inopinata_noise, genes_tropicalis_noise, genes_remanei_noise, genes_latens_noise, genes_tribulationis_noise, genes_briggsae_noise, genes_nigoni_noise)
      genes_noise <- data.frame()
      
      #Loop over the data frames
      for (i in 1:10) {
        if(!(is.null(list_noise[[i]]))){
          if(nrow(list_noise[[i]])> 15 ){
            new_df <- list_noise[[i]] %>% group_by(seqnames) %>% arrange(desc(start)) %>% slice(1:15)
            genes_noise <- rbind(genes_noise, new_df)
          }else{
            list_noise[[i]]$Orthogroup <- as.character(list_noise[[i]]$Orthogroup)
            genes_noise <- rbind(genes_noise, list_noise[[i]])
          }
        }
      }
      genes_noise <- genes_noise %>%  mutate(filtered =  "(filtered)") %>% unite(seq_id, c("seq_id", "filtered"), sep = " ")
      
      # create a list of all the data frames
      df_list <- list(genes_bovis, genes_becei, genes_panamensis, genes_inopinata, genes_elegans, genes_tropicalis, genes_remanei, genes_latens, genes_tribulationis, genes_briggsae, genes_nigoni)
      df_reconfig <- list(reconfigured_singletons_bovis, reconfigured_singletons_becei, reconfigured_singletons_panamensis, reconfigured_singletons_inopinata, reconfigured_singletons_elegans, reconfigured_singletons_tropicalis, reconfigured_singletons_remanei, reconfigured_singletons_latens, reconfigured_singletons_tribulationis, reconfigured_singletons_briggsae, reconfigured_singletons_nigoni)
      species_name <-  c("bovis", "becei", "panamensis", "inopinata","elegans", "tropicalis", "remanei", "latens", "tribulationis", "briggsae", "nigoni")
      
      # initialize an empty data frame to store the combined data
      genes <- data.frame()
      reconfigured <- data.frame()
      
      # loop over the data frames
      for (i in 1:11) {
        if(!(is.null(df_list[[i]]))){
          # check the number of rows in the data frame
          if (nrow(df_list[[i]]) > 15) {
            # Find the region interval in that species taking the min genes and the length of elegans region
            new_df <- df_list[[i]] %>% filter(start == min(df_list[[i]]$start)) %>% mutate(end = start + region$width, length = region$width) %>% select(start, end, length, seqnames)
            # Keep the overlapping genes in that region
            selection <- subset_by_overlap(df_list[[i]], new_df) %>% select(seq_id, seqnames, start, end, length, Orthogroup)
            genes <- rbind(genes, selection)
            
            # Same for reconfigured 
            reconfig <- selection %>% group_by(seqnames) %>% mutate(start = min(start), end = max(end)) %>% dplyr::slice(1) %>% mutate(seq_id = paste0(species_name[i],"_",seqnames)) %>% select(seq_id, seqnames, start, end) %>% mutate(length = (end - start)+1)
            reconfigured <- rbind(reconfigured, reconfig)
          } else {
            # bind the data frame to the existing data frame if it has 100 or fewer rows
            genes <- rbind(genes, df_list[[i]])
            reconfigured <- rbind(reconfigured,  df_reconfig[[i]])
            
          }
        }
      }
      genes <- rbind(genes, genes_noise)
      genes <- genes %>% mutate(seq_id = str_replace(seq_id, "_"," "))
      
    }else{
      
      # create a list of all the data frames
      df_list <- list(genes_bovis, genes_becei, genes_panamensis, genes_inopinata, genes_elegans, genes_tropicalis, genes_remanei, genes_latens, genes_tribulationis, genes_briggsae, genes_nigoni)
      df_reconfig <- list(reconfigured_singletons_bovis, reconfigured_singletons_becei, reconfigured_singletons_panamensis, reconfigured_singletons_inopinata, reconfigured_singletons_elegans, reconfigured_singletons_tropicalis, reconfigured_singletons_remanei, reconfigured_singletons_latens, reconfigured_singletons_tribulationis, reconfigured_singletons_briggsae, reconfigured_singletons_nigoni)
      species_name <-  c("bovis", "becei", "panamensis", "inopinata","elegans", "tropicalis", "remanei", "latens", "tribulationis", "briggsae", "nigoni")
      
      # initialize an empty data frame to store the combined data
      genes <- data.frame()
      reconfigured <- data.frame()
      
      # loop over the data frames
      for (i in 1:11) {
        if(!(is.null(df_list[[i]]))){
          # check the number of rows in the data frame
          if (nrow(df_list[[i]]) > 15) {
            # Find the region interval in that species taking the min genes and the length of elegans region
            new_df <- df_list[[i]] %>% filter(start == min(df_list[[i]]$start)) %>% mutate(end = start + region$width, length = region$width) %>% select(start, end, length, seqnames)
            # Keep the overlapping genes in that region
            selection <- subset_by_overlap(df_list[[i]], new_df) %>% select(seq_id, seqnames, start, end, length, Orthogroup)
            genes <- rbind(genes, selection)
            
            # Same for reconfigured 
            reconfig <- selection %>% group_by(seqnames) %>% mutate(start = min(start), end = max(end)) %>% dplyr::slice(1) %>% mutate(seq_id = paste0(species_name[i],"_",seqnames)) %>% select(seq_id, seqnames, start, end) %>% mutate(length = (end - start)+1)
            reconfigured <- rbind(reconfigured, reconfig)
          } else {
            # bind the data frame to the existing data frame if it has 100 or fewer rows
            genes <- rbind(genes, df_list[[i]])
            reconfigured <- rbind(reconfigured,  df_reconfig[[i]])
            
          }
        }
      }
      genes <- genes %>% mutate(seq_id = str_replace(seq_id, "_"," "))
      
    }
  }else{
    genes <- NULL
    reconfigured <- rbind(reconfigured_singletons_bovis, reconfigured_singletons_becei, reconfigured_singletons_panamensis, reconfigured_singletons_inopinata, reconfigured_singletons_elegans, reconfigured_singletons_tropicalis, reconfigured_singletons_remanei, reconfigured_singletons_latens, reconfigured_singletons_tribulationis, reconfigured_singletons_briggsae, reconfigured_singletons_nigoni) 
  }
  
  #Seqs normal
  #reconfigured <- rbind(reconfigured_singletons_bovis, reconfigured_singletons_becei, reconfigured_singletons_panamensis, reconfigured_singletons_inopinata, reconfigured_singletons_elegans, reconfigured_singletons_tropicalis, reconfigured_singletons_remanei, reconfigured_singletons_latens, reconfigured_singletons_tribulationis, reconfigured_singletons_briggsae, reconfigured_singletons_nigoni) 
  reconfigured$seqnames <- as.character(reconfigured$seqnames)
  
  seqs <- reconfigured %>% mutate(seq_id = str_replace(seq_id, "_"," "))
  
  #Seqs noise
  if(!(is.null(genes_noise))){
    seqs_noise <- genes_noise  %>% select(seq_id, seqnames,start, end, length) %>% mutate(seq_id = str_replace(seq_id, "_"," "))
  }else{
    seqs_noise <-  NULL
  }
  
  #Seqs total
  seqs <- rbind(seqs, seqs_noise)
  
  print("after 10")
  return(list(seqs = seqs, genes = genes, links = links))
  
}
