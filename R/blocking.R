#' Create blocks (sub-designs) from a given design.
#'
#' This function breaks down a design output from \code{Modfed} or \code{CEA} into a specified number of blocks
#' while aiming to maintain balance in levels frequency across the resulting blocks.
#' 
#' The argument \code{n.blocks} specifies the number of blocks to create. The algorithm strives to distribute the choice sets of the design evenly among the blocks, while maintaining level balance across them. The choice sets are assigned sequentially to the blocks, aiming to maintain the closest possible level balance among them up to that stage in the sequence. Hence, the algorithm runs different iterations, during each of which the choice sets in the design are shuffled randomly. The argument \code{blocking.iter} specifies the maximum number of these iterations.
#' 
#' If the design has a no.choice alternative then \code{no.choice} should be set to \code{TRUE}. Additionally, \code{asc.col} should indicate the number of alternative specific constants that are included in the design, if any.
#' 
#' This functionality is also available as an argument (\code{n.blocks}) when creating an efficient design using 
#' \code{\link{Modfed}} or \code{\link{CEA}}.
#' 
#' Note: To make sure the code works well, the names of the variables in the provided
#' design should be aligned with variable names that the function \code{Profiles} produces. For
#' example, if attribute 1 is a dummy variable of 3 levels then its corresponding columns
#' should have numbered names such as: var11 and var12, or (if labelled) price1 and price2, for instance.
#' 
#' @param des The design to be distributed into blocks (sub-designs).
#' @param n.blocks A numeric value indicating the desired number of blocks to
#'   create out of the provided design.
#' @param n.alts The number of alternatives in each choice set.
#' @param blocking.iter A numeric value indicating the maximum number of iterations 
#'   for optimising the level balance in the blocks. The default value is 50.
#' @param no.choice A logical value indicating whether a no choice alternative 
#'   is added to each choice set in the provided design. The default is \code{FALSE}.
#' @param alt.cte A binary vector indicating for each alternative whether an 
#'   alternative specific constant is present in the design. The default is \code{NULL}.
#' @return A list of blocks from the original design is returned. Additionally, the 
#'   frequency of every level in each block is returned.
#' @export
#'
#' @examples
#' \donttest{
#' # DB-efficient designs
#' # 3 Attributes with 3 levels, all dummy coded. 1 alternative specific constant = 7 parameters
#' cand.set <- Profiles(lvls = c(3, 3, 3), coding = c("D", "D", "D"))
#' mu <- c(0.5, 0.8, 0.2, -0.3, -1.2, 1.6, 2.2) # Prior parameter vector
#' v <- diag(length(mu)) # Prior variance.
#' set.seed(123) 
#' pd <- MASS::mvrnorm(n = 10, mu = mu, Sigma = v) # 10 draws.
#' p.d <- list(matrix(pd[,1], ncol = 1), pd[,2:7])
#' design <- Modfed(cand.set = cand.set, n.sets = 8, n.alts = 2, 
#'        alt.cte = c(1, 0), parallel = FALSE, par.draws = p.d)
#' Blocks(design$BestDesign$design, n.blocks = 2, n.alts = 2, alt.cte = c(1, 0))
#' }
Blocks <- function(des, n.blocks, n.alts, blocking.iter = 50,
                   no.choice = FALSE, alt.cte = NULL) {
  
  if (is.null(alt.cte)) {
    alt.cte <- rep(0L, n.alts)
  }
  #init
  n.cte <- length(which(alt.cte == 1))
  ### Errors
  #handling alt.cte
  if (length(alt.cte) != n.alts) {
    stop("'n.alts' does not match the 'alt.cte' vector")
  }
  if (!all(alt.cte %in% c(0, 1))) {
    stop("'alt.cte' should only contain zero or ones.")
  }
  #if no.choice
  if (!is.logical(no.choice)) {
    stop("'no.choice' should be TRUE or FALSE")
  }
  if (no.choice) {
    if (!isTRUE(all.equal(alt.cte[n.alts], 1))) {
      stop("if 'no.choice' is TRUE, alt.cte[n.alts] should equal 1.")
    }
  } 
  
  #if no.choice
  if (!is.logical(no.choice)) {
    stop("'no.choice' should be TRUE or FALSE")
  }
  
  ###
  
  `%nin%` = Negate(`%in%`)
  design2 = des[,(n.cte + 1):(ncol(des))] #to exclude asc columns
  design_counts = lvl.freq2(design2)
  
  n_sets = nrow(des) / n.alts
  
  if (no.choice) {
    ncsek <- seq(n.alts, (n_sets * n.alts), n.alts) 
  } else {
    ncsek <- NULL
  }
  
  design_countperset <- lapply(design_counts, function(counts) as.vector(counts)/n_sets)
  sets_per_block <- ceiling(n_sets / n.blocks)
  final_blocks <- vector("list", length = n.blocks)
  
  for(it in 1:blocking.iter) {
    
    #initialize empty data frames
    blocks = vector("list", length = n.blocks)
    
    for (i in 1:n.blocks) {
      blocks[[i]] <- data.matrix(data.frame(matrix(NA, nrow = 0, ncol = ncol(des))))
      colnames(blocks[[i]]) <- colnames(des)
    }
    
    #shuffle sets for random assignment
    shuffled_sets_indices <- sample(1:n_sets)
    
    for (index in shuffled_sets_indices) {
      # index=shuffled_sets_indices[1]
      set_rows <- (index - 1) * n.alts + (1:n.alts)
      set <- des[set_rows, , drop = FALSE]
      
        # Calculate the criterion for each block in case the current set is added to each block individually
        #in our implementation it is (the count in the block - the count in the desing)^2, summed for all attributes
        criterion <- sapply(blocks, function(block) {
          combined_block <- rbind(block, set)
          
            freq <- lvl.freq2(combined_block[,(n.cte + 1):(ncol(des))])
            block_count <- unlist(freq) #unlisting the count in the current block
            total_count <- unlist(lvl.freq2(design2)) #unlisting the count in the original design
            
            #not all levels might be present in the current block so far
            present_lvls <- which(names(total_count) %in% names(block_count))
            non_present_lvls <- which(names(total_count) %nin% names(block_count))
            total_count[present_lvls] <- block_count #replace the present lvls with the count in the block
            total_count[non_present_lvls] <- 0 #replace the non-present lvls with zeros
            block_countperset <- total_count/(nrow(combined_block)/n.alts ) #adjust the lvls by the number of sets in the block (lvl count per set)
            sum(((block_countperset) - unlist(design_countperset))^2 ) #calculate levels between the block counts per set and the original design counts per set
        })
        
        #assign set to the block with the minimum criterion, but only to the blocks with least amount of sets so far
        smallest_blocks <- which(unlist(lapply(blocks,nrow)) == min(unlist(lapply(blocks,nrow))))
        min_criterion_block <- which.min(criterion[smallest_blocks])
        chosen_block <- smallest_blocks[min_criterion_block]
        blocks[[chosen_block]] <- rbind(blocks[[chosen_block]], set)
      
    }
    
    #calculate overall criterion for this iteration, after all sets has been assigned to the blocks
    
    overall_criterion <- mean(sapply(blocks, function(block) {
        
    block_count <- unlist(lvl.freq2(block[,(n.cte + 1):(ncol(des))])) #unlisting the count in the current block so far
    total_count <- unlist(lvl.freq2(design2)) #unlisting the count in the original design
        
    #not all levels might be present in the current block so far
    present_lvls <- which(names(total_count) %in% names(block_count))
    non_present_lvls <- which(names(total_count) %nin% names(block_count))
    total_count[present_lvls] <- block_count #replace the present lvls with the count in the block
    total_count[non_present_lvls] <- 0 #replace the non-present lvls with zeros
    block_countperset <- total_count/(nrow(block)/n.alts ) #adjust the lvls by the number of sets in the block (lvl count per set)
    sum(((block_countperset) - unlist(design_countperset))^2 ) #calculate ssd between the block counts per set and the original design counts per set
      }))
     
    
    if (it == 1 || (!is.nan(overall_criterion) && overall_criterion < current_criterion)) {
      current_criterion <- overall_criterion
      final_blocks <- blocks
    }
  }
  
  #reorder the sets in each block
  for (i in 1:n.blocks) {
    if (no.choice) {
      ncsek_block <- ncsek[1:(nrow(final_blocks[[i]])/n.alts)]
      order_full <- sort(rownames(final_blocks[[i]]), index.return = TRUE)$ix
      ordered_indices <- order_full[-which(order_full %in% ncsek_block)] #without the no.choice indices
      for (nc in (ncsek_block)) {
        ind <- which(nc == ncsek_block)
        ordered_indices = append(ordered_indices, nc, after=n.alts*ind-1)
      }
    } else {
        ordered_indices <- sort(rownames(final_blocks[[i]]), index.return = TRUE)$ix
    }
    
    final_blocks[[i]] <- final_blocks[[i]][ordered_indices, ]
  }

  #calculate level frequencies in each block
  lvl.freq <- lapply(final_blocks, function(x) {
    if (no.choice) {
      ncsek_block <- ncsek[1:(nrow(final_blocks[[i]])/n.alts)]
      lvl.freq(x[-ncsek_block,(n.cte + 1):(ncol(des))])
    } else {
       lvl.freq(x[,(n.cte + 1):(ncol(des))]) 
    }
  })   
  lvl.freq_final <- lvl.freq[[1]]
  n <- length(lvl.freq_final)
  for (a in 1:n) {
    for (i in 2:n.blocks) {
      temp <- dplyr::bind_rows(lvl.freq_final[[a]], lvl.freq[[i]][[a]])
      temp[is.na(temp)] <- 0 #replace NAs with zeros
      lvl.freq_final[[a]] <- temp
    }
    rownames(lvl.freq_final[[a]]) <- paste0("block", 1:n.blocks)
  }
  
  #we add a note for non-divisable cases
  if (n_sets %% n.blocks != 0) {
    block_sizes <- unlist(lapply(final_blocks, function(x) nrow(x) / n.alts))
    message(sprintf("%d sets is not divisible by %d. The resulting blocks have %s sets.",
                    n_sets, n.blocks, paste(block_sizes, collapse = ", ")))
  }
  
  ##
  return(list("Blocks"= final_blocks, "levels.count(blocks)" = lvl.freq_final))
}


##levels frequency
lvl.freq <- function(des){
  attributes_index <- categorize_variables(colnames(des))
  num_attributes <- length(unique(attributes_index))
  freq <- vector('list')
  for(a in unique(attributes_index) ) {
    index <- which( a == as.numeric(attributes_index) )
    att_filtered <- des[,index,drop=F]
    mod_att_filtered <- att_filtered[,1]
    
    if (ncol(att_filtered) > 1){
      for (i in 2:ncol(att_filtered)) {
        mod_att_filtered <- paste(mod_att_filtered, att_filtered[,i],sep=",")
      }
      mod_att_filtered
    }
    
    level_names <- paste0("(",sort(unique(mod_att_filtered)),")")
    att_filtered_freq <- table(mod_att_filtered)
    names(attr(att_filtered_freq,"dimnames")) <-NULL #to clear the header from the frequency table
    attr(att_filtered_freq,"names") <- level_names
    freq[[a]] <- as.data.frame(t(data.matrix(att_filtered_freq)))
  }
  names(freq) <- paste0("Attribute", seq_along(freq))
  return(freq)
}


#this a slightly altered function from the previous one, that works better to construct the blocks
lvl.freq2 <- function(des){
  attributes_index <- categorize_variables(colnames(des))
  num_attributes <- length(unique(attributes_index))
  freq <- vector('list')
  for(a in unique(attributes_index) ) {
    index <- which( a == as.numeric(attributes_index) )
    temp <- des[,index,drop=F]
    temp2 <- temp[,1]
    
    if (ncol(temp) > 1) {
      for (i in 2:ncol(temp)) {
        temp2 <- paste(temp2, temp[,i],sep="_")
      }
    } 
    freq[[a]] <- table(temp2)
  }
  return(freq)
}

# counting overlap
lvl.overlap <- function(des, n.alts){
  
  attributes_index <- categorize_variables(colnames(des))
  output <- data.frame(set = numeric(), overlap.count = numeric())
  
  for (i in 1:nrow(des)) {
    set_index <- ((i - 1) %/% n.alts) + 1
    if (i > 1 && (set_index == output[nrow(output), 'set'])) next #from the 2nd iteration, if we're in the same set then skip to the next one
    set_row_start <- (set_index - 1) * n.alts + 1
    set_row_end <- set_index * n.alts
    
    choice_set <- des[set_row_start:set_row_end,,drop=F]
    overlap_count <- 0
    for(a in unique(attributes_index) ) {
      index <- which( a == as.numeric(attributes_index) )
      comparisons <- rep(TRUE, nrow(choice_set)) 
      for (r in 1:nrow(choice_set) ) {
        comparisons[r] <- all(choice_set[1, index] == choice_set[r, index])
        }
      overlapping_attrib <- sum(all(comparisons == TRUE))
      overlap_count <- overlap_count + overlapping_attrib
    }
    output[set_index,'set'] <- set_index
    output[set_index,'overlap.count'] <- overlap_count
  }
  return(output)
}
