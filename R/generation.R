

#' Profiles generation.
#' 
#' Function to generate all possible combinations of attribute levels (i.e. all 
#' possible profiles).
#' 
#' Valid arguments for \code{coding} are \code{C}, \code{D} and \code{E}. When
#' using \code{C} the attribute will be treated as continuous and no coding will
#' be applied. All possible levels should then be specified in \code{c.lvls}. If
#' \code{D} (dummy coding) is used \code{\link{contr.treatment}} will be applied
#' to that attribute. For \code{E} (effect coding) \code{\link{contr.sum}} will
#' be applied.
#' 
#' @param lvls  A numeric vector which contains for each attribute the number
#'   of levels.
#' @param coding Type of coding that needs to be used for each attribute.
#' @param c.lvls A list containing numeric vectors with the attribute levels for
#'   each continuous attribute. The default is \code{NULL}.
#' @return A numeric matrix which contains all possible profiles.
#' @examples 
#' # Without continuous attributes
#' at.lvls <- c(3, 4, 2) # 3 Attributes with respectively 3, 4 and 2 levels. 
#' c.type <- c("E", "E", "E") # All Effect coded.
#' Profiles(lvls = at.lvls, coding = c.type) # Generate profiles.
#' 
#' # With continuous attributes 
#' at.lvls <- c(3, 4, 2) # 3 attributes with respectively 3, 4 and 2 levels. 
#' # First attribute is dummy coded, second and third are continuous. 
#' c.type <- c("D", "C", "C") 
#' # Levels for continuous attributes, in the same order. 
#' con.lvls <- list(c(4, 6, 8, 10), c(7, 9))
#' Profiles(lvls = at.lvls, coding = c.type, c.lvls = con.lvls)
#' @export
Profiles <- function(lvls, coding, c.lvls = NULL) {
  # Continuous attributes. 
  contins <-  which(coding == "C")
  n.contins <-  length(contins)
  # error continuous levels 
  if (!is.null(c.lvls)){ 
    if(!is.list(c.lvls)){stop("'c.lvls' should be a list.")}
    if(!is.numeric(unlist(c.lvls))){stop("components of 'c.lvls' should be numeric")}
  }
  # Error correct coding types.
  codings.types <- c("E", "D", "C")
  if (!all(coding %in% codings.types) || (length(coding) != length(lvls))) {
    stop("coding argument is incorrect.")
  } 
  # Error lvls vector.
  if (length(lvls) < 2 || (!(is.numeric(lvls)))){
    stop("lvls argument is incorrect.")
  }
  # Error continuous specified and NULL.
  if (length(contins) > 0 && is.null(c.lvls)) {
    stop("when 'coding' contains C, 'c.lvls' should be specified")
  }
  # Error continuous levels specification. 
  if (!is.null(c.lvls)) {
    if (length(c.lvls) != n.contins) {
      stop("length of 'c.lvls' does not match number of specified continuous attributes in 'coding'")
    }
    # Error c.lvls same number of levels. 
    if (!isTRUE(all.equal(lvls[contins], lengths(c.lvls)))) {
      stop("the number of levels provided in 'c.lvls' does not match the expected based on 'lvls'")
    }
  }
  # Change into correct coding. 
  coding <- dplyr::recode(coding, D = "contr.treatment", E = "contr.sum")
  # Create all combinations of attribute levels.
  levels.list <- lapply(X = as.list(lvls), function(x) (1:x))
  # Replace continuous.
  levels.list[contins] <- c.lvls
  # Create grid. 
  dgrid <- as.data.frame(expand.grid(levels.list))
  # Apply coding to non continuous. 
  cn <- names(dgrid)
  if (!is.null(c.lvls)) {
    cn <- cn[-contins]
  }
  # Create factors.
  dgrid[, cn] <- apply(dgrid[, cn, drop = FALSE], 2, factor)
  # coding 
  con <- as.list(stats::setNames(coding, names(dgrid)))
  con[which(con == "C")] <- NULL
  cgrid <- as.data.frame(stats::model.matrix(~., dgrid, contrasts = con))
  # Delete intercept.
  cgrid <- cgrid[, -1]
  # Return profiles.
  return(as.matrix(cgrid))
}

# Create alternative specific coding.
Altspec <- function(alt.cte, n.sets) {
  if(!any(alt.cte == 0)){
    stop("'alt.cte' should at least contain 1 zero")
  }
  # create matrix
  mat <- diag(length(alt.cte))
  n.zero <- which(alt.cte == 0)
  mat[n.zero, n.zero] <- 0
  # delete zero columns
  del.col <- c(which(apply(mat, 2,   function(x) all(x == 0))))
  mat <- mat[, -del.col]
  #rbind for full design 
  mat <- as.matrix(mat)
  cte.mat <- do.call(rbind, replicate(n.sets, mat, simplify = FALSE)) 
  #return
  return(cte.mat)
}

# Create row and column names for designs 
Rcnames <- function(n.sets, n.alts, alt.cte, no.choice) {
  # rownames
  r.s <- rep(1:n.sets, each = n.alts)
  r.a <- rep(1:n.alts, n.sets)
  r.names <- paste(paste("set", r.s, sep = ""), paste("alt", r.a, sep = ""), sep = ".")
  if(no.choice){
    ncsek <- seq(n.alts, (n.sets * n.alts), n.alts)  
    r.names[ncsek] <- "no.choice"
  }
  # colnames alternative specific constants
  if(sum(alt.cte) > 0.2){
    cte.names <- paste(paste("alt", which(alt.cte == 1), sep = ""), ".cte", sep = "") 
  } else {
    cte.names <- NULL
  }
  # return
  return(list(r.names, cte.names))
}

# Lattice multivariate standard normal distribution.
# 
# Generates a grid of points coming from a multivariate standard normal
# distribution.
# @param K Numeric value indicating the dimensionality of the grid.
# @param b Numeric value indicating the base.
# @param m Numeric value. Number of draws = b^m.
# @return Matrix of lattice points drawn from a multivariate standard normal
#   distribution. Each row is a draw.
Lat <- function(K, b, m) {
  base <- function(num){
    a1 <- c1 <- rep(0, m)
    a1[m] <- c1[m] <- num %% b
    for(j in seq(m - 1, 1, by = -1)){
      tem <- (num - c1[j + 1]) / (b^(m - j))
      a1[j] <- tem %% b
      c1[j] <- c1[j + 1] + a1[j] * (b^(m - j))
    }
    return(a1)
  }
  a <- 1571
  N <- b^m   #number of lattice points
  u <- stats::runif(K)
  av <- rep(1, K)
  for (i in 2:K) {
    av[i] <- (a * av[i - 1]) %% N
  }
  e <- matrix(0, N, K)
  seq <- rep(0, m)
  kk <- -m
  for (k in 1:m) {
    seq[k] <- b^kk
    kk <- kk + 1
  }
  for (i in 1:N) {
    ei <- c(crossprod(seq, base(i - 1))) * av + u
    e[i, ] <- ei - floor(ei)
  }
  latt <- matrix(0, N, K)
  for (i in 1:K) {
    for (j in 1:N) {
      latt[j, i] <- 1 - abs(2 * e[j, i] - 1)
    }
  }
  latt <- stats::qnorm(latt)
  return(latt)
}

# Lattice multivariate t-distribution.
# 
# Generates a grid of points coming from a multivariate t-distribution.
# @param mean Numeric vector indicating the multivariate mean.
# @param cvar A matrix which specifies the covariance matrix.
# @param df Numeric value indicating the degrees of freedom for the multivariate
#   t-distribution.
# @param m Numeric value. Number of draws = b^m.
# @param b Numeric value indicating the base (default = 2).
# @return Matrix of lattice points drawn from a multivariate t-distribution.
#   Each row is a draw.
Lattice_mvt <- function (mean, cvar, df, m, b=2) {
  # Dimension
  dim <- length(mean)
  # Generate lattice from standard normal
  lattice <- Lat(K = dim, b, m)
  mean <- t(mean)
  X <- matrix(NA, nrow(lattice), dim)
  A <- chol(cvar)
  # Transform to multivariate t distribution
  for (i in 1:nrow(lattice)) {
    inv <- stats::rgamma(1, df / 2)
    invw <- inv / (df / 2)
    W <- 1 / invw
    Z <- lattice[i, ]
    r <- mean + sqrt(W) * (Z %*% t(A))
    X[i, ] <- t(r)
  }
  return (X)
}

# Lattice multivariate normal distribution.
# 
# Generates a grid of points coming from a multivariate normal distribution.
# @param mean Numeric vector indicating the multivariate mean.
# @param cvar A matrix which specifies the covariance matrix.
# @param m Numeric value. Number of draws = b^m.
# @param b Numeric value indicating the base (default = 2).
# @return Matrix of lattice points drawn from a multivariate normal
#   distribution. Each row is a draw.
Lattice_mvn <- function(mean, cvar, m, b=2) {
  dim <- length(mean)
  l <- Lat(K = dim, b, m)
  mean <- t(mean)
  X <- matrix(NA, nrow(l), dim)
  A <- chol(cvar)
  for (i in 1:nrow(l)) {
    Z <- l[i, ]
    r <- mean + (Z %*% t(A))
    X[i, ] <- t(r)
  }
  return(X)
}


## generate grid for truncated distribution
Lattice_trunc <- function (n, mean, cvar, lower, upper, df) {
  # Dimension
  dim <- length(mean)
  # Generate lattice from standard normal
  left <- n
  mean <- t(mean)
  A <- chol(cvar)
  XX <- NULL
  while(left > 0.2){
    m = ceiling(log(left)/log(2))
    if(m < 2){m = 2}
    lattice <- Lat(K = dim, b = 2, m = m)
    X <- matrix(NA, nrow(lattice), dim)
    # Transform to multivariate t distribution
    for (i in 1:nrow(lattice)) {
      Z <- lattice[i, ]
      r <- mean + (Z %*% t(A))
      X[i, ] <- t(r)
    }
    XS <- matrix(unlist(apply(X, 1, function(X) {if(all(lower < X) && all(X < upper)) return(X)})), ncol = dim, byrow = TRUE)
    XX <- rbind(XX, XS)
    left <- (n - nrow(XX))
  }
  return(XX[1:n, ])
}

##list of all possible choice sets 
Fullsets <- function(cand.set, n.alts, no.choice, reduce = TRUE){
  
  if(!is.null(no.choice)){
    n.alts <- n.alts - 1
  }
  full.comb <- utils::combn(1:nrow(cand.set), n.alts, FUN = function(x)  cand.set[x, ], simplify = FALSE)
  #reduce
  if (reduce){
    m <- stats::rnorm(ncol(cand.set))
    inf <-list()
    for(i in 1:length(full.comb)){
      inf[[i]] <- round(InfoDes(m, full.comb[[i]], n.alts), digits = 3)
    }
    t <- array(unlist(inf), dim = c(length(m), length(m), length(inf))) 
    full.comb <- full.comb[!duplicated(t, MARGIN = 3)]
  }
  if(!is.null(no.choice)){
    full.comb <- lapply(full.comb, Inchoice, no.choice = no.choice)
  }
  return(full.comb)
}

# when opt.out = TRUE in modfed 
# Optout <- function(des, n.alts, n.sets){
#   optdes <- cbind(rep(0, n.alts * n.sets), des)
#   no.choice <- c(1, rep(0, ncol(optdes) - 1))
#   design.opt <- vector(mode = "list")
#   for (s in 1: n.sets){
#     design.opt[[s]] <- rbind(optdes[ (((s-1) * n.alts) + 1) : (s * n.alts) , ], no.choice)
#   }
#   optdes <- do.call(rbind, design.opt)
#   colnames(optdes)[1] <- "no.choice.cte"
#   return(optdes)
# }

Optout <- function(des, n.alts, alt.cte, n.sets){
  
  n.cte <- sum(alt.cte)
  names <- colnames(des)
  names <- append(names, "no.choice.cte", n.cte)
  
  if(n.cte > 0.2){
    stripdes <- des[ , -(1:n.cte)]
  } else {
    stripdes <- des 
  }
  
  no.choice <-  rep(0, ncol(stripdes))
  design.opt <- vector(mode = "list")
  for (s in 1: n.sets){
    design.opt[[s]] <- rbind(stripdes[ (((s-1) * n.alts) + 1) : (s * n.alts) , ], no.choice)
  }
  optdes <- do.call(rbind, design.opt)
  cte.des <- Altspec(c(alt.cte, 1), n.sets)
  optdes <- cbind(cte.des, optdes)
  optdes
  colnames(optdes) <- names
  return(optdes)
}

## include no.choice
Inchoice <- function(X , no.choice) {
  if(no.choice > nrow(X)){
    X <- rbind(X, rep(0, ncol(X)))
  } else {
    X <- rbind(X, rep(0, ncol(X)))
    X[seq(no.choice + 1, nrow(X)), ] <- X[seq(no.choice, nrow(X) - 1), ]
    X[no.choice, ] <- rep(0, ncol(X))
  }
  return(X)
}


#helper function to specify to which variable each column belongs to, in a matrix resulting from the function Profiles
categorize_variables <- function(original_vector) {
  # Initialize the result vector with the original vector
  result <- original_vector
  names(result) <- original_vector
  
  #if the names do not start with Var as is the output of Profiles, then the names are converted
  if (!all(grepl("Var", original_vector,ignore.case = TRUE))) {
    prefixes <- gsub("[0-9]+$", "", original_vector)
    unique_prefixes <- unique(prefixes)
    prefix_map <- stats::setNames(paste0("var", seq_along(unique_prefixes)), unique_prefixes)
    suffixes <- gsub("^[A-Za-z]+", "", original_vector)
    original_vector <- paste0(prefix_map[prefixes], suffixes)
  }
  
  # Loop through each element of the original vector
  for (i in 1:length(original_vector)) {
    # Remove any trailing '.1' (for continuous variables that can result from the Profiles function)
    vector <- gsub("\\.1$", "", original_vector[i])
    
    # Extract the numeric part after 'Var'
    num_part <- as.numeric(gsub("Var", "", vector, ignore.case = TRUE))
    
    # If the number part is longer than the index, truncate the last digit
    if (nchar(i) < nchar(num_part)) {
      num_part <- as.numeric(substr(num_part, 1, nchar(num_part) - 1))
    }
    
    # Process based on the number part
    if (num_part > 1) {
      # If the difference with the previous variable number is more than 1
      if (num_part - as.numeric(result[i - 1]) > 1) {
        # Increment the previous variable number
        result[i] <- as.numeric(result[i - 1]) + 1
        # The variable number should correspond to the truncated number part, if not then take the number for the previous variable
        result[i] <- ifelse(as.numeric(result[i]) - as.numeric(substr(num_part, 1, nchar(num_part) - 1)), 
                            as.numeric(result[i]) - 1, result[i])
      } else {
        result[i] <- num_part
      }
    } else {
      result[i] <- num_part
    }
  }
  return(result)
}

#helper function to create a choice set with specified overlap
create_overlapping_cs <- function(cand.set, overlap, n.alts) {
  choice_set = matrix(nrow = n.alts, ncol = ncol(cand.set), dimnames = list(c(rep(0,n.alts)), c(colnames(cand.set)))) 
  
  #first we select a random profile for the first alternative
  first_profile <- cand.set[sample(nrow(cand.set), 1), , drop = FALSE]
  choice_set[1, ] <- first_profile
  rownames(choice_set) = rep(rownames(first_profile)[1],n.alts)
  
  #here we determine which attributes to overlap
  num_attributes <- as.numeric(max(categorize_variables(colnames(cand.set))))
  overlap_var <- sample(num_attributes, overlap)
  overlap_indices <- which(categorize_variables(colnames(cand.set)) %in% overlap_var) #extend the filtering to all the relevant columns for the overlapped variables (in case of dummy/effect treatments with more than 2 levels)
  
  #Filter cand.set for profiles matching the first profile in overlapping attributes
  filtered_cand_set = cand.set
  for (index in overlap_indices) {
    filtered_cand_set = filtered_cand_set[filtered_cand_set[, index] == first_profile[index], , drop=FALSE] #drop=FALSE preserves the df structure
  }
  
  #If filtered_cand_set is empty or too small, then we choose to go back to the original cand.set
  if (nrow(filtered_cand_set) < n.alts - 1) {
    filtered_cand_set = cand.set
  }
  
  #lastly we fill the rest of the choice set from the filtered candidate set
  for (alt in 2:n.alts) {
    r = sample(nrow(filtered_cand_set), 1)
    choice_set[alt, ] = filtered_cand_set[ r , ,drop=F]
    rownames(choice_set)[alt] = rownames(filtered_cand_set[ r , ,drop=F])
  }
  
  return(choice_set)
}


#constraints####
#the functions here use regexp, the meanings of some terms can be consulted on the webpage:
#https://cran.r-project.org/web/packages/stringr/vignettes/regular-expressions.html

#helper function to evaluate a design against user-specified constraints
check_constraints <- function(desje, constraints, set_row_start, n.cte) {
  
  #this is just to save the original constraint to be used later for checking spacing between the terms
  original_constraints = constraints
  
  #translate logical operators
  constraints = gsub("AND", "&", constraints)
  constraints = gsub("OR", "|", constraints)
  constraints = gsub(" = ", " %in% ", constraints)
    
  # if `if` statement is used, then first we break it down into two parts and
  # handle them separately
  if (grepl("^if", constraints, ignore.case = TRUE)) {
    # Split constraints into 'if' part and 'then' part
    parts <- strsplit(constraints, "then", fixed = TRUE)[[1]]
    if_part <- trimws(parts[1])
    then_part <- trimws(parts[2])
    
    # Remove the initial 'if' from the if_part
    if_part2 <- sub("^if ", "", if_part, ignore.case = TRUE)
    
    if (identical(if_part, if_part2)) {
      stop(paste0("Please use spaces between the terms in the following constraint: ",if_part2,"."))
    }
    
    # Evaluate the 'if' part, and if it is met, then the 'then' part is evaluated
    if_result <- check_constraints(desje, if_part2, set_row_start, n.cte)
    if (if_result) {
      return(check_constraints(desje, then_part, set_row_start, n.cte))
    } else {
      # If the 'if' condition is not met, return TRUE without evaluating the 'then' part (that is because the condition for the constraint is not met)
      return(TRUE)  
    }
  } else {
  # Split the condition based on spaces to handle expressions
  parts <- unlist(strsplit(constraints, " "))
  
  # added a check for the proper spacing between terms, there should be at least
  # 1 operator in the terms, and an extra operator per each 'AND'/'OR' in the constraint
  needed_operators = 1 + sum(parts %in% c("&", "|"))
  actual_operators = sum(parts %in% c("%in%", "!=", "<", "<=",">",">="))
  if (actual_operators < needed_operators) {
    stop(paste0("Please use spaces between the terms and the operators in the following constraint: ",original_constraints,"."))
  }
  
  
  # initialize an expression string
  expression_string <- ""
  # this is used in case there are "=" or "!=" in the constraint
  matching <- NULL
  # Loop through parts to replace AltX.AttY with the correct indices
  for (part in parts) {
    if (grepl("&|\\|", part)) {
      
      if (matching) {
        expression_string <- paste(expression_string, ")", sep = "")
      }
      intermediate <- (sapply(expression_string, function(x) eval(parse(text=x))))
      expression_string <- paste(intermediate,part, sep = " ")
      matching = FALSE
      
    } else if (grepl("Alt", part)) {
      #here we extract the alt and att numbers
      alt <- as.numeric(sub(".*Alt(\\d+).*", "\\1", part)) + set_row_start - 1
      att <- as.numeric(sub(".*Att(\\d+).*", "\\1", part))
      #we use an index to extract the affected columns based on the 'att' number (since multiple columns could belong to 1 attribute)
      index <- as.character(which(categorize_variables(colnames(desje[,(n.cte+1):ncol(desje)]) ) %in% att) + n.cte)
      #concatenate indices into a single string
      index <- paste(index, collapse = ",")
      
      #append the correct index reference to the expression string
      expression_string <- paste(expression_string, paste0("desje[",alt,", c(", index, ")]"), sep = " ")
      
    } else if (grepl("%in%|!=", part))  {
      #to handle %in% or !=, we replace them with the function match_values as defined above
      # and we do the evaluation for this part right away, to avoid adding match_values again if 
      #these operators occur again in the constraints expression
      
      #find all occurrences of "&" or "|"
      all_separators <- gregexpr("&|\\|", expression_string)
      
      # Get the last occurrence of "&" or "|"
      last_separator <- lapply(all_separators, function(x) {if (length(x) > 0) max(x) else 0})
      last_separator <- unlist(last_separator)
      # Insert "match_values(" after the last "&" or "|"
      for (i in seq_along(last_separator) ) {
        if (last_separator[i] > 0) {
          expression_string[i] <- paste0(
          substr(expression_string[i], 1, last_separator[i]),          # Everything up to and including the "&" or "|"
          " match_values(",                                                                 # Insert 'match_values('
          substr(expression_string[i], last_separator[i] + 1, nchar(expression_string[i]))  # Everything after the "&" or "|"
          ,",")
        } else {
          expression_string[i] <- paste("match_values(",expression_string[i],",", sep = "")
        }
      }
        matching = TRUE
      
        if (grepl("!=", part)) {
          #here we turn match_values(...) into !match_values(...), to negate the matching
          #the (?<!!) expression is to check that 'match_values' is not already preceded by an exclamation mark,
          #so as to avoid adding it again
          expression_string <- sub("(?<!!)match_values", "!match_values", expression_string, perl=TRUE)
        }
    } else {
      # Append operators and numbers directly to the expression string
      expression_string <- paste(expression_string, part, sep = " ")
    }
   }
  }
  
  #add the last ')' if there is matching -> to produce the function match_values(...)
  if (!is.null(matching) && matching) {
    expression_string <- paste(expression_string, ")", sep = "")
  }
  #finally we evaluate the constructed expression
  return(all(sapply(expression_string, function(x) eval(parse(text=x)))))
}

#function to check value constraints (similar to %in%, but works for vectors as well)
match_values <- function(check_attribute, possible_values) {
  if (length(check_attribute) > 1) {
    if (!is.list(possible_values)) {
      possible_values <- list(possible_values)
    }
    
    any(sapply(possible_values, function(vec) all(vec == check_attribute)))
  } else {
    check_attribute %in% possible_values
  }
}


#helping function to extract the number of all the attributes included in constraints
ext_constrained_atts <- function(constraints) {
  attribute_numbers <- numeric()
  for (constraint in constraints) {
    #find all occurrences of 'Att' followed by one or more digits
    matches <- regmatches(constraint, gregexpr("Att(\\d+)", constraint))
    # Extract matches and convert them to numeric and append the result
    attr_nums <- as.numeric(unlist(lapply(matches, function(x) sub("Att", "", x))))
    attribute_numbers <- c(attribute_numbers, attr_nums)
  }
  # the unique attribute numbers are returned
  return(unique(attribute_numbers))
}

#helping function to extract the number of all the alternative included in constraints
ext_constrained_alts <- function(constraints) {
  alternative_numbers <- numeric()
  for (constraint in constraints) {
    #find all occurrences of 'Alt' followed by one or more digits
    matches <- regmatches(constraint, gregexpr("Alt(\\d+)", constraint))
    # Extract matches and convert them to numeric and append the result
    alt_nums <- as.numeric(unlist(lapply(matches, function(x) sub("Alt", "", x))))
    alternative_numbers <- c(alternative_numbers, alt_nums)
  }
  # the unique alternative numbers are returned
  return(unique(alternative_numbers))
}

#helper function to translate the dummy coding (for example 1A or 3C) to the corresponding dummy/effect coding in constraints.
translate_coded_constraints <- function(constraint, coding = NULL, lvls = NULL, cand.set = NULL) {
  #in case Modfed is used, then coding and lvls are deduced
  if (is.null(coding) && is.null(lvls)) {
    coding = vector()
    lvls = vector()
    # categorize_variables(colnames(desje[,(n.cte+1):ncol(desje)]) )
    categorized_vars = as.numeric(categorize_variables(colnames(cand.set) ))
    as.numeric(categorize_variables(colnames(cand.set) ))
    
    for (v in 1:max(categorized_vars) ) {
      i <- which(v == categorized_vars)
    
      if (all(cand.set[,i] %in% c(-1,0,1) ) ) {
        
        if (!any(cand.set[,i] %in% -1)) {
          coding[v] = "D"
          lvls[v] = length(cand.set[1,i]) + 1
        } else {
          coding[v] = "E"
          lvls[v] = length(cand.set[1,i]) + 1
        }
      } else {
        coding[v] = "C"
        lvls[v] = length(unique(cand.set[,i]))
        }
    }
  }

  #If there are no dummy/effect coded variables, then do nothing
  '%nin%' = Negate('%in%')
  if (all(coding %nin% c("D", "E")) ) {
    return(constraint)
  }
  
  #here we use functions that work with regex to extract the contents (levels) 
  #in the specified list in the constraint
  cleaned_constraint = regmatches(constraint, gregexpr("list\\(([^)]+)\\)", constraint))
  #we remove everything except the contents inside list()
  cleaned_constraint = unlist(lapply(cleaned_constraint, function(x) gsub("list\\(|\\)", "", x)))
  #if nothing is found then do not change the constraint expression
  if (length(cleaned_constraint) == 0) {
    return(constraint)
  }
  
  # cleaned_constraint <- gsub(".*list\\((.*)\\).*", "\\1", constraint)
  # Split the cleaned string by commas to get the individual parts
  parts = vector()
  for (p in 1:length(cleaned_constraint)) {
    parts = c(parts, strsplit(cleaned_constraint[p], ",\\s*")[[1]])
  }
  
  # Initialize
  coding_list <- list()
  
  # Generate coding for each specified level
  for (part in parts) {
    # part = parts[1]
  # for (level in specified_levels) {
    
    #extract numbers first
    num_part = gsub("[^0-9]", "", part)
    #extract letters
    letter_part = gsub("[0-9]", "", part)
    
    # some checks for some invalid scenarios (when num or letter parts are empty)
    if (any(c(num_part,letter_part) %in% "")) {
      next
    }
    if (num_part == 0) {
      next
    }
    if (letter_part %nin% c(LETTERS,paste0(LETTERS,LETTERS),paste0(LETTERS,LETTERS,LETTERS))) {
      next
    }
    
    #
    var_coding = coding[as.numeric(num_part)]
    var_levels = lvls[as.numeric(num_part)]
    
    #in case there are more than 26 levels we extend to levels AA, AB,..., AAA, AAB...ZZZ
    if (var_levels > 26) {
      levels <- c(LETTERS,paste0(LETTERS,LETTERS),paste0(LETTERS,LETTERS,LETTERS))[1:var_levels]
      } else {
      levels <- c(LETTERS,paste0(LETTERS,LETTERS))[1:var_levels]
    }

    # we put a check here that the chosen level is indeed within the possible levels of the variable
    if (letter_part %nin% levels) {
      stop("Invalid levels specified in constraints. Possibly, the level specified is higher than the maximum level of the constrained variable.")
    } else {
      level_index <- which(levels == letter_part) - 1
    }
    
    #inisitailize a verctor of zeros
    coding_vector = numeric(var_levels - 1)
    
    # Set the corresponding index in coding_vector to 1
    if (!is.na(level_index) && length(level_index) > 0) {
      if (var_coding == "D") {
        coding_vector[level_index] <- 1
      } else if (var_coding == "E") {
          if (level_index == 0) {
            coding_vector <- coding_vector - 1
          } else {
            coding_vector[level_index] <- 1
          }
      }
    }
    
    # Add the coding vector to the list with correct naming
    coding_list[[part]] <- paste0("c(",paste(coding_vector,collapse = ","),")")
  }
  
  # Check if the list is valid
  if (length(coding_list) == 0) {
    return(constraint)
  }
  
  # Apply replacements
  expression <- constraint
  for (level in names(coding_list) ) {
    pattern <- paste0("\\b", level, "\\b")  # \b for word boundary to match whole word only
    expression <- gsub(pattern, coding_list[[level]], expression, perl = TRUE)
  }
  
  return(expression)
}


