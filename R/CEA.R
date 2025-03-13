
#' Coordinate Exchange algorithm for MNL models.
#' 
#' The algorithm improves an initial start design by considering changes on an
#' attribute-by-attribute basis. By doing this, it tries to minimize the chosen error 
#' (A(B) or D(B)-error) based on a multinomial logit model. This routine is repeated for
#' multiple starting designs.
#' 
#' Each iteration will loop through all profiles from the initial design, 
#' evaluating the change in A(B) or D(B)-error (as specified) for every level 
#' in each attribute. The algorithm stops when an iteration occurred without 
#' replacing a profile or when \code{max.iter} is reached.
#' 
#' By specifying a numeric vector in \code{par.draws}, the A- or D-error will be 
#' calculated and the design will be optimised locally. By specifying a matrix, 
#' in which each row is a draw from a multivariate distribution, the AB/DB-error 
#' will be calculated, and the design will be optimised globally. Whenever there
#' are alternative specific constants, \code{par.draws} should be a list 
#' containing two matrices: The first matrix containing the parameter draws for
#' the alternative specific constant parameters. The second matrix containing
#' the draws for the rest of the parameters.
#' 
#' The AB/DB-error is calculated by taking the mean over A- / D-errors, respectively. 
#' It could be that for some draws the design results in an infinite error. 
#' The percentage of draws for which this was true for the final design can be 
#' found in the output \code{inf.error}.
#' 
#' Alternative specific constants can be specified in \code{alt.cte}. The length
#' of this binary vector should equal \code{n.alts}, were \code{0} indicates the
#' absence of an alternative specific constant and \code{1} the opposite.
#' 
#' \code{start.des} is a list with one or several matrices corresponding to 
#' initial start design(s). In each matrix, each row is a profile. The number of rows
#' equals \code{n.sets * n.alts}, and the
#' number of columns equals the number of columns of the design matrix + the
#' number of non-zero elements in \code{alt.cte}. Consider that for a 
#' categorical attribute with *p* levels, there are *p - 1* columns in the design
#' matrix, whereas for a continuous attribute there is only one column. If
#' \code{start.des = NULL}, \code{n.start} random initial designs will be
#' generated. If start designs are provided, \code{n.start} is ignored.
#' 
#' Note: To make sure the code works well, the names of the variables in the starting 
#' design should be aligned with variable names that the function \code{Profiles} produces. For
#' example, if attribute 1 is a dummy variable of 3 levels then its corresponding columns
#' should have numbered names such as: var11 and var12, or (if labelled) price1 and price2, for instance.
#' 
#' If \code{no.choice} is \code{TRUE}, in each choice set an alternative with
#' one alternative specific constant is added. The return value of the
#' A(B) or D(B)-error is however based on the design without the no choice option.
#' 
#' When \code{parallel} is \code{TRUE}, \code{\link[parallel]{detectCores}} will
#' be used to decide upon the number of available cores. That number minus 1 
#' cores will be used to search for efficient designs. The computation time will
#' decrease significantly when \code{parallel = TRUE}.
#' 
#' **Partial profiles/overlapping attributes**
#' 
#' If \code{overlap} is set to 1 or more, then partial profiles will be used in 
#' the resulting efficient designs. The value of \code{overlap} determines the minimum 
#' number of attributes to overlap in each choice set. The optimising algorithm will 
#' enforce this constraint across all choice sets. Note that the running time may increase 
#' significantly, as the algorithm searches through all possible (combinations of) attributes 
#' to achieve optimisation.
#' 
#' **Blocking**
#' 
#' If the value of \code{n.blocks} is more than 1, a new list with the specified number of blocks of the best design (one with the least A(B)- or D(B)-error) will be added to the output. The algorithm strives to distribute the choice sets of the best design evenly among the blocks, while maintaining level balance across them. The choice sets are assigned sequentially to the blocks, aiming to maintain the closest possible balance among them up to that stage in the sequence. Hence, the algorithm runs different iterations, during each of which the choice sets in the design are shuffled randomly. The argument \code{blocking.iter} specifies the maximum number of these iterations.
#' This functionality is also available as a separate function in \code{\link{Blocks}} that works with a given design.
#' 
#' 
#' **Adding constraints to the design**
#' 
#' The argument \code{constraints} can be used to determine a list of constraints to be enforced on the resulting efficient design.
#' The package offers flexibility in the possible constraints. The basic syntax for the constraint should determine an attribute Y within an alternative X (\code{AltX.AttY}) and an operator to be applied on that attribute followed by a list of values or another attribute. In addition to this basic syntax, conditional If statements can be included in the conditions as will be shown in the examples below.
#' The following operators can be used:
#' 
#'  - \code{=}
#'  - \code{!=}
#'  - \code{<} or \code{<=}
#'  - \code{>} or \code{>=}
#'  - \code{AND}
#'  - \code{OR}
#'  - +, -, *, / operations for continuous attributes.
#' 
#' For example, if attributes 1, 2 and 3 are continuous attributes, then possible constraints include:
#'  * \code{"Alt2.Att1 = list(100, 200)"}: restrict values of attribute 1 in alternative 2 to 100 and 200.
#'  * \code{"Alt1.Att1 > Alt2.Att1"}: enforce that attribute 1 in alternative 1 to be higher than the attribute's value in alternative 2.
#'  * \code{"Alt1.Att1 + Alt1.Att2 < Alt1.Att3"}: enforce that the sum of attributes 1 and 2 to be less than the value of attribute 3 in alternative 1.
#'  * \code{"Alt1.Att1 > Alt1.Att3 OR Alt1.Att2 > Alt1.Att3"}: either attribute 1 or attribute 2 should be higher than attribute 3 in alternative 1.
#' 
#' For dummy and effect coded attributes, the levels are indicated with the number of the attribute followed by a letter from the alphabet. For example \code{1A} is the first level of attribute 1 and \code{3D} is the fourth level of attribute 3. Examples on constraints with dummy/effect coded variables:
#'  * \code{"Alt2.Att1 = list(1A,1B)"}: restrict attribute 1 in alternative 2 to the reference level (A) and the second level (B).
#'  * \code{"Alt1.Att1 = list(1B,1C) AND Alt2.Att2 != list(2A, 2E)"}: restrict attribute 1 in alternative 1 to the second and third levels, and at the same time, attribute 2 in alternative 2 cannot be the first and fifth levels of the attribute.
#' 
#' Additionally, and as aforementioned, conditional If statements can be included in the conditions. Examples:
#'  * \code{"if Alt1.Att1 != Alt2.Att1 then Alt2.Att2 = list(100,200)"}
#'  * \code{"if Alt1.Att1 = Alt2.Att1 OR Alt1.Att1 = 0 then Alt2.Att1 > 3"}
#' 
#' Lastly, more than one constraint can be specified at the same time. For example: \code{constraints = list("if Alt1.Att1 != Alt2.Att1 then Alt2.Att2 = list(100,200)", "Alt1.Att3 = list (3A, 3C)")}. 
#' 
#' To ensure the best use of constraints in optimising designs, please keep in mind the following:
#'  - Proper spacing should be respected between the terms, to make sure the syntax translates properly into an \code{R} code. To clarify, spaces should be placed before and after the operators listed above. Otherwise, the console will return an error.
#'  - Lists should be used for constrained values as shown in the examples above.
#'  - Constraints should not be imposed on the \code{no.choice} alternative because it is fixed with zeros for all attributes. The \code{no.choice} alternative, if included, will be the last alternative in every choice set in the design. Therefore, if \code{no.choice} is \code{TRUE} and the no.choice alternative number (= \code{n.alts}) is included in the constraints, the console will return an Error.
#'  - Attention should be given when a starting design that does not satisfy the constraint is provided. It is possible that the algorithm might not find a design that is more efficient and, at the same time, that satisfies the constraints.
#'  -	With tight constraints, the algorithm might fail to find a design that satisfies all the specified constraints.
#' 
#' 
#' @param lvls  A numeric vector which contains for each attribute the number
#'   of levels.
#' @param coding Type of coding that needs to be used for each attribute.
#' @param c.lvls A list containing numeric vectors with the attribute levels for
#'   each continuous attribute. The default is \code{NULL}.
#' @param n.sets Numeric value indicating the number of choice sets.
#' @param n.alts Numeric value indicating the number of alternatives per choice 
#'   set.
#' @param par.draws A matrix or a list, depending on \code{alt.cte}.
#' @param alt.cte A binary vector indicating for each alternative whether an 
#'   alternative specific constant is desired. The default is \code{NULL}.
#' @param no.choice A logical value indicating whether a no choice alternative 
#'   should be added to each choice set. The default is \code{FALSE}.
#' @param start.des A list containing one or more matrices corresponding to 
#' initial start design(s). The default is \code{NULL}.
#' @param parallel Logical value indicating whether computations should be done 
#'   over multiple cores. The default is \code{TRUE}.
#' @param max.iter A numeric value indicating the maximum number allowed 
#'   iterations. The default is \code{Inf}.
#' @param n.start A numeric value indicating the number of random start designs
#'   to use. The default is 12.
#' @param optim A character value to choose between "D" and "A" optimality. The 
#'    default is \code{"D"}.
#' @param overlap A numeric value indicating the minimum number of attributes to overlap 
#'    in every choice sets to create partial profiles. The default is \code{NULL}.
#' @param n.blocks A numeric value indicating the desired number of blocks to
#'   create out of the most efficient design.
#' @param blocking.iter A numeric value indicating the maximum number of iterations 
#'   for optimising the blocks. The default value is 50.
#' @param constraints A list of constraints to enforce on the attributes and alternatives 
#'   in every choice set. The default is \code{NULL}. 
#' @return 
#'   Two lists of designs and statistics are returned: First, the list \code{BestDesign} contains the design with the lowest A(B)- or D(B)- error. The method \code{print} can be used to return this list. 
#'   Second, the list \code{AllDesigns} contains the results of all (provided) start designs. The method \code{summary} can be used to return this list.
#'   \item{design}{A numeric matrix wich contains an efficient design.} 
#'   \item{optimality}{\code{"A"} or \code{"D"}, depending on the chosen optimality criteria.} 
#'   \item{inf.error}{Numeric value indicating the percentage of draws for which the D-error was \code{Inf}.} 
#'   \item{probs}{Numeric matrix containing the probabilities of each alternative 
#'   in each choice set. If a sample matrix was provided in \code{par.draws}, 
#'   this is the average over all draws.}
#'   \item{AB.error}{Numeric value indicating the A(B)-error of the design.} 
#'   \item{DB.error}{Numeric value indicating the D(B)-error of the design.} 
#'   \item{SD}{The standard deviations of the parameters. Calculated by taking the diagonal of the varcov matrix, averaged over all draws if a sample matrix was provided in \code{par.draws}.} 
#'   \item{level.count}{The count of all levels of each attribute in the design.} 
#'   \item{level.overlap}{The count of overlapping levels accross alternatives in every choice set in the design.} 
#'   \item{Orthogonality}{Numeric value indicating the degree of orthogonality of the design. The closer the value to 1, the more orthogonal the design is.} 
#'   \item{Blocks}{A list showing the created blocks of the best design, along with the level counts in each block. For more details, see function \code{\link{Blocks}}.} 
#' @examples
#' \donttest{
#' # DB-efficient designs
#' # 3 Attributes, all dummy coded. 1 alternative specific constant = 7 parameters
#' mu <- c(1.2, 0.8, 0.2, -0.3, -1.2, 1.6, 2.2) # Prior parameter vector
#' v <- diag(length(mu)) # Prior variance.
#' set.seed(123) 
#' pd <- MASS::mvrnorm(n = 10, mu = mu, Sigma = v) # 10 draws.
#' p.d <- list(matrix(pd[,1], ncol = 1), pd[,2:7])
#' CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), par.draws = p.d,
#' n.alts = 2, n.sets = 8, parallel = FALSE, alt.cte = c(0, 1))
#' # Or AB-efficient design
#' set.seed(123) 
#' CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), par.draws = p.d,
#' n.alts = 2, n.sets = 8, parallel = FALSE, alt.cte = c(0, 1), optim = "A")
#' 
#' # DB-efficient design with categorical and continuous factors
#' # 2 categorical attributes with 4 and 2 levels (effect coded) and 1 
#' # continuous attribute (= 5 parameters)
#' mu <- c(0.5, 0.8, 0.2, 0.4, 0.3) 
#' v <- diag(length(mu)) # Prior variance.
#' set.seed(123) 
#' pd <- MASS::mvrnorm(n = 3, mu = mu, Sigma = v) # 10 draws.
#' CEA(lvls = c(4, 2, 3), coding = c("E", "E", "C"), par.draws = pd,
#' c.lvls = list(c(2, 4, 6)), n.alts = 2, n.sets = 6, parallel = FALSE)
#' # The same can be done if A-optimality is chosen
#' set.seed(123)
#' CEA(lvls = c(4, 2, 3), coding = c("E", "E", "C"), par.draws = pd,
#' c.lvls = list(c(2, 4, 6)), n.alts = 2, n.sets = 6, parallel = FALSE, optim = "A")
#' 
#' # DB-efficient design with start design provided.  
#' # 3 Attributes with 3 levels, all dummy coded (= 6 parameters).
#' mu <- c(0.8, 0.2, -0.3, -0.2, 0.7, 0.4) 
#' v <- diag(length(mu)) # Prior variance.
#' sd <- list(example_design)
#' set.seed(123)
#' ps <- MASS::mvrnorm(n = 10, mu = mu, Sigma = v) # 10 draws.
#' CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), par.draws = ps,
#' n.alts = 2, n.sets = 8, parallel = FALSE, start.des = sd)
#' 
#' # DB-efficient design with partial profiles
#' # 3 Attributes, all dummy coded. = 6 parameters
#' mu <- c(1.2, 0.8, 0.2, -0.3, -1.2, 1.6) # Prior parameter vector
#' v <- diag(length(mu)) # Prior variance.
#' set.seed(123) 
#' pd <- MASS::mvrnorm(n = 10, mu = mu, Sigma = v) # 10 draws.
#' CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), par.draws = pd,
#' n.alts = 2, n.sets = 8, parallel = FALSE, alt.cte = c(0, 0), overlap = 1)
#' # The same function but asking for blocks (and no overlap)
#' set.seed(123)
#' CEA(lvls = c(3, 3, 3), coding = c("D", "D", "D"), par.draws = pd,
#' n.alts = 2, n.sets = 8, parallel = FALSE, alt.cte = c(0, 0), n.blocks = 2)
#' 
#' # AB-efficient design with constraints
#' # 2 dummy coded attributes, 1 continuous attribute and 1 effect coded
#' # attribute (with 4 levels). = 8 parameters
#' mu <- c(1.2, 0.8, 0.2, 0.5, -0.3, -1.2, 1, 1.6) # Prior parameter vector
#' v <- diag(length(mu)) # Prior variance.
#' set.seed(123) 
#' pd <- MASS::mvrnorm(n = 10, mu = mu, Sigma = v) # 10 draws.
#' constraints <- list("Alt2.Att1 = list(1A,1B)",
#'                     "if Alt1.Att3 = list(4) then Alt2.Att4 = list(4C, 4D)")
#' CEA(lvls = c(3, 3, 2, 4), coding = c("D", "D", "C", "E"), c.lvls = list(c(4,7)), par.draws = pd,
#' n.alts = 2, n.sets = 8, parallel = FALSE, alt.cte = c(0, 0), optim = "A", constraints = constraints)
#' 
#' 
#'}
#' @importFrom Rdpack reprompt
#' @export

CEA <- function(lvls, coding, c.lvls = NULL, n.sets, n.alts, par.draws, optim = "D", 
                alt.cte = NULL, no.choice = FALSE, start.des = NULL, 
                parallel = TRUE, max.iter = Inf, n.start = 12, overlap = NULL,
                n.blocks = 1, blocking.iter = 50, constraints = NULL) {
  
  ### Error handling for creating initial random design
  # Error lvls vector. At least 2 attributes
  if (length(lvls) < 2 || (!(is.numeric(lvls)))) {
    stop("lvls argument is incorrect.")
  }
  # Error correct coding types.
  codings.types <- c("E", "D", "C")
  if (!all(coding %in% codings.types) || (length(coding) != length(lvls))) {
    stop("coding argument is incorrect.")
  } 
  # Continuous attributes. 
  contins <-  which(coding == "C")
  n.contins <-  length(contins)
  # error continuous levels (Not a list)
  if (!is.null(c.lvls) && !is.list(c.lvls)) { 
    stop('c.lvls should be a list.')
  }
  # Error continuous specified and NULL.
  if (n.contins > 0 && is.null(c.lvls)) {
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
  
  ### Error handling for design specifications
  # If no alternative constant is given, create the variable as a vector of 0s
  if (is.null(alt.cte)) {
    alt.cte <- rep(0L, n.alts)
  }
  #init
  n.cte <- length(which(alt.cte == 1))
  
  # If only one draw is given, transform it to a matrix
  if (!is.list(par.draws)) {
    if (is.vector(par.draws)) {
      par.draws <- matrix(par.draws, nrow = 1)
    }
  }
  # Check if Optimality parameter is A or D only
  if (!(optim %in% c("A","a","D","d")) ) {
    stop("'optim' can only be 'A' or 'D'")
  }
  
  # Alternative constant errors
  if (length(alt.cte) != n.alts) {
    stop("'n.alts' does not match the 'alt.cte' vector")
  }
  if (!all(alt.cte %in% c(0, 1))) {
    stop("'alt.cte' should only contain zero or ones.")
  }
  
  # No choice errors
  if (!is.logical(no.choice)) {
    stop("'no.choice' should be TRUE or FALSE")
  }
  if (no.choice) {
    if (!isTRUE(all.equal(alt.cte[n.alts], 1))) {
      stop("if 'no.choice' is TRUE, the last alternative constant should equal 1.")
    }
    ncsek <- seq(n.alts, (n.sets * n.alts), n.alts)  # Rows in the design matrix
    # to be assigned as zeros in all attributes
  } else {
    ncsek <- NULL
  }
  
  # Handling par.draws with alternative specific constants.
  # This conditional is when there is only one alternative constant
  if (isTRUE(all.equal(n.cte, 1))) {
    if (!(is.list(par.draws))) {
      stop("par.draws should be a list")
    }
    if (!isTRUE(all.equal(length(par.draws), 2))) {
      stop("'par.draws' should contain two components")
    }
    # If only draws for the constant are given in a vector, transform it to a
    # matrix
    if (is.vector(par.draws[[1]])) {
      par.draws[[1]] <- matrix(par.draws[[1]], ncol = 1) 
    }
    if (!(all(unlist(lapply(par.draws, is.matrix))))) {
      stop("'par.draws' should contain two matrices")
    }
    dims <-  as.data.frame(lapply(par.draws, dim))
    if (!isTRUE(all.equal(dims[1, 1], dims[1, 2]))) { 
      stop("the number of rows in the components of 'par.draws' should be equal")
    }
    par.draws  <- do.call("cbind", par.draws) # Transform draws to a matrix
  } else {
    if (n.cte > 1.2) {
      if (!(is.list(par.draws))) { 
        stop("par.draws should be a list")
      } 
      if (!isTRUE(all.equal(length(par.draws), 2))) {
        stop("'par.draws' should contain two components")
      }
      if (!(all(unlist(lapply(par.draws, is.matrix))))) {
        stop("'par.draws' should contain two matrices")
      }
      if (!isTRUE(all.equal(ncol(par.draws[[1]]), n.cte))) {
        stop("the first component of 'par.draws' should contain the same number of columns as there are non zero elements in 'alt.cte'")
      }
      dims <-  as.data.frame(lapply(par.draws, dim))
      if (!isTRUE(all.equal(dims[1, 1], dims[1, 2]))) {
        stop("the number of rows in the components of 'par.draws' should be equal")
      }
      par.draws  <- do.call("cbind", par.draws) # Transform draws to a matrix
    }
  }
  
  # Error identifying model.
  if (n.sets < ncol(par.draws)) {
    stop("Model is unidentified. Increase the number of choice sets or decrease parameters to estimate.")
  }
  
  ### Create all levels for each attribute
  # Change into correct coding. 
  constraints_coding <- coding #meant to be used for constraints latter without modification
  coding <- dplyr::recode(coding, D = "contr.treatment", E = "contr.sum")
  # Create all combinations of attribute levels.
  levels.list <- lapply(X = as.list(lvls), function(x) (1:x))
  # Replace continuous.
  levels.list[contins] <- c.lvls
  # Transform to matrix
  levels.list[contins] <- lapply(X = levels.list[contins], as.matrix)
  # Transform categorical attributes
  categ <-  which(coding %in% c("contr.treatment", "contr.sum"))
  n.categ <- length(categ)
  if (n.categ > 0) {
    levels.list[categ] <- lapply(X = levels.list[categ], factor)
    # Apply coding
    if (n.contins > 0) {
      for (i in 1:length(lvls)) {
        if (!(i %in% contins)) {
          stats::contrasts(levels.list[[i]]) <- coding[i]
        }
      }
    } else {
      for (i in 1:length(lvls)) {
        stats::contrasts(levels.list[[i]]) <- coding[i]
      }
    }
    # Compute all possible values for each categorical attribute
    levels.list[categ] <- lapply(X = levels.list[categ], stats::contrasts)
  }
  # Set colnames for the design matrix
  c.nam <- list()
  for (i in 1:length(lvls)) {
    if (coding[i] == "contr.treatment") {
      c.nam[[i]] <- paste("Var", i, 2:lvls[i], sep = "")
    } else {
      if (coding[i] == "contr.sum") {
        c.nam[[i]] <- paste("Var", i, 1:(lvls[i] - 1), sep = "")
      } else {
        c.nam[[i]] <- paste("Var", i, sep = "")
      }
    }
  }
  c.nam <- unlist(c.nam)
  
  # Create alternative specific design.
  cte.des <- Altspec(alt.cte = alt.cte, n.sets = n.sets)
  
  # Count the number of columns of the design matrix without constants
  ncol.des.noconst <- sum(unlist(lapply(levels.list,ncol)))
  #ncol.des.noconst <- sum(unlist(lapply(levels.list,function(x){
  #  if (is.matrix(x)) { ncol(x) }
  #  else{if (is.numeric(x)) { 1 }}
  #})))
  
  if (!identical(as.integer(ncol(par.draws)), 
                 as.integer(n.cte + ncol.des.noconst))) {
    stop("The sum of the number of columns in the components of 'par.draws' should equal the number of columns of design matrix (including alternative specific constants)")
  }
  
  ### Random initial design.
  num_attributes <- as.numeric(length(lvls))  #calculate the number of attributes
  
  if (is.null(start.des)) {
    #create start designs
    nr.starts <- n.start
    start.des <- vector(mode = 'list', length = nr.starts)
    okstart <- FALSE
    while (okstart == FALSE) {
      for (i in 1:nr.starts) {
        
        # r is to know which levels to take in each attribute
        r <- NULL
        start <- NULL
        
        if (!is.null(overlap) && overlap > 0) {
          # Sample a random variable for overlapping in the starting des
          overlap_var <- sample(num_attributes, overlap)
        }
        
        for (j in 1:length(lvls)) {
          
          if(!is.null(overlap) && j %in% overlap_var ){
            #with overlapping
            r <- round(unlist(replicate(n.sets,rep( stats::runif(1, 1, lvls[j]),n.alts ), simplify = F))) ##QI
          } else {
            #without overlapping
            r <- round(stats::runif((n.sets * n.alts), 1, lvls[j]))
          }
          start <- cbind(start, levels.list[[j]][r,])
        }
        colnames(start) <- c.nam
        start.des[[i]] <- cbind(cte.des, start)
        if (no.choice) {
          start.des[[i]][ncsek, (ncol(cte.des) + 1):(ncol(cte.des) + ncol.des.noconst)] <- c(rep(0, ncol.des.noconst))
        }
        
        #constraints####
        if (!is.null(constraints)) {
          sets_indeces <- ((c(1:nrow(start.des[[1]])) - 1) %/% n.alts) + 1 #get the choice set number
          sets_row_start <- unique((sets_indeces - 1) * n.alts + 1) #get the start row in the choice sets
          sets_row_end <- unique(ifelse(rep(no.choice,length(sets_indeces)), sets_indeces * n.alts - 1, sets_indeces * n.alts)) #get the end row in the choice sets
          
          #first translate the constraints on coded variables
          constraints <- lapply(constraints, function(x) {translate_coded_constraints(x, coding = constraints_coding , lvls = lvls, cand.set = NULL)} )
          
          #extract the attributes in the constraints to alter them
          constrained_atts <- unique(sort(unlist(lapply(constraints, ext_constrained_atts))))
          
          #check that the no.choice alternative is included in any of the constraints
          if (no.choice) {
            constrained_alts <- unique(sort(unlist(lapply(constraints, ext_constrained_alts))))
            if (n.alts %in% constrained_alts) {
              stop(paste("Alternative", n.alts, "has the no.choice alternative and should not be included in the constraints."))
            }
          }
            
          #then we do the checking for every choice set  
            for (s in 1:n.sets) {
              check <- all(unlist(lapply(constraints, function(x) {check_constraints(start.des[[i]],x,sets_row_start[s],n.cte)})))
              counter <- 0
              while (check == FALSE) {
                #loop over the possible variables in the constrained attributes until the constraints are satisfied in the set s
                for (j in constrained_atts) { #3
                  cols.j <- which(categorize_variables(colnames(start.des[[i]][,(n.cte+1):ncol(start.des[[i]])] )) %in% j) + n.cte
                  if (!is.null(overlap) && j %in% overlap_var) {
                    #with overlap
                    r <- round(unlist(stats::runif(1, 1, lvls[j]),n.alts )) ##QI
                    for (alt in sets_row_start[s]:sets_row_end[s]) {
                      start.des[[i]][alt, cols.j ] <- levels.list[[j]][r,]   
                    }
                  } else {
                    #without overlap
                    draws <- ifelse(no.choice, n.alts-1, n.alts)
                    r <- round(stats::runif(draws, 1, lvls[j]))
                    start.des[[i]][c(sets_row_start[s]:sets_row_end[s]), cols.j ] <- levels.list[[j]][r,]
                  }
                } #3
              
                #check the constraint again
                check <- all(unlist(lapply(constraints, function(x) {check_constraints(start.des[[i]],x,sets_row_start[s],n.cte)})))
                counter <- counter + 1
                #we leave a message advising to stop the run, if the constraints are too tight.
                if (counter == 1000) {
                  message("Exceeded 1000 trials to satisfy the constraints. Consider stopping the run and modifying them.")
                }
              }
            } #end of loop for set s
        } #end of loop for constraints
      }
      # Compute D or A-optimality for each design and each draw
      if(optim %in% c("D","d") ) {
        d.start <- lapply(start.des, StartDB, par.draws, n.alts) #DB-error
      } else {
        d.start <- lapply(start.des, StartAB, par.draws, n.alts) #AB-error
      }
      # If the DB or AB-optimality of any starting design is finite, continue
      if (any(is.finite(unlist(lapply(d.start, mean, na.rm = TRUE))))) {
        okstart <- TRUE
      } 
    }
  } else {
    if (!is.list(start.des)) {
      stop("'start.des' should be a list")
    }
    if (!(all(unlist(lapply(start.des, is.matrix))))) {
      stop("'start.des' should contain matrices as components")
    }
    # Save the dimension of each starting design
    dimstart <- as.matrix(lapply(start.des, dim)) 
    # Save the number of random starts given
    nr.starts <- length(dimstart)
    if (nr.starts > 1.5) {
      if (!isTRUE(all.equal(length(unique(unlist(dimstart))), 2))) {
        stop("start designs have different dimensions")
      }
    }
    if (!isTRUE(all.equal(n.alts * n.sets, unique(unlist(dimstart))[1]))) {
      stop("number of rows of start design(s) does not match with 'n.alts' * 'n.sets'")
    }
    if (!isTRUE(all.equal(as.integer(n.cte + ncol.des.noconst), 
                          unique(unlist(dimstart))[2]))) {
      stop("number of columns of start design(s) does not match with the number of columns in the design matrix")
    }
    # Compute D or A-optimality for each design and each draw
    if(optim %in% c("D","d") ) {
      d.start <- lapply(start.des, StartDB, par.draws, n.alts) #DB-error
    } else {
      d.start <- lapply(start.des, StartAB, par.draws, n.alts) #AB-error
    }
    
    if (!any(is.finite(unlist(lapply(d.start, mean, na.rm = TRUE))))) {
      stop("One or more of the provided start designs resulted in an invalid db-error.")
    }
  }
  
  ### Improving the initial design
  if (parallel) {
    no_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(no_cores)
    parallel::clusterExport(cl, c("n.sets", "par.draws", "n.alts", "n.cte", 
                                  "alt.cte", "no.choice", "max.iter","ncsek",
                                  "overlap", "num_attributes", "optim", "constraints"),
                            envir = environment())
    deslist <- parallel::parLapply(cl, start.des, CEAcore_ucpp, par.draws, 
                                   levels.list, n.alts, n.sets, n.cte, alt.cte, 
                                   no.choice, max.iter, ncsek, 
                                   overlap, num_attributes, optim, constraints)
    parallel::stopCluster(cl)
  } else {
    deslist <- lapply(start.des, CEAcore_ucpp, par.draws, levels.list, 
                      n.alts, n.sets, n.cte, alt.cte, no.choice, 
                      max.iter = max.iter, ncsek, overlap, num_attributes, optim, constraints) 
  }                                 

  #measures/stats ####
  #1-errors
  ab.error <- unlist(lapply(deslist, function(x) {
    mean(apply(par.draws, 1, Aerr_ucpp, des = x$design, n.alts = n.alts), na.rm = TRUE)}))
  names(ab.error) <- paste0("des", seq_along(start.des))
  
  db.error <- unlist(lapply(deslist, function(x) {
    mean(apply(par.draws, 1, Derr_ucpp, des <- x$design, n.alts = n.alts), na.rm = TRUE)}))
  names(db.error) <- paste0("des", seq_along(start.des))
  
  #2- SD of the parameter estimates
  SD <- lapply(deslist, function(x) {
    rowMeans( apply(par.draws, 1, SD, des = x$design, n.alts = n.alts), na.rm = TRUE)
  })
  SD.mat <- do.call(rbind, SD)
  rownames(SD.mat) <- paste0("des", seq_along(start.des))
  
  #3-orthogonality
  orthogonality <- unlist(lapply(deslist, function(x) {
    mean(apply(par.draws, 1, Orthogonality, des = x$design, n.alts = n.alts), na.rm = TRUE)}))
    # Orthogonality(x$design, n.alts)}))
  names(orthogonality) <- paste0("des", seq_along(start.des))
  
  #4-level overlaps
  if (no.choice) {
    lvl.overlap <- lapply(deslist, function(x) {lvl.overlap(x$design[-ncsek,(ncol(cte.des) + 1):(ncol(cte.des) + ncol.des.noconst)], n.alts-1)})   
  } else {
    lvl.overlap <- lapply(deslist, function(x) {lvl.overlap(x$design[,(ncol(cte.des) + 1):(ncol(cte.des) + ncol.des.noconst)], n.alts)})
  }
  lvl.overlap_final <- lvl.overlap[[1]][, 1, drop=FALSE]
  for (i in seq_along(start.des)) {
    lvl.overlap_final <- cbind(lvl.overlap_final,lvl.overlap[[i]][, 2])
  }
  colnames(lvl.overlap_final)[2:(length(start.des)+1)] <- paste0("des", 1:length(start.des))
  
  #5-level frequencies
  if (no.choice) {
    lvl.freq <- lapply(deslist, function(x) {lvl.freq(x$design[-ncsek,(ncol(cte.des) + 1):(ncol(cte.des) + ncol.des.noconst)]) } )   
  } else {
    lvl.freq <- lapply(deslist, function(x) {lvl.freq(x$design[,(ncol(cte.des) + 1):(ncol(cte.des) + ncol.des.noconst)]) } )   
  }
  lvl.freq_final <- lvl.freq[[1]]
  n <- length(lvl.freq_final)
  for (a in 1:n) {
    if (length(start.des)>1) {
      for (i in 2:length(start.des)) {
        temp <- dplyr::bind_rows(lvl.freq_final[[a]], (lvl.freq[[i]][[a]]))
        temp[is.na(temp)] <- 0 #replace NAs with zeros
        lvl.freq_final[[a]] <- temp   #cbind(lvl.overlap_final,lvl.overlap[[i]][, 2])
        }
    }
    rownames(lvl.freq_final[[a]]) <- paste0("des",1:length(start.des))
  }
  #
  stat <- list(AB.error = ab.error, DB.error = db.error, SD = SD.mat, level.count = lvl.freq_final,
              level.overlap = lvl.overlap_final, Orthogonality = orthogonality)
  
  #best design
  bestdes_ind <- ifelse(optim == "D", which.min(db.error), which.min(ab.error))
  bestdes <- deslist[[bestdes_ind]]
  bestdes_stat <- list(
    AB.error = stat$AB.error[bestdes_ind],
    DB.error = stat$DB.error[bestdes_ind],
    SD = stat$SD[bestdes_ind,],
    level.count = lapply(stat$level.count, function(x) x[bestdes_ind,]),
    level.overlap = stat$level.overlap[,c(1,bestdes_ind+1)],
    Orthogonality = stat$Orthogonality[bestdes_ind]
  )
  #add blocking (if any) and AllDesigns
  if (n.blocks > 1) {
    bestdes_blocks <- Blocks(bestdes$design, n.blocks, n.alts, blocking.iter, 
                            no.choice, alt.cte) 
    output_design <- list(BestDesign = c(bestdes[-length(bestdes)], bestdes_stat, bestdes_blocks), 
                          AllDesigns = c(lapply(deslist, function(x) x[-length(x)]), stat)) #lapply is used to remove $error from each list to avoid duplication in stat
  } else {
      output_design <- list(BestDesign = c(bestdes[-length(bestdes)], bestdes_stat), 
                            AllDesigns = c(lapply(deslist, function(x) x[-length(x)]), stat)) #lapply is used to remove $error from each list to avoid duplication in stat
  }
  
  ##
  class(output_design) <- 'design_list'
  return(output_design)
}

# Core of Coordinate Exchange
CEAcore_ucpp <- function(des, par.draws, levels.list, n.alts, n.sets, n.cte, 
                          alt.cte, no.choice, max.iter, ncsek, overlap, num_attributes, optim, constraints) {
  converge <- FALSE # Boolean for convergence
  change <- FALSE # Boolean for a change in the attribute
  it <- 1 # Indicator of number of iterations
  n.samples <- nrow(par.draws)  # Number of samples from distribution of betas
  n.par <- ncol(des) # Number of parameters

  # Beginning of algorithm
  while (!converge & it <= max.iter) {
    # Compute DB or AB-optimality for initial design
    if(optim %in% c("D","d") ) {
      db.start <- mean(apply(par.draws, 1, Derr_ucpp, des = des, 
                             n.alts = n.alts), na.rm = TRUE) #DB-error
    } else {
      db.start <- mean(apply(par.draws, 1, Aerr_ucpp, des = des, 
                             n.alts = n.alts), na.rm = TRUE) #AB-error
    }
    
    it <- it + 1 # Increase iteration
    # Save design before iteration.
    iter.des <- des
    # For every row in the design.
    sek <- 1:nrow(des)
    if (no.choice) {
      sek <- sek[-ncsek] # If no.choice is given, then it is not improved
    }
    # Loop for each row of the initial design
    for (i in sek) {
      
      #Initialize the set index and the starting and ending row (to be used for overlapping/constraints)
      set_index <- ((i - 1) %/% n.alts) + 1
      set_row_start <- (set_index - 1) * n.alts + 1
      set_row_end <- ifelse(no.choice, set_index * n.alts - 1, set_index * n.alts)
      
      if (!is.null(overlap) && overlap > 0) {
        #get all combinations of the variables to overlap
        attribute_combinations <- utils::combn(num_attributes, overlap, simplify = FALSE)
        
        for (combo in attribute_combinations) {
          
          overlap_var <- combo
          
          # Loop for each attribute in the row
          for (j in sample(1:length(levels.list))) { ##added sample to shuffle the order of variables so as to eliminate any ordering effect
            #%%%%%%%
            # This can be improved by removing mods object and just replace the 
            # initial design with the new attributes to compute the db.error and 
            # then if one of those is better than the initial db error. Replace 
            # finally the initial design (as it is in modfed)
            # Initialize mods object to track changes of the design in each attribute
            mods <- vector("list", nrow(levels.list[[j]])) 
            mods <- lapply(mods,function(x){des})
            # Initialize db object for db (or ab) error of each level
            db <- numeric(nrow(levels.list[[j]]))
            # Indicator of columns modified
            ncol.ok <- sum(unlist(lapply(levels.list[(1:j - 1)],ncol)))
            # Indicator of columns to change
            ncol.mod <- ifelse( j == 1, (n.cte + j), 
                                (n.cte + 1 + ncol.ok))
            # ncol.left tracks which columns are left to be improved in the row
            ncol.left <- sum(unlist(lapply(levels.list[-(1:j)],ncol)))
            # Loop for each level in the attribute
            for (k in 1:nrow(levels.list[[j]])) {
              # Update design with new attribute
              
              ## QI
              if (j %in% overlap_var) {
                # Apply the level change to all alternatives in the same choice set
                for (alt in set_row_start:set_row_end) {
                  mods[[k]][alt, ncol.mod:(n.par - ncol.left)] <- levels.list[[j]][k,]
                }
                
              } else {
                mods[[k]][i, ncol.mod:(n.par - ncol.left) ] <- levels.list[[j]][k,]
              }
              
              # Check if the constraints are met, and skip to the next row if negative
              if (!is.null(constraints)) {
                print("checking")
                print(unlist(constraints))
                check  <- all(unlist(lapply((constraints), function(x) {check_constraints(mods[[k]],x,set_row_start,n.cte)})))
                print(paste("checked", check))
                if (!check) {
                  db[k] <- NA
                  next
                } 
              }
              
              # Calculate D or A-optimality for each draw
              if (optim %in% c("D","d")) {
                d.errors <- apply(par.draws, 1, Derr_ucpp, des = mods[[k]], 
                                  n.alts = n.alts) #DB-error
              } else {
                d.errors <- apply(par.draws, 1, Aerr_ucpp, des = mods[[k]], 
                                  n.alts = n.alts) #AB-error
              }
              
              # Compute DB or AB-optimality 
              db[k] <- mean(d.errors, na.rm = TRUE)
            }
            
            pr <- which.min(db)
            db <- suppressWarnings(min(db, na.rm = TRUE)) #if all errors are equal to NaN, db will be Inf. We suppress this warning here.
            
            # Change if lower db (or ab) error.
            if (!is.na(db) && !is.na(db.start)) {
              if (db < db.start) {
                
                ###check that the overlapping attributes number is still respected with the new design
                choice_set  <- mods[[pr]][set_row_start:set_row_end,]
                attributes_index <- categorize_variables(colnames(des[,(n.cte+1):(ncol(des))]))
                overlap_count <- 0
                for(a in unique(attributes_index) ) {
                  index <- which( a == as.numeric(attributes_index) ) + n.cte
                  comparisons <- rep(TRUE, nrow(choice_set))
                  for (r in 1:nrow(choice_set) ) {
                    comparisons[r] <- all(choice_set[1, index] == choice_set[r, index])
                  }
                  overlapping_attrib <- sum(all(comparisons == TRUE))
                  overlap_count <- overlap_count + overlapping_attrib
                }
                ####
                if (overlap_count >= overlap) {
                  des <- mods[[pr]]
                  db.start <- db
                }
              }
            }
          } # End loop for each attribute in the row
        } # End loop for each combo
      } else {
        ## if no Overlap then proceed as usual
        for (j in 1:length(levels.list)) {
          #%%%%%%%
          # This can be improved by removing mods object and just replace the 
          # initial design with the new attributes to compute the db.error and 
          # then if one of those is better than the initial db error. Replace 
          # finally the initial design (as it is in modfed)
          # Initialize mods object to track changes of the design in each attribute
          mods <- vector("list", nrow(levels.list[[j]])) 
          mods <- lapply(mods,function(x){des})
          # Initialize db object for db (or ab) error of each level
          db <- numeric(nrow(levels.list[[j]]))
          # Indicator of columns modified
          ncol.ok <- sum(unlist(lapply(levels.list[(1:j - 1)],ncol)))
          # Indicator of columns to change
          ncol.mod <- ifelse( j == 1, (n.cte + j), 
                              (n.cte + 1 + ncol.ok))
          # ncol.left tracks which columns are left to be improved in the row
          ncol.left <- sum(unlist(lapply(levels.list[-(1:j)],ncol)))
          # Loop for each level in the attribute
          for (k in 1:nrow(levels.list[[j]])) {
            # Update design with new attribute
            mods[[k]][i, ncol.mod:(n.par - ncol.left) ] <- levels.list[[j]][k,]
            
            #check for constraints if any
            if (!is.null(constraints)) {
              check  <- all(unlist(lapply((constraints), function(x) {check_constraints(mods[[k]],x,set_row_start,n.cte)})))
              if (!check) {
                db[k] <- NA
                next
              } 
            }
            
            # Calculate D or A-optimality for each draw
            if (optim %in% c("D","d")) {
              d.errors <- apply(par.draws, 1, Derr_ucpp, des = mods[[k]], 
                                n.alts = n.alts) #DB-error
            } else {
              #A-errors
              d.errors <- apply(par.draws, 1, Aerr_ucpp, des = mods[[k]], 
                                n.alts = n.alts) #AB-error
            }
            
            # Compute DB or AB-optimality 
            db[k] <- mean(d.errors, na.rm = TRUE)
          }
          pr <- which.min(db)
          db <- suppressWarnings(min(db, na.rm = TRUE)) 
          # Change if lower db error.
          if (!is.na(db) && !is.na(db.start)) {
            if (db < db.start) {
              des <- mods[[pr]]
              db.start <- db
            }
          }
        } # End loop for each attribute in the row
      }
    } # End loop for the whole row
    converge <- isTRUE(all.equal(des, iter.des)) # Convergence if no profile is swapped this iteration.
  } # End while (after convergence)
  
  # calculate percentage NA values.
  na.percentage <- 0
  if(optim %in% c("D","d") ) {
    d.errors <- apply(par.draws, 1, Derr_ucpp, des = des,  n.alts = n.alts)
  } else {
    d.errors <- apply(par.draws, 1, Aerr_ucpp, des = des,  n.alts = n.alts)
  }
  if (any(is.na(d.errors))) {
    na.percentage <- scales::percent(sum(is.na(d.errors)) / n.samples)
  } 
  # Utility balance.
  # ub <- apply(par.draws, 1, Utbal, des = des,  n.alts = n.alts)
  # ub2 <- apply(par.draws, 1, InfoDes2, des = des,  n.alts = n.alts, utbal = TRUE)
  # ub_ucpp <- apply(par.draws, 1, InfoDes_cpp, des = des,  n_alts = n.alts, 
  #                  utbal = TRUE)
  
  # Utility balance using c++ function
  ub <- apply(par.draws, 1, InfoDes_cpp, des = des,  n_alts = n.alts, 
              utbal = TRUE)
  pmat <- matrix(rowMeans(ub), ncol = n.alts, byrow = TRUE)
  rownames(pmat) <- paste("set", 1:n.sets, sep = "")
  colnames(pmat) <- paste(paste("Pr(", paste("alt", 1:n.alts, sep = ""), 
                                sep = ""), ")", sep = "")
  if (no.choice) {
    colnames(pmat)[n.alts] <- "Pr(no choice)"
  }
  # Rownames design. 
  des.names <- Rcnames(n.sets = n.sets, n.alts = n.alts, alt.cte = alt.cte, 
                       no.choice = no.choice)
  rownames(des) <- des.names[[1]]
  # Colnames alternative specific constants. 
  if (n.cte != 0 && !is.null(colnames(des))) {
    colnames(des)[1:n.cte] <- des.names[[2]]
  }
  # Return design, D(B)error, percentage NA's, utility balance. 
  return(list("design" = des, "optimality" =  optim, "inf.error" = na.percentage,
              "probs" = pmat, "error" =  db.start))
}


#' Sequential Coordinate Exchange algorithm for MNL model.
#' 
#' Selects the choice set that minimizes the DB-error when added to an initial 
#' design, given (updated) parameter values.
#' 
#' This algorithm is ideally used in an adaptive context. The algorithm will 
#' select the next DB-efficient choice set given parameter values and possible
#' previously generated choice sets. In an adaptive context these parameter
#' values are updated after each observed response.
#' 
#' Previously generated choice sets, which together form an initial design, can
#' be provided in \code{des}. When no design is provided, the algorithm will
#' select the most efficient choice set based on the fisher information of the
#' prior covariance matrix \code{prior.covar}.
#' 
#' If \code{alt.cte = NULL}, \code{par.draws} should be a matrix in which each 
#' row is a sample from the multivariate parameter distribution. In case that 
#' \code{alt.cte} is not \code{NULL}, a list containing two matrices should be 
#' provided to \code{par.draws}. The first matrix containing the parameter draws
#' for the alternative specific parameters. The second matrix containing the
#' draws for the rest of the parameters.
#' 
#' The list of potential choice sets is created by selecting randomly a level for
#' each attribute in an alternative/profile. \code{n.cs} controls the number of
#' potential choice sets to consider. The default is \code{
#' NULL}, which means that the number of possible choice sets is the product of
#' attribute levels considered in the experiment. For instance, an experiment 
#' with 3 attribute and 3 levels each will consider 3^3 = 27 possible choice sets. 
#' 
#' The \code{weights} argument can be used when the \code{par.draws} have 
#' weights. This is for example the case when parameter values are updated using
#' \code{\link{ImpsampMNL}}.
#' 
#' When \code{parallel} is \code{TRUE}, \code{\link[parallel]{detectCores}} will
#' be used to decide upon the number of available cores. That number minus 1 
#' cores will be used to search for the optimal choice set. For small problems 
#' (6 parameters), \code{parallel = TRUE} can be slower. For larger problems the
#' computation time will decrease significantly.
#' 
#' *Note:* this function is faster than \code{\link[idefix]{SeqMOD}}, but 
#' the output is not as stable. This happens because this function 
#' makes a random search to get the choice set, whereas 
#' \code{\link[idefix]{SeqMOD}} makes an exhaustive search.
#' @inheritParams CEA
#' @param par.draws A matrix or a list, depending on \code{alt.cte}. 
#' @param des A design matrix in which each row is a profile. If alternative 
#'   specific constants are present, those should be included as the first 
#'   column(s) of the design. Can be generated with \code{\link{Modfed}} or
#'   \code{\link{CEA}}
#' @param n.cs An integer indicating the number of possible random choice sets to 
#' consider in the search for the next best choice set possible. The default is
#'  \code{NULL}. 
#' @param prior.covar Covariance matrix of the prior distribution.
#' @param weights A vector containing the weights of the draws. Default is 
#'   \code{NULL}. See also \code{\link{ImpsampMNL}}.
#' @param parallel Logical value indicating whether computations should be done 
#'   over multiple cores.
#' @param no.choice An integer indicating the no choice alternative. The default
#'   is \code{NULL}.
#' @param reduce Logical value indicating whether the candidate set should be 
#'   reduced or not.
#' @return \item{set}{A matrix representing a DB efficient choice set.} 
#'   \item{error}{A numeric value indicating the DB-error of the whole 
#'   design.}
#' @importFrom Rdpack reprompt
#' @references \insertRef{idefix}{idefix}
#' @references \insertRef{ju}{idefix}
#' @references \insertRef{cea}{idefix}
#' @references \insertRef{cea_discrete}{idefix}
#' @examples 
#' # DB efficient choice set, given a design and parameter draws. 
#' # 3 attributes with 3 levels each
#' m <- c(0.3, 0.2, -0.3, -0.2, 1.1, 2.4) # mean (total = 6 parameters).
#' pc <- diag(length(m)) # covariance matrix
#' set.seed(123)
#' sample <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc)
#' # Initial design.
#' des <- example_design
#' # Efficient choice set to add.
#' SeqCEA(des = des, lvls = c(3, 3, 3), coding = c("D", "D", "D"), n.alts = 2,
#'        par.draws = sample, prior.covar = pc, parallel = FALSE)
#' 
#' # DB efficient choice set, given parameter draws. 
#' # with alternative specific constants 
#' des <- example_design2
#' ac <- c(1, 1, 0) # Alternative specific constants.
#' m <- c(0.3, 0.2, -0.3, -0.2, 1.1, 2.4, 1.8, 1.2) # mean
#' pc <- diag(length(m)) # covariance matrix
#' pos <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc)
#' sample <- list(pos[ , 1:2], pos[ , 3:8])
#' # Efficient choice set.
#' SeqCEA(des = des, lvls = c(3, 3, 3), coding = c("D", "D", "D"), n.alts = 3, 
#'       par.draws = sample, alt.cte = ac, prior.covar = pc, parallel = FALSE)
#' @export
SeqCEA <- function(des = NULL, lvls, coding, c.lvls = NULL, n.alts, par.draws, 
                   prior.covar, alt.cte = NULL, no.choice = NULL,
                   weights = NULL, parallel = TRUE, reduce = TRUE, n.cs = NULL) {
  # Error handling initial design
  if (is.null(des)) {
    n.sets <- 1L
  } else { 
    if (!is.matrix(des)) {
      stop("'des' should be a matrix or NULL")
    }
    if (!isTRUE(nrow(des) %% n.alts == 0)) {
      stop("'n.alts' does not seem to match with the number of rows in 'des'")
    }
    n.sets <- nrow(des) / n.alts
  }
  
  ### Error handling for design specifications
  # No choice errors
  if (!is.null(no.choice)) {
    if (!is.wholenumber(no.choice)) {
      stop("'no.choice' should be an integer or NULL")
    }
    if (any(isTRUE(no.choice > (n.alts + 0.2)), isTRUE(no.choice < 0.2))) {
      stop("'no.choice' does not indicate one of the alternatives")
    }
    if (is.null(alt.cte)) {
      stop("if there is a no choice alternative, 'alt.cte' should be specified")
    }
    if (!isTRUE(all.equal(alt.cte[no.choice], 1))) {
      stop("the no choice alternative should correspond with a 1 in 'alt.cte'")
    }
  }
  
  # Alternative constant errors
  if (!is.null(alt.cte)) {
    if (length(alt.cte) != n.alts) {
      stop("'n.alts' does not match the 'alt.cte' vector")
    }
    if (!all(alt.cte %in% c(0, 1))) {
      stop("'alt.cte' should only contain zero or ones.")
    }
    # alternative specific constants
    n.cte <- length(which(alt.cte == 1L))
    if (isTRUE(all.equal(n.cte, 0L))) {
      alt.cte <- NULL
      cte.des <- NULL
    }
    
    # Handling errors when there are alternative constants
    if (!is.null(alt.cte)) {
      if (n.cte > 0.2) {
        if (!is.list(par.draws)) {
          stop("'par.draws' should be a list when 'alt.cte' is not NULL")
        }
        if (!isTRUE(all.equal(length(par.draws), 2))) {
          stop("'par.draws' should contain two components")
        }
        # If there is only one specific constant and is a vector, then it is
        # transformed to a matrix
        if (isTRUE(all.equal(n.cte, 1))) {
          if (is.vector(par.draws[[1]])) {
            par.draws[[1]] <- matrix(par.draws[[1]], ncol = 1)
          }
        }
        if (!(all(unlist(lapply(par.draws, is.matrix))))) {
          stop("'par.draws' should contain two matrices")
        }
        if (!isTRUE(all.equal(ncol(par.draws[[1]]), n.cte))) {
          stop("the first component of 'par.draws' should contain the same number 
               of columns as there are non zero elements in 'alt.cte'")
        }
        dims <-  as.data.frame(lapply(par.draws, dim))
        if (!isTRUE(all.equal(dims[1, 1], dims[1, 2]))) { 
          stop("the number of rows in the components of 'par.draws' should be equal")
        }
        par.draws  <- do.call("cbind", par.draws) # Transform par.draws to a matrix
      }
      # Create alternative specific design.
      cte.des <- Altspec(alt.cte = alt.cte, n.sets = n.sets)
      cte.set <- matrix(cte.des[1:n.alts, ], ncol = n.cte, byrow = FALSE)
    } 
  } else {
    cte.des <- NULL
    n.cte <- 0
    # if no alternative constants 
    if (!is.matrix(par.draws)) {
      stop("'par.draws'should be a matrix when 'alt.cte' = NULL")
    }
  }
  
  # Weights errors
  n.par <- ncol(par.draws)
  if (!is.null(weights)) {
    if (!isTRUE(all.equal(length(weights), nrow(par.draws)))) {
      stop("length of 'weights' does not match number total number of rows in 'par.draws'")
    }
  } else {
    weights <- rep(1L, nrow(par.draws))
  }
  
  ### Create all levels for each attribute
  # Continuous attributes. 
  contins <-  which(coding == "C")
  n.contins <-  length(contins)
  # Change into correct coding. 
  coding <- dplyr::recode(coding, D = "contr.treatment", E = "contr.sum")
  # Create all combinations of attribute levels.
  levels.list <- lapply(X = as.list(lvls), function(x) (1:x))
  # Replace continuous.
  levels.list[contins] <- c.lvls
  # Transform to matrix
  levels.list[contins] <- lapply(X = levels.list[contins], as.matrix)
  # Transform categorical attributes
  categ <-  which(coding %in% c("contr.treatment", "contr.sum"))
  n.categ <- length(categ)
  if (n.categ > 0) {
    levels.list[categ] <- lapply(X = levels.list[categ], factor)
    # Apply coding
    if (n.contins > 0) {
      for (i in 1:length(lvls)) {
        if (!(i %in% contins)) {
          stats::contrasts(levels.list[[i]]) <- coding[i]
        }
      }
    }else {
      for (i in 1:length(lvls)) {
        stats::contrasts(levels.list[[i]]) <- coding[i]
      }
    }
    # Compute all possible values for each categorical attribute
    levels.list[categ] <- lapply(X = levels.list[categ], stats::contrasts)
  }
  # Set colnames for the design matrix
  c.nam = list()
  for (i in 1:length(lvls)) {
    if (coding[i] == "contr.treatment") {
      c.nam[[i]] <- paste("Var", i, 2:lvls[i], sep = "")
    } else {
      if (coding[i] == "contr.sum") {
        c.nam[[i]] <- paste("Var", i, 1:(lvls[i] - 1), sep = "")
      } else {
        c.nam[[i]] <- paste("Var", i, sep = "")
      }
    }
  }
  c.nam = unlist(c.nam)
  
  # Count the number of colums of the design matrix without constants
  ncol.des.noconst <- sum(unlist(lapply(levels.list,ncol)))
  
  # if (!identical(as.integer(n.par), as.integer(n.cte + ncol.des.noconst))) {
  #   stop("The sum of the number of columns in the components of 'par.draws' should equal the number of columns of design matrix (including alternative specific constants)")
  # }
  
  if (!identical(as.integer(ncol(prior.covar)), 
                 as.integer(n.cte + ncol.des.noconst))) {
    stop("number of columns of 'prior.covar' does not equal the number of columns of design matrix (including alternative specific constants)")
  }
  
  ## When a design is supplied
  if (!is.null(des)) {
    # Error par.draws
    if (!isTRUE(all.equal(ncol(des), n.par))) {
      stop("number of columns in 'par.draws' does not match the number of columns in 'des'")
    }
    # Error dimension of initial design and design matrix
    if (!identical(as.integer(ncol(des)), as.integer(n.cte + ncol.des.noconst))) {
      stop("number of columns in 'des' does not match the number of columns of design matrix (including alternative specific constants)")
    }
    # Starting and initializing values.
    i.cov <- solve(prior.covar)
    d.start <- apply(par.draws, 1, DerrC_ucpp, des = des, n.alts = n.alts, 
                     i.cov = i.cov)
    db.start <- mean(d.start, na.rm = TRUE)
    full.comb <- Newsets_ucpp(levels.list = levels.list, n.alts = n.alts, 
                              no.choice = no.choice, reduce = reduce, 
                              c.names = c.nam, n.cs = n.cs)
    # Adding alternative constants
    if (!is.null(cte.des)) {
      full.comb <- lapply(full.comb, function(x) cbind(cte.set, x))
    }
    
    # Adding these new choice sets to the initial design
    full.des <- lapply(full.comb, function(x) rbind(des, x))
    
    # Select the next best choice set
    if (parallel) {
      no_cores <- parallel::detectCores() - 1L
      cl <- parallel::makeCluster(no_cores)
      db.errors <- parallel::parLapply(cl, full.des, DBerrS.P_ucpp, par.draws, 
                                       n.alts, i.cov, weights)
      parallel::stopCluster(cl)
    } else {
      # For each potential set, select best. 
      db.errors <- lapply(full.des, DBerrS.P_ucpp, par.draws, n.alts, i.cov,
                          weights)
    }
    dbs <- unlist(db.errors, use.names = FALSE)
    set <- full.comb[[which.min(dbs)]]
    row.names(set) <- NULL
    colnames(set) <- colnames(des)
    db <- min(dbs)
    #return best set and db error design.
    return(list("set" = set, "error" = db))
  } # End if when a design is supplied
  else {
    # Starting and initializing values.
    i.cov <- solve(prior.covar)
    full.comb <- Newsets_ucpp(levels.list = levels.list, n.alts = n.alts, 
                              no.choice = no.choice, reduce = reduce, 
                              c.names = c.nam, n.cs = n.cs)
    # Adding alternative constants
    if (!is.null(cte.des)) {
      full.comb <- lapply(full.comb, function(x) cbind(cte.set, x))
    }
    # Select the next best choice set
    if (parallel) {
      no_cores <- parallel::detectCores() - 1L
      cl <- parallel::makeCluster(no_cores)
      db.errors <- parallel::parLapply(cl, full.comb, DBerrS.P_ucpp, par.draws, 
                                       n.alts, i.cov, weights)
      parallel::stopCluster(cl)
    } else {
      # For each potential set, select best. 
      db.errors <- lapply(full.comb, DBerrS.P_ucpp, par.draws, n.alts, i.cov,
                          weights)
    }
    dbs <- unlist(db.errors, use.names = FALSE)
    set <- full.comb[[which.min(dbs)]]
    row.names(set) <- NULL
    if (n.cte > 0) {
      colnames(set) <- c(paste("alt",1:n.cte,".cte", sep = ""), c.nam)
    }
    db <- min(dbs)
    #return best set and db error design.
    return(list("set" = set, "error" = db))
  }
}