

#' Modified Fedorov algorithm for MNL models.
#' 
#' The algorithm swaps every profile of an initial start design with candidate 
#' profiles. By doing this, it tries to minimize the D(B)-error, based on a 
#' multinomial logit model. This routine is repeated for multiple starting 
#' designs.
#' 
#' Each iteration will loop through all profiles from the initial design, 
#' evaluating the change in A(B) or D(B)-error (as specified) for every profile 
#' from \code{cand.set}. The algorithm stops when an iteration occurred without 
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
#' The AB/DB-error is calculated by taking the mean over A/D-errors, respectively. 
#' It could be that for some draws the design results in an infinite error. 
#' The percentage of draws for which this was true for the final design can be 
#' found in the output \code{inf.error}.
#' 
#' Alternative specific constants can be specified in \code{alt.cte}. The length
#' of this binary vector should equal \code{n.alts}, were \code{0} indicates the
#' absence of an alternative specific constant and \code{1} the opposite.
#' 
#' \code{start.des} is a list with one or several matrices  corresponding to 
#' initial start design(s). In each matrix, each
#' row is a profile. The number of rows equals \code{n.sets * n.alts}, and the
#' number of columns equals the number of columns of \code{cand.set} + the
#' number of non-zero elements in \code{alt.cte}. If \code{start.des
#' = NULL}, \code{n.start} random initial designs will be
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
#' @param cand.set A numeric matrix in which each row is a possible profile. The
#'   \code{\link{Profiles}} function can be used to generate this matrix.
#' @param n.sets Numeric value indicating the number of choice sets.
#' @param n.alts Numeric value indicating the number of alternatives per choice 
#'   set.
#' @param alt.cte A binary vector indicating for each alternative whether an 
#'   alternative specific constant is desired. The default is \code{NULL}.
#' @param par.draws A matrix or a list, depending on \code{alt.cte}.
#' @param no.choice A logical value indicating whether a no choice alternative 
#'   should be added to each choice set. The default is \code{FALSE}.
#' @param start.des A list containing one or more matrices corresponding to initial start design(s). The default is \code{NULL}.
#' @param parallel Logical value indicating whether computations should be done 
#'   over multiple cores. The default is \code{TRUE}.
#' @param max.iter A numeric value indicating the maximum number allowed 
#'   iterations. The default is \code{Inf}.
#' @param n.start A numeric value indicating the number of random start designs
#'   to use. The default is 12.
#' @param optim A character value to choose between "D" and "A" optimality. The default is \code{"D"}.
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
#'   \item{SD}{The standrad deviation of the parameters. Calculated by taking the diagonal of the varcov matrix, averaged over all draws if a sample matrix was provided in \code{par.draws}.} 
#'   \item{level.count}{The count of all levels of each attribute in the design.} 
#'   \item{level.overlap}{The count of overlapping levels accross alternatives in every choice set in the design.} 
#'   \item{Orthogonality}{Numeric value indicating the degree of orthogonality of the design. The closer the value to 1, the more orthogonal the design is.} 
#'   \item{Blocks}{A list showing the created blocks of the best design, along with the level counts in each block. For more details, see function \code{\link{Blocks}}.} 
#' @examples
#' \dontrun{
#' # DB-efficient designs
#' # 3 Attributes, all dummy coded. 1 alternative specific constant = 7 parameters
#' cand.set <- Profiles(lvls = c(3, 3, 3), coding = c("D", "D", "D"))
#' mu <- c(0.5, 0.8, 0.2, -0.3, -1.2, 1.6, 2.2) # Prior parameter vector
#' v <- diag(length(mu)) # Prior variance.
#' set.seed(123) 
#' pd <- MASS::mvrnorm(n = 10, mu = mu, Sigma = v) # 10 draws.
#' p.d <- list(matrix(pd[,1], ncol = 1), pd[,2:7])
#' Modfed(cand.set = cand.set, n.sets = 8, n.alts = 2, alt.cte = c(1, 0),
#'        parallel = FALSE, par.draws = p.d)
#' # Or AB-efficient design
#' set.seed(123)
#' Modfed(cand.set = cand.set, n.sets = 8, n.alts = 2, alt.cte = c(1, 0),
#'        parallel = FALSE, par.draws = p.d, optim = "A")
#' 
#' # DB-efficient design with start design provided.  
#' # 3 Attributes with 3 levels, all dummy coded (= 6 parameters).
#' cand.set <- Profiles(lvls = c(3, 3, 3), coding = c("D", "D", "D")) 
#' mu <- c(0.8, 0.2, -0.3, -0.2, 0.7, 0.4) # Prior mean (total = 5 parameters).
#' v <- diag(length(mu)) # Prior variance.
#' sd <- list(example_design)
#' set.seed(123)
#' ps <- MASS::mvrnorm(n = 10, mu = mu, Sigma = v) # 10 draws.
#' Modfed(cand.set = cand.set, n.sets = 8, n.alts = 2, 
#'        alt.cte = c(0, 0), parallel = FALSE, par.draws = ps, start.des = sd)
#'        
#' # DB-efficient design with partial profiles
#' # 3 Attributes, all dummy coded. = 5 parameters
#' cand.set <- Profiles(lvls = c(3, 3, 2), coding = c("D", "D", "D")) 
#' mu <- c(1.2, 0.8, 0.2, -0.3, -1.2) # Prior parameter vector
#' v <- diag(length(mu)) # Prior variance.
#' set.seed(123) 
#' pd <- MASS::mvrnorm(n = 10, mu = mu, Sigma = v) # 10 draws.
#' Modfed(cand.set = cand.set, par.draws = pd, n.alts = 2,
#' n.sets = 8, parallel = TRUE, alt.cte = c(0, 0), overlap = 1)
#' # The same function but asking for blocks (and no overlap)
#' set.seed(123)
#' Modfed(cand.set = cand.set, par.draws = pd, n.alts = 2,
#' n.sets = 8, parallel = TRUE, alt.cte = c(0, 0), n.blocks = 2)
#' 
#' # AB-efficient design with constraints
#' # 2 dummy coded attributes, 1 continuous attribute and 1 effect coded
#' # attribute (with 4 levels). = 8 parameters
#' cand.set <- Profiles(lvls = c(3, 3, 2, 4), coding = c("D", "D", "C", "E"),
#'      c.lvls = list(c(4,7))) 
#' mu <- c(1.2, 0.8, 0.2, 0.5, -0.3, -1.2, 1, 1.6) # Prior parameter vector
#' v <- diag(length(mu)) # Prior variance.
#' set.seed(123) 
#' pd <- MASS::mvrnorm(n = 10, mu = mu, Sigma = v) # 10 draws.
#' constraints <- list("Alt2.Att1 = list(1A,1B)",
#'                     "if Alt1.Att3 = list(4) then Alt2.Att4 = list(4C, 4D)")
#' Modfed(cand.set = cand.set, par.draws = pd, n.alts = 2, n.sets = 8, 
#' parallel = TRUE, alt.cte = c(0, 0), optim = "A", constraints = constraints)
#'}
#' @importFrom Rdpack reprompt 
#' @importFrom MASS mvrnorm
#' @references \insertRef{idefix}{idefix}
#' @export
Modfed <- function(cand.set, n.sets, n.alts, par.draws, optim = "D", alt.cte = NULL, 
                   no.choice = FALSE, start.des = NULL, parallel = TRUE, max.iter = Inf, 
                   n.start = 12, overlap = NULL, n.blocks = 1, blocking.iter = 50, constraints = NULL) {
  
  if (is.null(alt.cte)) {
    alt.cte <- rep(0L, n.alts)
  }
  #init
  n.cte <- length(which(alt.cte == 1))
  ### Errors
  if (!is.list(par.draws)) {
    if (is.vector(par.draws)) {
      par.draws <- matrix(par.draws, nrow = 1)
    }
  }
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
    ncsek <- seq(n.alts, (n.sets * n.alts), n.alts) 
  } else {
    ncsek <- NULL
  }
  
  # Handling par.draws with alternative specific constants.
  if (isTRUE(all.equal(n.cte, 1))) {
    if (!(is.list(par.draws))) {
      stop("par.draws should be a list")
    }
    if (!isTRUE(all.equal(length(par.draws), 2))) {
      stop("'par.draws' should contain two components")
    }
    if (is.vector(par.draws[[1]])) {
      par.draws[[1]] <- matrix(par.draws[[1]], ncol = 1)
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
    if (!identical((dims[2, 1] + dims[2, 2]), (n.cte + ncol(cand.set)))) { 
      stop("the sum of the number of columns in the components of 'par.draws' 
           should equal the number of columns of 'cand.set' + the number of non-zero elements in 'alt.cte'")
    }
    par.draws  <- do.call("cbind", par.draws)
  }
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
      stop("the first component of 'par.draws' should contain the same number 
           of columns as there are non zero elements in 'alt.cte'")
    }
    dims <-  as.data.frame(lapply(par.draws, dim))
    if (!isTRUE(all.equal(dims[1, 1], dims[1, 2]))) { 
      stop("the number of rows in the components of 'par.draws' should be equal")
    }
    if (!identical((dims[2, 1] + dims[2, 2]), (n.cte + ncol(cand.set)))) { 
      stop("the sum of the number of columns in the components of 'par.draws' 
           should equal the number of columns of 'cand.set' + the number of non-zero elements in 'alt.cte'")
    }
    par.draws  <- do.call("cbind", par.draws)
  }
  # Create alternative specific design.
  cte.des <- Altspec(alt.cte = alt.cte, n.sets = n.sets)
  # Error identifying model.
  if (n.sets < ncol(par.draws)) {
    stop("Model is unidentified. Increase the number of choice sets or decrease parameters to estimate.")
  }
  # Handling cand.set
  if (!all(is.finite(cand.set))) {
    stop("'cand.set' contains non finite values.")
  }
  # Error handling cte.des
  if (ncol(cand.set) + ncol(cte.des) != ncol(par.draws)) {
    stop("The number of parameters in the components of 'par.draws' does not match the number 
         of non-zero parameters in 'alt.cte' + the number of parameters in 'cand.set'.")
  }
  # Random start design.
  if (!is.null(start.des)) {
    if (!is.list(start.des)) {
      stop("'start.des' should be a list")
    }
    if (!(all(unlist(lapply(start.des, is.matrix))))) {
      stop("'start.des' should contain matrices as components")
    }
    dimstart <- as.matrix(lapply(start.des, dim))
    nr.starts <- length(dimstart)
    if (nr.starts > 1.5) {
      if (!isTRUE(all.equal(length(unique(unlist(dimstart))), 2))) {
        stop("start designs have different dimensions")
      }
    }
    if (!isTRUE(all.equal(n.alts * n.sets, unique(unlist(dimstart))[1]))) {
      stop("number of rows of start design(s) does not match with 'n.alts' * 'n.sets'")
    }
    if (!isTRUE(all.equal(sum(ncol(cand.set), ncol(cte.des)), 
                          unique(unlist(dimstart))[2]))) {
      stop("number of columns of start design(s) does not match with the number
           of columns in 'cand.set' + the non zero parameters in 'alt.cte'")
    }
    d.start <- lapply(start.des, StartDB, par.draws, n.alts)
    if (!any(is.finite(unlist(lapply(d.start, mean, na.rm = TRUE))))) {
      stop("One or more of the provided start designs resulted in an unvalid db-error.
           Make sure the utility of any alternative cannot exceed 700.")
    }
  } 
  
  if (is.null(start.des)) {
    #create start designs
    nr.starts <- n.start
    start.des <- vector(mode = 'list', length = nr.starts)
    okstart <- FALSE
    while (okstart == FALSE) {
      for (i in 1:nr.starts) {
        if (!is.null(overlap) && overlap > 0) {
          #with overlapping
          start.des[[i]] <- do.call(rbind, lapply(1:n.sets, function(x) create_overlapping_cs(cand.set, overlap, n.alts)))  
          start.des[[i]] <- cbind(cte.des, data.matrix(start.des[[i]]))
        } else {
          #without overlapping
          r <- round(stats::runif((n.sets * n.alts), 1, nrow(cand.set)))
          start.des[[i]] <- cbind(cte.des, data.matrix(cand.set[r, ]))
        }
        
        #add no choice alternative
        if (no.choice) {
          start.des[[i]][ncsek, (ncol(cte.des) + 1):(ncol(cte.des) + ncol(cand.set))] <- c(rep(0, ncol(cand.set)))
        }
      
      #constraints####
      if (!is.null(constraints)) {
        sets_indeces <- ((c(1:nrow(start.des[[1]])) - 1) %/% n.alts) + 1 #get the choice set number
        sets_row_start <- unique((sets_indeces - 1) * n.alts + 1) #get the start row in the choice sets
        sets_row_end <- unique(ifelse(rep(no.choice,length(sets_indeces)), sets_indeces * n.alts - 1, sets_indeces * n.alts)) #get the end row in the choice sets
        n.par <- ncol(start.des[[1]])
        
        #first translate the constraints on coded variables
        constraints <- lapply(constraints, function(x) {translate_coded_constraints(x, coding = NULL , lvls = NULL, cand.set = cand.set)} )
        
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
            if (!is.null(overlap) && overlap > 0) {
              if (no.choice) {
                overlapping_set <- do.call(rbind, lapply(1, function(x) create_overlapping_cs(cand.set, overlap, n.alts)))[-1,]  
              } else {
                overlapping_set <- do.call(rbind, lapply(1, function(x) create_overlapping_cs(cand.set, overlap, n.alts)))
              }
              start.des[[i]][c(sets_row_start[s]:sets_row_end[s]), (n.cte + 1):n.par] <- overlapping_set
            } else {
              draws <- ifelse(no.choice, n.alts-1, n.alts)
              r <- round(stats::runif(draws, 1, nrow(cand.set)))
              start.des[[i]][c(sets_row_start[s]:sets_row_end[s]), (n.cte + 1):n.par] <- data.matrix(cand.set[r, ])
            }
            check <- all(unlist(lapply(constraints, function(x) {check_constraints(start.des[[i]],x,sets_row_start[s],n.cte)})))
            counter <- counter + 1
            if (counter == 1000) {
              message("Exceeded 1000 trials to satisfy the constraints. Consider stopping the run and modifying them.")
            }
          }
        } #2
      }
    }
      ##
      
      if(optim %in% c("D","d") ) {
        d.start <- lapply(start.des, StartDB, par.draws, n.alts) #DB-error
      } else {
        d.start <- lapply(start.des, StartAB, par.draws, n.alts) #AB-error
      }
      
      if(!any(is.finite(unlist(lapply(d.start, mean, na.rm = TRUE))))){
        stop("One or more draws resulted in an invalid d-error. 
             Make sure the utility of any alternative cannot exceed 700.")
      } else {
        okstart <- TRUE
      } 
    }
  }
  if (parallel) {
    #####
    no_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(no_cores)
    parallel::clusterExport(cl, c("n.sets", "par.draws", "cand.set", "n.alts", "n.cte", "alt.cte", "no.choice", "max.iter","ncsek", "overlap", "optim", "constraints"), envir = environment())
    deslist <- parallel::parLapply(cl, start.des, Modfedje_ucpp, par.draws, cand.set, n.alts, n.sets, n.cte, alt.cte, no.choice, max.iter, ncsek, overlap, optim, constraints)
    parallel::stopCluster(cl)
    ####
  } else {
    deslist <- lapply(start.des, Modfedje_ucpp, par.draws, cand.set, n.alts, n.sets, n.cte, alt.cte, no.choice, max.iter, ncsek, overlap, optim, constraints)
  }
  
  #measures/stats####
  #1-errors
  ab.error <- unlist(lapply(deslist, function(x) {
    mean(apply(par.draws, 1, Aerr_ucpp, des = x$design, n.alts = n.alts), na.rm = TRUE)  } ))
  names(ab.error) <-paste0("des",1:length(start.des))
  db.error <- unlist(lapply(deslist, function(x) {
    mean(apply(par.draws, 1, Derr_ucpp, des = x$design, n.alts = n.alts), na.rm = TRUE)  } ))
  names(db.error) <- paste0("des",1:length(start.des))
  
  #2- SD of the parameter estimates
  SD <- lapply(deslist, function(x) {
    rowMeans(apply(par.draws, 1, SD, des = x$design, n.alts = n.alts), na.rm = TRUE)
  })
  SD.mat <- do.call(rbind, SD)
  rownames(SD.mat) <- paste0("des",1:length(start.des))
  
  #3-orthogonality
  orthogonality <- unlist(lapply(deslist, function(x) {
    mean(apply(par.draws, 1, Orthogonality, des = x$design, n.alts = n.alts), na.rm = TRUE)}))
  names(orthogonality) <- paste0("des",1:length(start.des))
  
  #4-level overlaps
  if (no.choice) {
    lvl.overlap <- lapply(deslist, function(x) {lvl.overlap(x$design[-ncsek,(ncol(cte.des) + 1):(ncol(cte.des) + ncol(cand.set))], n.alts-1) } )   
  } else {
    lvl.overlap <- lapply(deslist, function(x) {lvl.overlap(x$design[,(ncol(cte.des) + 1):(ncol(cte.des) + ncol(cand.set))], n.alts) } )   
  }
  lvl.overlap_final <- lvl.overlap[[1]][, 1, drop=FALSE]
  for (i in 1:length(start.des)) {
    lvl.overlap_final <- cbind(lvl.overlap_final,lvl.overlap[[i]][, 2])
  }
  colnames(lvl.overlap_final)[2:(length(start.des)+1)] = paste0("des",1:length(start.des))
  
  #5-level frequencies
  if (no.choice) {
    lvl.freq <- lapply(deslist, function(x) {lvl.freq(x$design[-ncsek,(ncol(cte.des) + 1):(ncol(cte.des) + ncol(cand.set))]) } )   
    } else {
    lvl.freq <- lapply(deslist, function(x) {lvl.freq(x$design[,(ncol(cte.des) + 1):(ncol(cte.des) + ncol(cand.set))]) } )   
    }
  lvl.freq_final <- lvl.freq[[1]]
  n <- length(lvl.freq_final)
  for (a in 1:n) {
    if (length(start.des)>1) {
      for (i in 2:n.start) {
        temp <- dplyr::bind_rows(lvl.freq_final[[a]], (lvl.freq[[i]][[a]]))
        temp[is.na(temp)] <- 0 #replace NAs with zeros
        lvl.freq_final[[a]] <- temp 
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
                          AllDesigns = c( lapply(deslist, function(x) x[-length(x)] ) , stat)) #lapply is used to remove $error from each list to avoid duplication in stat
  } else {
    output_design <- list(BestDesign = c(bestdes[-length(bestdes)], bestdes_stat), 
                        AllDesigns = c( lapply(deslist, function(x) x[-length(x)] ) , stat)) #lapply is used to remove $error from each list to avoid duplication in stat
  }
  ##
  class(output_design) <- 'design_list'
  return(output_design)
}


# Core of the Modfed algorithm
Modfedje_ucpp <- function(desje, par.draws, cand.set, n.alts, n.sets, n.cte, alt.cte,
                          no.choice, max.iter, ncsek, overlap, optim, constraints){
  
  converge <- FALSE
  change <- FALSE
  it <- 1
  n.samples <- nrow(par.draws)
  n.par <- ncol(desje)
  num_attributes <- as.numeric(max(categorize_variables(colnames(cand.set))))
  ###
  while (!converge & it <= max.iter) {
    if(optim %in% c("D","d") ) {
      db.start <- mean(apply(par.draws, 1, Derr_ucpp, des = desje,  n.alts = n.alts), na.rm = TRUE) #DB-error
    } else {
      db.start <- mean(apply(par.draws, 1, Aerr_ucpp, des = desje,  n.alts = n.alts), na.rm = TRUE) #AB-error
    }
    it <- it + 1
    # save design before iteration.
    iter.des <- desje
    # For every row in the design.
    sek <- 1:nrow(desje)
    if (no.choice) {
      sek <- sek[-ncsek]
    }
    for (r in sek) {
      db <- numeric(nrow(cand.set))
      
      #Initialize the set index and the starting and ending row (to be used for overlapping/constraints)
      set_index <- ((r - 1) %/% n.alts) + 1 #get the choice set number
      set_row_start <- (set_index - 1) * n.alts + 1 #get the start row in the choice set
      set_row_end <- ifelse(no.choice, set_index * n.alts - 1, set_index * n.alts) #get the end row in the choice set (and avoid altering the no choice row)
        
      ### beginning of loop for overlap
      if(!is.null(overlap) && overlap > 0) {
        
        # Get all combinations of the variables to overlap
        attribute_combinations <- utils::combn(num_attributes, overlap, simplify = FALSE)
        
        
        for (combo in attribute_combinations) {
          # combo <- attribute_combinations[[1]]
          overlap_indices <- combo
          overlap_indices <- which(categorize_variables(colnames(cand.set)) %in% overlap_indices) + n.cte #extend the filtering to all the relevant columns for the overlapped variables (in case of dummy/effect treatments with more than 2 levels)
          
          temp_desje <- list()
          for (c in 1:nrow(cand.set)) {
            temp_desje[[c]] <- desje
            temp_desje[[c]][r, (n.cte + 1):n.par] <- cand.set[c, ]
            
            # Apply the overlapping values to all alternatives in the same choice set
            for (alt in set_row_start:set_row_end) {
              if (alt != r) {
                temp_desje[[c]][alt, overlap_indices] <- temp_desje[[c]][r, overlap_indices]
              }
            }
            
            # Check if the constraints are met, and skip to the next row if negative
            if (!is.null(constraints)) {
              check <- all(unlist(lapply((constraints), function(x) {check_constraints(temp_desje[[c]], x, set_row_start, n.cte)})))
              if (!check) {
                db[c] <- NA
                next
              }
            }
            
            # Calculate D-errors
            if(optim %in% c("D","d") ) {
              # temp_score <- mean(apply(par.draws, 1, Derr_ucpp, des = temp_desje[[c]], n.alts = n.alts), na.rm = TRUE)
              db[c] <- mean(apply(par.draws, 1, Derr_ucpp, des = temp_desje[[c]], n.alts = n.alts), na.rm = TRUE)
            } else {
              # temp_score <- mean(apply(par.draws, 1, Aerr_ucpp, des = temp_desje[[c]], n.alts = n.alts), na.rm = TRUE)
              db[c] <- mean(apply(par.draws, 1, Aerr_ucpp, des = temp_desje[[c]], n.alts = n.alts), na.rm = TRUE)
            }
          }
          
          pr <- which.min(db)
          db <- suppressWarnings(min(db, na.rm = TRUE)) #if all errors are equal to NaN, db will be Inf. We suppress this warning here.
          
          # Check if this design is better
          # if (!is.na(db[c]) && !is.na(db.start)) {
          if (!is.na(db) && !is.na(db.start)) {
            if (db < db.start) {
              desje <- temp_desje[[pr]]
              db.start <- db
            }
          }
        }
        ### end of loop for overlap
      } else {
      # Switch with every row in candidate set. 
        for (c in 1:nrow(cand.set)) {
          desje[r, (n.cte + 1):n.par ] <- cand.set[c, ]
          # Check if the constraints are met, and skip to the next row if negative
          if (!is.null(constraints)) {
            check <- all(unlist(lapply((constraints), function(x) {check_constraints(desje,x,set_row_start,n.cte)})))
            if (!check) {
              db[c] <- NA #since this iteration does not satisfy the constraint
              next
            } 
          }
        
          # Calculate D-errors or A-errors.
          if(optim %in% c("D","d") ) {
            d.errors <- apply(par.draws, 1, Derr_ucpp, des = desje,  n.alts = n.alts)
          } else {
            d.errors <- apply(par.draws, 1, Aerr_ucpp, des = desje,  n.alts = n.alts)
          }
        
          # AB/DB-error. 
          db[c] <- mean(d.errors, na.rm = TRUE)
        }
      
        pr <- which.min(db)
        db <- suppressWarnings(min(db, na.rm = TRUE)) #if all errors are equal to NaN, db will be Inf. We suppress this warning here.
      
        # Change if lower db error.
        if (!is.na(db) && !is.na(db.start)) {
          if (db < db.start) {
            best.row <- as.numeric(cand.set[pr, ])
            db.start <- db
            change <- TRUE
          }
        }
        # Replace with best profile if change.
        if (change) {
          desje[r, (n.cte + 1):n.par] <- best.row
        } else {
          desje[r, ] <- iter.des[r, ]
        }
        # Initialize variables again. 
      }
      change <- FALSE
      na.percentage <- 0
    }
    converge <- isTRUE(all.equal(desje, iter.des)) # Convergence if no profile is swapped this iteration.
  }
  # calculate percentage NA values.
  if(optim %in% c("D","d") ) {
    d.errors <- apply(par.draws, 1, Derr_ucpp, des = desje,  n.alts = n.alts)
  } else {
    d.errors <- apply(par.draws, 1, Aerr_ucpp, des = desje,  n.alts = n.alts)
  }
  if (any(is.na(d.errors))) {
    na.percentage <- scales::percent(sum(is.na(d.errors)) / n.samples)
  } 
  # Utility balance.
  # ub <- apply(par.draws, 1, Utbal, des = desje,  n.alts = n.alts)
  # ubcpp <- apply(par.draws, 1, InfoDes_cpp, des = desje,  n_alts = n.alts, 
  #                utbal = TRUE)
  
  # Utility balance using c++ function
  ub <- apply(par.draws, 1, InfoDes_cpp, des = desje,  n_alts = n.alts, 
              utbal = TRUE)
  pmat <- matrix(rowMeans(ub), ncol = n.alts, byrow = TRUE)
  rownames(pmat) <- paste("set", 1:n.sets, sep = "")
  colnames(pmat) <- paste(paste("Pr(", paste("alt", 1:n.alts, sep = ""), 
                                sep = ""), ")", sep = "")
  if (no.choice) {
    colnames(pmat)[n.alts] <- "Pr(no choice)"
  }
  # Rownames design. 
  des.names <- Rcnames(n.sets = n.sets, n.alts = n.alts, alt.cte = alt.cte, no.choice = no.choice)
  rownames(desje) <- des.names[[1]]
  # Colnames alternative specific constants. 
  if (n.cte != 0 && !is.null(colnames(desje))) {
    colnames(desje)[1:n.cte] <- des.names[[2]]
  }
  # Return design, D(B)error, percentage NA's, utility balance. 
  return(list("design" = desje, "inf.error" = na.percentage, "probs" = pmat, "optimality" =  optim,
              "error" =  db.start
              ))
}




#' Sequential modified federov algorithm for MNL model.
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
#' select te most efficient choice set based on the fisher information of the
#' prior covariance matrix \code{prior.covar}.
#' 
#' If \code{alt.cte = NULL}, \code{par.draws} should be a matrix in which each 
#' row is a sample from the multivariate parameter distribution. In case that 
#' \code{alt.cte} is not \code{NULL}, a list containing two matrices should be 
#' provided to \code{par.draws}. The first matrix containing the parameter draws
#' for the alternative specific parameters. The second matrix containing the
#' draws for the rest of the parameters.
#' 
#' The list of potential choice sets are created using 
#' \code{\link[utils]{combn}}. If \code{reduce} is \code{TRUE}, 
#' \code{allow.rep = FALSE} and vice versa. Furthermore, the list of 
#' potential choice sets will be screaned in order to select only those choice 
#' sets with a unique information matrix. If no alternative specific constants are used, 
#' \code{reduce} should always be \code{TRUE}. When alternative specific 
#' constants are used \code{reduce} can be \code{TRUE} so that the algorithm 
#' will be faster, but the combinations of constants and profiles will not be 
#' evaluated exhaustively.
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
#' *Note:* this function is more stable than \code{\link[idefix]{SeqCEA}}, but 
#' it takes more time to get the output. This happens because this function 
#' makes an exhaustive search to get the choice set, whereas 
#' \code{\link[idefix]{SeqCEA}} makes a random search.
#' 
#' @inheritParams Modfed
#' @param par.draws A matrix or a list, depending on \code{alt.cte}. 
#' @param des A design matrix in which each row is a profile. If alternative 
#'   specific constants are present, those should be included as the first 
#'   column(s) of the design. Can be generated with \code{\link{Modfed}} or \code{\link{CEA}}.
#' @param prior.covar Covariance matrix of the prior distribution.
#' @param weights A vector containing the weights of the draws. Default is 
#'   \code{NULL}, See also \code{\link{ImpsampMNL}}.
#' @param parallel Logical value indicating whether computations should be done 
#'   over multiple cores.
#' @param no.choice An integer indicating the no choice alternative. The default
#'   is \code{NULL}.
#' @param reduce Logical value indicating whether the candidate set should be 
#'   reduced or not.
#' @param allow.rep Logical value indicating whether repeated choice sets are
#' allowed in the design.
#' @return \item{set}{A matrix representing a DB efficient choice set.} 
#'   \item{error}{A numeric value indicating the DB-error of the whole 
#'   design.}
#' @importFrom Rdpack reprompt
#' @references \insertRef{idefix}{idefix}
#' @references \insertRef{ju}{idefix}
#' @examples 
#' # DB efficient choice set, given a design and parameter draws. 
#' # Candidate profiles 
#' cs <- Profiles(lvls = c(3, 3, 3), coding = c("D", "D", "D"))
#' m <- c(0.3, 0.2, -0.3, -0.2, 1.1, 2.4) # mean (total = 6 parameters).
#' pc <- diag(length(m)) # covariance matrix
#' set.seed(123)
#' sample <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc)
#' # Initial design.
#' des <- example_design 
#' # Efficient choice set to add. 
#' SeqMOD(des = des, cand.set = cs, n.alts = 2, par.draws = sample, 
#'            prior.covar = pc, parallel = FALSE)
#' 
#' # DB efficient choice set, given parameter draws. 
#' # with alternative specific constants 
#' des <- example_design2 
#' cs <- Profiles(lvls = c(3, 3, 3), coding = c("D", "D", "D"))
#' ac <- c(1, 1, 0) # Alternative specific constants. 
#' m <- c(0.3, 0.2, -0.3, -0.2, 1.1, 2.4, 1.8, 1.2) # mean 
#' pc <- diag(length(m)) # covariance matrix
#' pos <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc)
#' sample <- list(pos[ , 1:2], pos[ , 3:8])
#' # Efficient choice set. 
#' SeqMOD(des = des, cand.set = cs, n.alts = 3, par.draws = sample, alt.cte = ac, 
#'            prior.covar = pc, parallel = FALSE)
#' @export
SeqMOD <- function(des = NULL, cand.set, n.alts, par.draws, prior.covar, 
                  alt.cte = NULL, no.choice = NULL, weights = NULL, 
                  parallel = TRUE, reduce = TRUE, 
                  allow.rep = FALSE) {
  #init
  if (is.null(des)) {
    n.sets <- 1L
    allow.rep <- TRUE
  } else { 
    if (!is.matrix(des)) {
      stop("'des' should be a matrix or NULL")
    }
    if (!isTRUE(nrow(des) %% n.alts == 0)) {
      stop("'n.alts' does not seem to match with the number of rows in 'des'")
    }
    n.sets <- nrow(des) / n.alts
  }
  # if alternative constants 
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
  } else {
    n.cte <- 0
  }
  #if no.choice
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
  if (!is.null(alt.cte)) {
    #prior.covar
    if (!isTRUE(all.equal(ncol(prior.covar), (ncol(cand.set) + n.cte)))) {
      stop("number of columns of 'prior.covar' does not equal 
           the number of columns in 'cand.set' + nonzero elements in 'alt.cte'")
    }
    if (isTRUE(all.equal(n.cte, 1))) {
      if (!is.list(par.draws)) {stop("'par.draws' should be a list when 'alt.cte' is not NULL")}
      if (!isTRUE(all.equal(length(par.draws), 2))) {
        stop("'par.draws' should contain two components")
      }
      if (is.vector(par.draws[[1]])) {
        par.draws[[1]] <- matrix(par.draws[[1]], ncol = 1)
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
      if (!identical((dims[2, 1] + dims[2, 2]), (n.cte + ncol(cand.set)))) { 
        stop("the sum of the number of columns in the components of 'par.draws' 
             should equal the number of columns of 'cand.set' + the number of non-zero elements in 'alt.cte'")
      }
      par.draws  <- do.call("cbind", par.draws)
    }
    if (n.cte > 1.2) {
      if (!(is.list(par.draws))) {stop("'par.draws' should be a list when 'alt.cte' is not NULL")} 
      if (!isTRUE(all.equal(length(par.draws), 2))) {
        stop("'par.draws' should contain two components")
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
      if (!identical((dims[2, 1] + dims[2, 2]), (n.cte + ncol(cand.set)))) { 
        stop("the sum of the number of columns in the components of 'par.draws' 
             should equal the number of columns of 'cand.set' + the number of non-zero elements in 'alt.cte'")
      }
      par.draws  <- do.call("cbind", par.draws)
    }
    # Create alternative specific design.
    cte.des <- Altspec(alt.cte = alt.cte, n.sets = n.sets)
    cte.set <- matrix(cte.des[1:n.alts, ], ncol = n.cte, byrow = FALSE)
  } else {cte.des = NULL}
  # if no alternative constants 
  if (!is.matrix(par.draws)) {
    stop("'par.draws'should be a matrix when 'alt.cte' = NULL")
  }
  #init
  n.par <- ncol(par.draws)
  if (!is.null(weights)) {
    if (!isTRUE(all.equal(length(weights), nrow(par.draws)))) {
      stop("length of 'weights' does not match number total number of rows in 'par.draws'")
    }
  }
  # If no weights, equal weights.
  if (is.null(weights)) {
    weights <- rep(1L, nrow(par.draws))
  }
  ## whenever a design is supplied 
  if (!is.null(des)) {
    # Error par.draws
    if (!isTRUE(all.equal(ncol(des), n.par))) {
      stop("Numbers of columns in 'par.draws' does not match the number of columns in 'des'")
    }
    # Starting and initializing values.
    i.cov <- solve(prior.covar)
    d.start <- apply(par.draws, 1, DerrC_ucpp, des = des,  n.alts = n.alts, 
                     i.cov = i.cov)
    db.start <- mean(d.start, na.rm = TRUE)
    full.comb <- Fullsets_ucpp(cand.set = cand.set, n.alts = n.alts, 
                               no.choice = no.choice, reduce = reduce, 
                               allow.rep = allow.rep, des = des, n.cte = n.cte)
    #if alt.cte
    if (!is.null(cte.des)) {
      full.comb <- lapply(full.comb, function(x) cbind(cte.set, x))
    }
    full.des <- lapply(full.comb, function(x) rbind(des, x))
    # For each potential set, select best.
    ##### parallel #####
    if (parallel) {
      no_cores <- parallel::detectCores() - 1L
      cl <- parallel::makeCluster(no_cores)
      # New line to copy DerrS.P from .GlobalEnv to the cluster
      # https://stackoverflow.com/questions/12023403/using-parlapply-and-clusterexport-inside-a-function
      #clusterExport(cl=cl,varlist=c("DerrS.P_ucpp"))
      db.errors <- parallel::parLapply(cl, full.des, DBerrS.P_ucpp, par.draws, 
                                       n.alts, i.cov, weights)
      parallel::stopCluster(cl)
      ##### parallel #####
    } else {
      # For each potential set, select best. 
      db.errors <- lapply(full.des, DBerrS.P_ucpp, par.draws, n.alts, i.cov, weights)
    }
    dbs <- unlist(db.errors, use.names = FALSE)
    set <- full.comb[[which.min(dbs)]]
    row.names(set) <- NULL
    db <- min(dbs)
    #return best set and db error design.
    return(list("set" = set, "error" = db))
  }
  
  if (is.null(des)) {
    
    # Starting and initializing values.
    i.cov <- solve(prior.covar)
    full.comb <- Fullsets_ucpp(cand.set = cand.set, n.alts = n.alts, 
                               no.choice = no.choice, reduce = reduce, 
                               allow.rep = allow.rep, des = des)
    #if alt.cte
    if (!is.null(cte.des)) {
      full.comb <- lapply(full.comb, function(x) cbind(cte.set, x))
    } else {
      if (!isTRUE(all.equal(ncol(cand.set), ncol(par.draws)))) {
        stop("number of  columns of 'par.draws' and 'cand.set' should be equal when 'alt.cte = NULL'")
      }
    }
    # For each potential set, select best.
    
    ##### parallel #####
    if (parallel) {
      no_cores <- parallel::detectCores() - 1L
      cl <- parallel::makeCluster(no_cores)
      # New line to copy DerrS.P from .GlobalEnv to the cluster
      # https://stackoverflow.com/questions/12023403/using-parlapply-and-clusterexport-inside-a-function
      # Becareful with cpp functions in cluster export
      # Check https://stackoverflow.com/questions/38518387/using-rcpp-functions-inside-of-rs-parapply-functions-from-the-parallel-package
      # https://stackoverflow.com/questions/25606733/using-rcpp-function-in-parlapply-on-windows/25606950
      #clusterExport(cl=cl,varlist=c("DerrS.P_ucpp"))
      db.errors <- parallel::parLapply(cl, full.comb, DBerrS.P_ucpp, par.draws, 
                                       n.alts, i.cov, weights)
      parallel::stopCluster(cl)
      ##### parallel #####
    } else {
      # For each potential set, select best. 
      db.errors <- lapply(full.comb, DBerrS.P_ucpp, par.draws, n.alts, i.cov,
                          weights)
    }
    dbs <- unlist(db.errors, use.names = FALSE)
    set <- full.comb[[which.min(dbs)]]
    row.names(set) <- NULL
    db <- min(dbs)
    #return best set and db error design.
    return(list("set" = set, "error" = db))
  }
}


#' Sequential Kullback-Leibler based algorithm for the MNL model.
#' 
#' Selects the choice set that maximizes the Kullback-Leibler divergence between
#' the prior parameter values and the expected posterior, assuming a MNL model.
#' 
#' This algorithm is ideally used in an adaptive context. The algorithm selects 
#' the choice set that maximizes the Kullback-Leibler 
#' divergence between prior and expected posterior. Otherwisely framed the 
#' algorithm selects the choice set that maximizes the expected information 
#' gain.
#' 
#' If \code{alt.cte = NULL}, \code{par.draws} should be a matrix in which each 
#' row is a sample from the multivariate parameter distribution. In case that 
#' \code{alt.cte} is not \code{NULL}, a list containing two matrices should be 
#' provided to \code{par.draws}. The first matrix containing the parameter draws
#' for the alternative specific parameters. The second matrix containing the
#' draws for the rest of the parameters.
#'
#' The list of potential choice sets are created using 
#' \code{\link[utils]{combn}}. The \code{weights} argument can be used when the
#'  \code{par.draws} have 
#' weights. This is for example the case when parameter values are updated using
#' \code{\link{ImpsampMNL}}.
#' 
#' 
#' @inheritParams SeqMOD
#' @param alt.cte A binary vector indicating for each alternative if an
#'   alternative specific constant is desired.
#' @return \item{set}{Numeric matrix containing the choice set that maximizes the expected KL divergence.}
#' \item{kl}{Numeric value which is the Kullback leibler divergence.}
#' @importFrom Rdpack reprompt
#' @references 
#' \insertRef{crabbe}{idefix}
#' @examples 
#' # KL efficient choice set, given parameter draws. 
#' # Candidate profiles 
#' cs <- Profiles(lvls = c(3, 3), coding = c("E", "E"))
#' m <- c(0.3, 0.2, -0.3, -0.2) # Prior mean (4 parameters).
#' pc <- diag(length(m)) # Prior variance
#' set.seed(123)
#' ps <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc) # 10 draws.
#' # Efficient choice set to add. 
#' SeqKL(cand.set = cs, n.alts = 2, alt.cte = NULL, par.draws = ps, weights = NULL)
#' 
#' # KL efficient choice set, given parameter draws. 
#' # Candidate profiles 
#' cs <- Profiles(lvls = c(3, 3), coding = c("C", "E"), c.lvls = list(c(5,3,1)))
#' m <- c(0.7, 0.3, -0.3, -0.2) # Prior mean (4 parameters).
#' pc <- diag(length(m)) # Prior variance
#' set.seed(123)
#' ps <- MASS::mvrnorm(n = 10, mu = m, Sigma = pc) # 10 draws.
#' sample <- list(ps[ , 1], ps[ , 2:4])
#' ac <- c(1, 0) # Alternative specific constant. 
#' # Efficient choice set to add. 
#' SeqKL(cand.set = cs, n.alts = 2, alt.cte = ac, par.draws = sample, weights = NULL)
#' @export
SeqKL <- function(des = NULL, cand.set, n.alts, par.draws, alt.cte = NULL, 
                  no.choice = NULL, weights = NULL, allow.rep = FALSE) {
  # Handling error initial design
  if (is.null(des)) {
    n.sets <- 1L
    allow.rep <- TRUE
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
  
  # Error alternative specific constants.
  if (!is.null(alt.cte)) {
    if (length(alt.cte) != n.alts) {
      stop("n.alts does not match the alt.cte vector")
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
    cte.des <- Altspec(alt.cte = alt.cte, n.sets = 1)  
    cte.set <- matrix(cte.des[1:n.alts, ], ncol = n.cte, byrow = FALSE)
    # Error handling cte.des
    if (ncol(cand.set) + ncol(cte.des) != ncol(par.draws)) {
      stop("dimension of par.draws does not match the dimension of alt.cte + cand.set.")
    }
  } else {
    n.cte <- 0
    cte.des <- NULL
    # Error handling cte.des
    if (ncol(cand.set) != ncol(par.draws)) {
      stop("dimension of par.draws does not match the dimension of alt.cte + cand.set.")
    }
  }
  
  # Check number of columns in par.draws and weights
  n.par <- ncol(par.draws)
  if (!is.null(weights)) {
    if (!isTRUE(all.equal(length(weights), nrow(par.draws)))) {
      stop("length of 'weights' does not match number total number of rows in 'par.draws'")
    }
  }
  # All choice sets.
  # full.comb <- gtools::combinations(n = nrow(cand.set), r = n.alts, 
  #                                   repeats.allowed = !reduce)
  full.comb <- Fullsets_ucpp(cand.set = cand.set, n.alts = n.alts, 
                                 no.choice = no.choice, reduce = FALSE, 
                                 allow.rep = allow.rep, des = des, n.cte = n.cte)
  
  # If no weights, equal weights.
  if (is.null(weights)) {
    weights <- rep(1, nrow(par.draws))
  }
  # Add alternative specific constants if necessary
  if (!is.null(cte.des)) {
    full.comb <- lapply(full.comb, function(x) cbind(cte.set, x))
  }
  # Calculate KL for each set.
  #kl.infos <- apply(full.comb, 1, KLs, par.draws, cte.des, cand.set, weights)
  kl.infos <- lapply(full.comb, KL, par.draws, weights)
  
  # Select maximum.
  #comb.nr <- as.numeric(full.comb[which.max(kl.infos), ])
  #set <- cand.set[comb.nr, ]
  set <- full.comb[[which.max(kl.infos)]]
  # Add alternative specific constants if necessary
  #if (!is.null(cte.des)) {
  #  set <- cbind(cte.des, set)
  #}
  row.names(set) <- NULL
  # return.
  return(list(set = set, kl = max(unlist(kl.infos))))
}
