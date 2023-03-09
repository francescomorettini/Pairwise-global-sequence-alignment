#' A Reference Class to represent the nodes used in the FOGSAA algorithm.
#' 
#' The class \code{node} is used to represent objects in the FOGSAA algorithm,
#' this is the data type on which the algorithm is implemented.
#' 
#' @name node
#' @field P1 It is the pointer to the characters in sequence 1.
#' @field P2 It is the pointer to the characters in sequence 2.
#' @field PrS Present score of the node.
#' @field Fmin Minimum of the future score.
#' @field Fmax Maximum of the future score.
#' @field Tmin Minimum of the Fitted score computed as the sum of PrS and Fmin.
#' @field Tmax Maximum of the Fitted score computed as the sum of PrS and Fmax.
#' @field v1 Value of sequence 1 at which P1 is pointing.
#' @field v2 Value of sequence 2 at which P2 is pointing.
#' @field x1 Defined as n - P2 where n is the length of sequence 2.
#' @field x2 Defined as m - P1 where m is the length of sequence 1.
#' @field code Indicates how the instance was originated, 11 if incrementing 
#' both P1 and P2 from the previous node, 01 if only P2 is incremented and 10 if 
#' only P1 is incremented.
#' @export
#' @importFrom methods setRefClass
#' @importFrom methods new
node <- setRefClass("node", 
                    fields = list(P1 = 'numeric',P2 = 'numeric',PrS = 'numeric',
                                  Fmin = 'numeric', Fmax = 'numeric',
                                  Tmin = 'numeric',Tmax = 'numeric', 
                                  v1 = 'character', v2 = 'character',
                                  x1 = 'numeric', x2 = 'numeric', 
                                  code = 'character'))


node$methods(equals = function(node1){
  if(node1$v1 == v1 & node1$v2 == v2){
    return(TRUE)
  }else{
    return(FALSE)
  }
})