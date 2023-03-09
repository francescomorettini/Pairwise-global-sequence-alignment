#' FOGSAA algorithm for global sequence alignment.
#' 
#' Computes the optimal global alignment of two sequences using the FOGSAA 
#' algorithm.
#' 
#' @usage fogsaa_alignment(s1, s2, match, mismatch, gap)
#' @param s1 A string representing sequence 1.
#' @param s2 A string representing sequence 2.
#' @param match An integer which indicates the score for a match between 
#' the two sequences.
#' @param mismatch An integer which indicates the score for a mismatch
#' between the two sequences.
#' @param gap An integer which indicates the score for a gap inserted 
#' in one of the two sequences.
#' @return A list containing the alignment of \code{s1} and \code{s2}
#' @references \url{https://pubmed.ncbi.nlm.nih.gov/23624407/}\cr
#' @examples 
#' fogsaa_alignment('ACGGTTGC','AGCGTC', 1, -1, -2)
#' @export
fogsaa_alignment <- function(s1, s2, match, mismatch, gap){
  m <- nchar(s1)
  n <- nchar(s2)
  ALPHABET <- c('A', 'C', 'G', 'T', '-')   #Possible characters in the string
  
  s1 <- unlist(strsplit(s1,''))
  s2 <- unlist(strsplit(s2,''))
  
  #Check if all the characters of the string are present in the reference alphabet
  if(any(!sapply(s1, is.element, ALPHABET))){
    stop('The inserted string contains characters not allowed')
  }
  
  if(any(!sapply(s2, is.element, ALPHABET))){
    stop('The inserted string contains characters not allowed')
  }
  
  #Stop if the gap score is greater than 0
  if(gap > 0){
    stop('You have inserted a gap score greater than 0')
  }
  
  #Stop if the match score is smaller than the mismatch one
  if(match < mismatch){
    stop('Your mismatch score is greater than your match score')
  }
  
  #Create the first node which is root, so both pointers are set to 0 as well as
  #the present score
  current_node <- node(P1 = 0, P2= 0)
  current_node$v1 <- ''
  current_node$v2 <- ''
  current_node$x1 <- n - current_node$P2
  current_node$x2 <- m - current_node$P1
  current_node$PrS <- 0
  
  #Computation of the Fmin and Fmax score 
  if(current_node$x2 < current_node$x1){
    current_node$Fmin <- current_node$x2 * mismatch +
      gap*(current_node$x1-current_node$x2)
  }else{
    current_node$Fmin <- current_node$x1 * mismatch +
      gap*(current_node$x2-current_node$x1)
  }
  
  if(current_node$x2 < current_node$x1){
    current_node$Fmax <- current_node$x2 * match +
      gap*(current_node$x1-current_node$x2)
  }else{
    current_node$Fmax <- current_node$x1 * match +
      gap*(current_node$x2-current_node$x1)
  }
  
  #Computation of Tmin and Tmax score
  current_node$Tmin <- current_node$PrS + current_node$Fmin
  current_node$Tmax <- current_node$PrS + current_node$Fmax
  
  
  
  current_branch <- list(current_node)  #The branch we are currently expanding
  queue <- list()    #The queue containing the nodes not yet expanded
  top_branch <- list()   #The current best branch
  lowerBound <- current_node$Tmin    #The lowerBound which is the score top_branch
  expanded_nodes <- list()  #A list containing all the nodes expanded
  
  while(current_node$P1 <= m-1 & current_node$P2 <= n-1){
    
    #Check if we have already expanded a better node than the current one
    if(any(sapply(expanded_nodes, compare, current_node))){
      #If we have check if there is any other node in the queue
      if(length(queue) != 0){
        #If there is take the first from the queue and use it as current node
        current_node <- queue[[1]]                           
        queue <- queue[-1]     #remove it from the queue                              
        current_branch <- list(current_node)   #Clean the current_branch
      }else{
        #If the queue is empty stop
        break
      }
      #If we update the current node we go to the next iteration and do the control
      #again
      next
    }else{
      #If there aren't better nodes go on and expand it
      expanded_nodes <- append(expanded_nodes, current_node)
    }
    
    #Three possible way to expand the current node
    #Increasing all the pointer or only one of them
    node_11 <- node(P1 = current_node$P1 + 1, P2 = current_node$P2 + 1,
                    PrS = current_node$PrS, code = '11')
    node_10 <- node(P1 = current_node$P1 + 1, P2 = current_node$P2,
                    PrS = current_node$PrS, code = '10')
    node_01 <- node(P1 = current_node$P1, P2 = current_node$P2 + 1,
                    PrS = current_node$PrS, code = '01')
    
    #If one pointer is greater than the length of the string insert a -
    #Compute v1 for every node
    if(node_11$P1 > m){
      node_11$v1 <- '-'
    }else{
      node_11$v1 <- s1[node_11$P1]
    }
    
    if(node_10$P1 > m){
      node_10$v1 <- '-'
    }else{
      node_10$v1 <- s1[node_10$P1]
    }
    
    node_01$v1 <- '-'
    
    #Compute v2 for every node
    
    if(node_11$P2 > n){
      node_11$v2 <- '-'
    }else{
      node_11$v2 <- s2[node_11$P2]
    }
    
    if(node_01$P2 > n){
      node_01$v2 <- '-'
    }else{
      node_01$v2 <- s2[node_01$P2]
    }
    
    node_10$v2 <- '-'
    
    #Compute x1 and x2 values for the new nodes
    node_11$x1 <- n - node_11$P2
    node_10$x1 <- n - node_10$P2
    node_01$x1 <- n - node_01$P2
    
    node_11$x2 <- m - node_11$P1
    node_10$x2 <- m - node_10$P1
    node_01$x2 <- m - node_01$P1
    
    #Compute the score for all the three nodes that we have generated
    compute_scores(node_11, current_node, match, mismatch, gap)
    compute_scores(node_10, current_node, match, mismatch, gap)
    compute_scores(node_01, current_node, match, mismatch, gap)
    
    #The new current node which will be expanded is the one with the highest score
    current_node <- max_score(node_11, node_10, node_01)
    
    #The other two nodes are inserted in the queue
    #The method equals look at v1 and v2
    if(current_node$equals(node_11)){
      queue <- append(queue, node_01)
      queue <- append(queue, node_10)
    }else if(current_node$equals(node_10)){
      queue <- append(queue, node_11)
      queue <- append(queue, node_01)
    }else{
      queue <- append(queue, node_11)
      queue <- append(queue, node_10)
    }
    
    #Update the current branch with the new current node 
    current_branch <- append(current_branch, current_node)
    
    #Oreder the queue based on the value of Tmax. In case of any draw we check
    #the value of Tmin
    queue <- queue[order(sapply(queue, function(x) x$Tmax),
                         (sapply(queue, function(x) x$Tmin)),
                         decreasing = TRUE)]
    
    #Remove from the queue all the nodes which don't have the potential anymore
    #to have a greater score than the current Lower bound
    queue <- queue[sapply(queue, function(x) x$Tmax >= lowerBound)]
    
    #If we are arrived at the end of the tree
    if(current_node$P1 == m & current_node$P2 == n){ 
      #and if the score of this branch is better than the previous top branch
      if(current_node$Tmax >= lowerBound){             
        lowerBound <- current_node$Tmax   #Update the lower bound
        top_branch <- current_branch    #Update the top branch
      }
      #If the queue is not empty select the best node from it and start a new branch
      if(length(queue) != 0){
        current_node <- queue[[1]]
        current_branch <- list(current_node)
        queue <- queue[-1]
      }else{
        #If the queue is empty break
        break
      }
      
    }
  }
  
  #Reordering the expanded nodes list
  expanded_nodes <- 
    expanded_nodes[order(sapply(expanded_nodes, function(x) x$Tmax),
                         (sapply(expanded_nodes, function(x) x$Tmin)), 
                         decreasing = TRUE)]
  
  #Arrived at this point in top branch we have saved only the last part of the tree
  #To retrieve the previous one up to root we have to use the attribute code of the nodes.
  #If the attribute of a node is '11' it means that it has been originated incrementing
  #both the pointers. The same apply for code '10' and '01'. Since we have reordered
  #the list we will retrieve the node from which we arrived
  while(top_branch[[1]]$P1 > 0 & top_branch[[1]]$P2 > 0){
    if(top_branch[[1]]$code == '11'){
      index_1 <- top_branch[[1]]$P1 - 1
      index_2 <- top_branch[[1]]$P2 - 1
      possible_nodes <- expanded_nodes[sapply(expanded_nodes, function(x) 
        x$P1 == index_1 & x$P2 == index_2)]
      top_branch <- append(possible_nodes[[1]], top_branch)
    }else if(top_branch[[1]]$code == '10'){
      index_1 <- top_branch[[1]]$P1 - 1
      index_2 <- top_branch[[1]]$P2
      possible_nodes <- expanded_nodes[sapply(expanded_nodes, function(x) 
        x$P1 == index_1 & x$P2 == index_2)]
      top_branch <- append(possible_nodes[[1]], top_branch)
    }else{
      index_1 <- top_branch[[1]]$P1
      index_2 <- top_branch[[1]]$P2 - 1
      possible_nodes <- expanded_nodes[sapply(expanded_nodes, function(x) 
        x$P1 == index_1 & x$P2 == index_2)]
      top_branch <- append(possible_nodes[[1]], top_branch)
    }
    
    
  }
  
  #Removing root
  top_branch <- top_branch[-1]
  #Return the two aligned strings
  aligned1 <- paste(sapply(top_branch, function(x) x$v1),collapse="")
  aligned2 <- paste(sapply(top_branch, function(x) x$v2),collapse="")
  return(list(aligned1, aligned2))
}


#' Computation of present and future score of each node.
#' 
#' It computes both the present and future score of each new node, finally it
#' sums the two to obtain the Fitness score which is used to identify the best
#' node
#' 
#' @param new_vertex The node for which we compute the scores
#' @param current_vertex The previous node from which we retrieve the PrS score
#' @param match An integer which indicates the score for a match between 
#' the two sequences.
#' @param mismatch An integer which indicates the score for a mismatch
#' between the two sequences.
#' @param gap An integer which indicates the score for a gap inserted 
#' in one of the two sequences.
#' @return The node in which all the attributes reguarding the different scores
#' are computed
compute_scores <- function(new_vertex, current_vertex, match, mismatch, gap){
  
  if(new_vertex$v1 == new_vertex$v2 & new_vertex$code == '11'){
    new_vertex$PrS <- current_vertex$PrS + match
  }else if(new_vertex$v1 != new_vertex$v2 & new_vertex$code == '11'){
    new_vertex$PrS <- current_vertex$PrS + mismatch
  }else{
    new_vertex$PrS <- current_vertex$PrS + gap
  }
  
  
  if(new_vertex$x2 < new_vertex$x1){
    new_vertex$Fmin <- new_vertex$x2 * mismatch + 
      gap*(new_vertex$x1-new_vertex$x2)
  }else{
    new_vertex$Fmin <- new_vertex$x1 * mismatch + 
      gap*(new_vertex$x2-new_vertex$x1)
  }
  
  if(new_vertex$x2 < new_vertex$x1){
    new_vertex$Fmax <- new_vertex$x2 * match + gap*(new_vertex$x1-new_vertex$x2)
  }else{
    new_vertex$Fmax <- new_vertex$x1 * match + gap*(new_vertex$x2-new_vertex$x1)
  }
  
  new_vertex$Tmin <- new_vertex$PrS + new_vertex$Fmin
  new_vertex$Tmax <- new_vertex$PrS + new_vertex$Fmax
  
  
}

#' Max score for a node.
#'
#' It computes the node which has the highest score among the three
#'
#' @param vertex_11 A node originated from the current node incrementing 
#' both P1 and P2
#' @param vertex_10 A node originated from the current node incrementing
#' only P1
#' @param vertex_01 A node originated from the current node incrementing
#' only P2
#' @return The node with the highest Fitness score
max_score <- function(vertex_11, vertex_10, vertex_01){
  score_11 <- vertex_11$Tmax
  score_10 <- vertex_10$Tmax
  score_01 <- vertex_01$Tmax
  max_score <- max(score_11, score_10, score_01)
  if(max_score == score_11){
    return(vertex_11)
  }else if(max_score == score_10){
    return(vertex_10)
  }else{
    return(vertex_01)
  }
  
}

#' Comparison of two nodes.
#'
#' It compares two nodes based on their pointers P1 and P2 and checking the Tmax
#' score
#' 
#' @param new_vertex A node
#' @param current_vertex A node
#' @return True if the nodes are equal in terms of P1 and P2 and the new node
#' has an higher score, false in any other case
compare <- function(new_vertex, current_vertex){
  if(new_vertex$P1 == current_vertex$P1 & new_vertex$P2 == current_vertex$P2){
    if(new_vertex$Tmax >= current_vertex$Tmax){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }else
    return(FALSE)
}
