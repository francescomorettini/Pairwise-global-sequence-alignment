#' Needlman-Wunsch global sequence alignment algorithm.
#' 
#' Computes the global sequence alignment using the Needleman-Wunsh algorithm.
#' 
#' @usage nw_alignment(s1, s2, match, mismatch, gap)
#' @param s1 A string representing sequence 1.
#' @param s2 A string representing sequence 2.
#' @param match An integer which indicates the score for a match between 
#' the two sequences.
#' @param mismatch An integer which indicates the score for a mismatch
#' between the two sequences.
#' @param gap An integer which indicates the score for a gap inserted 
#' in one of the two sequences.
#' @return A list containing the alignment of \code{s1} and \code{s2}
#' @references \url{https://pubmed.ncbi.nlm.nih.gov/5420325/}\cr
#' @examples 
#' 
#' nw_alignment('ACGGTTGC','AGCGTC', 1, -1, -2)
#' 
#' @export
#' @importFrom stringi stri_reverse
nw_alignment <- function(s1, s2, match, mismatch, gap){
  len1 <- nchar(s1)
  len2 <- nchar(s2)
  ALPHABET <- c('A', 'C', 'G', 'T', '-')     #Alphabet of the character admitted
  
  
  score_matrix <- matrix(0, nrow= len2 + 1, ncol = len1 + 1) #Scoring matrix initialisation
  trace_back_matrix <- matrix(0, nrow= len2 + 1, ncol = len1 + 1)  #Trace back matrix initialisation
  s1 <- unlist(strsplit(c('-',s1),""))  #Unlist both strings and adding a '-' as first character
  s2 <- unlist(strsplit(c('-',s2),""))
  
  #Check if all the characters of the string are present in the reference alphabet
  if(any(!sapply(s1, is.element, ALPHABET))){
    stop('The inserted string contains characters not allowed')
  }
  
  if(any(!sapply(s2, is.element, ALPHABET))){
    stop('The inserted string contains characters not allowed')
  }
  
  #Return a warning if the gap score is greater than 0
  if(gap > 0){
    warning('You have inserted a gap score greater than 0')
  }
  
  #Return a warning if the mismatch score is greater than the match score
  if(match < mismatch){
    warning('Your mismatch score is greater than your match score')
  }
  
  #Initialising the first row and the first column of both matrices.
  score_matrix[1,] <- seq(0, len1*gap, gap)
  score_matrix[,1] <- seq(0, len2*gap, gap)
  trace_back_matrix[1,-1] <- rep('left', len1)
  trace_back_matrix[-1,1] <- rep('up', len2)
  
  #Compute the score for the (i,j) cell
  for (i in seq(2,len2 + 1)){
    for (j in seq(2,len1 + 1)){
      from_left <- score_matrix[i, j-1] + gap   #Score coming from the left hand cell
      
      from_up <- score_matrix[i-1, j] + gap  #Score from the upper cell
      
      #Computation of the score coming from the diagonal cell
      if (s1[j] == s2[i]){  #Looking for a match or mismatch in the nucleotides 
        x = match           #in position i and j.
      }else {
        x = mismatch
      }
      
      from_diag <- score_matrix[i-1, j-1] + x   #Score from the diagonal cell
      
      max_value = max(from_left, from_up, from_diag)   #Taking the highest score of the three computed 
      score_matrix[i,j] = max_value                    #before
      
      #Adding the origin of the (i,j) cell's score in the trace back matrix
      if (max_value == from_diag){         #In case of any draw we prefer the diagonal path
        trace_back_matrix[i,j] = 'diag'    #because it means less gaps
      } else if (max_value == from_up){
        trace_back_matrix[i,j] = 'up'
      } else {
        trace_back_matrix[i,j] = 'left'
      }
      
      
    }
  }
  
  #Initialising the aligned strings
  aligned1 <- "" 
  aligned2 <- ""
  
  i = len1 + 1
  j = len2 + 1
  
  #Trace-back phase which continues until we are in the upper left cell
  while(i != 1 & j != 1){
    current <- trace_back_matrix[j,i]
    
    #From diag we add a character to both string
    if(current == 'diag'){
      aligned1 <- paste(aligned1, s1[i], sep = "")
      aligned2 <- paste(aligned2, s2[j], sep = "")
      i = i-1
      j = j-1
    }else if(current == 'up'){                          
      aligned1 <- paste(aligned1, "-", sep = "")     #One string will have a gap
      aligned2 <- paste(aligned2, s2[j], sep = "")
      j = j-1
    }else{
      aligned1 <- paste(aligned1, s1[i], sep = "")
      aligned2 <- paste(aligned2, "-", sep = "")
      i = i-1
    }
    
  }
  
  #Reverse the string in order to obtain the two aligned string
  ali1 <- stri_reverse(aligned1)
  ali2 <- stri_reverse(aligned2)
  return(list(ali1,ali2))
}