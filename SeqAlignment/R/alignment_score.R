#' Compute the alignment score.
#' 
#' @usage alignment_score(al1, al2, match, mismatch, gap)
#' @param al1 Aligned sequence 1.
#' @param al2 Aligned sequence 2.
#' @param match An integer which indicates the score for a match between 
#' the two sequences.
#' @param mismatch An integer which indicates the score for a mismatch
#' between the two sequences.
#' @param gap An integer which indicates the score for a gap inserted 
#' in one of the two sequences
#' @return The score of the alignment
#' @examples 
#' 
#' alignment_score('ACGGTTGC','A-GCGT-C', 1, -1, -2)
#' 
#' @export
alignment_score <- function(al1, al2, match, mismatch, gap){
  al1 <- unlist(strsplit(al1,""))
  al2 <- unlist(strsplit(al2,""))
  score <- 0
  
  #Check on the gap score
  if(gap > 0){
    warning('You have inserted a gap score greater than 0')
  }
  
  #Check that mismatch score is greater than match score
  if(match < mismatch){
    warning('Your mismatch score is greater than your match score')
  }
  
  #Compute the score of the two aligned sequences
  for (i in seq(1,length(al1))){
    if(al1[i]==al2[i]){
      score <- score + match
    }else if(al1[i]!=al2[i] & al1[i] != '-' & al2[i] != '-'){
      score <- score + mismatch
    }else{
      score <- score + gap
    }
  }
  return(score)
}