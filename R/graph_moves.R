n.edge = function(A){
  length(which(A[lower.tri(A)] == 1 | t(A)[lower.tri(A)] == 1))
}


check_decomp <- function(A){
  g <- igraph::graph_from_adjacency_matrix(A, mode="undirected", diag=FALSE)
  return(igraph::is_chordal(g)$chordal)
}

# Two possible actions

actions = c("iu", "du")

# Insert an undirected edge

iu = function(A, nodes){
  x = nodes[1]
  y = nodes[2]
  A[x,y] = A[y,x] = 1
  return(A)
}

# Delete an undirected edge

du = function(A, nodes){
  x = nodes[1]
  y = nodes[2]
  A[x,y] = A[y,x] = 0
  return(A)
}


move_decomposable <- function(A) {
  
  q <- nrow(A)
  A <- (A != 0) * 1
  
  valid_moves <- list()
  
  for (i in 2:q) {
    for (j in 1:(i-1)) {
      
      A_new <- A
      
      if (A[i,j] == 0) {
        A_new[i,j] <- A_new[j,i] <- 1
        if (check_decomp(A_new))
          valid_moves[[length(valid_moves)+1]] <- list(type="add", nodes=c(i,j))
      }
      
      if (A[i,j] == 1) {
        A_new[i,j] <- A_new[j,i] <- 0
        if (check_decomp(A_new))
          valid_moves[[length(valid_moves)+1]] <- list(type="del", nodes=c(i,j))
      }
    }
  }
  
  if (length(valid_moves) == 0)
    stop("No valid decomposable move")
  
  m <- valid_moves[[ sample(length(valid_moves),1) ]]
  
  A_new <- A
  i <- m$nodes[1]
  j <- m$nodes[2]
  
  if (m$type == "add") {
    A_new[i,j] <- A_new[j,i] <- 1
  } else {
    A_new[i,j] <- A_new[j,i] <- 0
  }
  
  return(list(
    A_new = A_new,
    n_moves = length(valid_moves), # needed because the proposed graph could have a different number of moves wrt the current one 
    nodes    = c(i, j),           # (u, v) involved in the move
    type     = m$type             # "add" or "del"
    ))
}