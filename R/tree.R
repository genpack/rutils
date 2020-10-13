# package gener:

# Header
# Filename:     tree.R
# Description:  Defines a standard Tree structure used for various purposes

# Version   Date               Action
# -----------------------------------
# 0.0.1     08 April 2017      Initial issue
# 0.0.2     28 August 2017     Function treeApply() added

# TREE is a list
#' @export
TREE = function(...){
  tree = list(...)
  nms = names(tree)
  if(is.null(nms) & (length(tree) > 0)){
    nms = paste('element', sequence(length(tree)), sep = '.')
  }
  for (i in seq(tree)){
    if(nms[i] == ''){nms[i] = paste('element', i, sep = '.')}
    if (inherits(tree[[i]], 'TREE')){tree[[i]]$name = nms[i]}
  }
  names(tree) <- nms
  class(tree) <- 'TREE'
  return(tree)
}

# Whatever arguments you put in ..., will appear in all the trees(branches)
df2tree.old = function(df, id_cols, var_cols, func = 'mean', name = 'root', ...){
  if(is.empty(df)){return(name)}

  tr = TREE(tree_name = name, ...)
  A  = df %>% aggr(id_cols[1], var_cols, func = func)
  A  = A[!is.na(A[, id_cols[1]]), ] %>% removeRownames %>% column2Rownames(id_cols[1])
  elements = rownames(A)
  pass = length(id_cols) > 1
  if (pass){
    pass = sum(!is.na(df[, id_cols[2]])) != 0
  }
  if (!pass){
    for(e in elements){
      tr[[e]] <- list(e)
      for (vc in var_cols){
        tr[[e]][[vc]] = A[e, vc]
      }
    }
    return(tr)
  }

  for (e in elements){
    tr[[e]] <- df2tree(df = df[df[, id_cols[1]] == e, , drop = F], id_cols = id_cols[-1], var_cols = var_cols, name = e, ...)
    for (vc in var_cols){
      tr[[e]][[vc]] = A[e, vc]
    }
  }
  return(tr)
}



#' @export
df2tree = function(df, id_cols, var_cols, func = 'mean', name = 'root', level = id_cols[1], ...){
  tr = TREE(tree_name = name, level = level, ...)
  if(is.empty(df)){return(tr)}

  A  = df %>% aggr(id_cols[1], var_cols, func = func)
  A  = A[!is.na(A[, id_cols[1]]), ] %>% removeRownames %>% column2Rownames(id_cols[1])
  elements = rownames(A)


  for (e in elements){
    df2 = df[df[,id_cols[1]] == e, ]
    pass = length(id_cols) > 1
    if (pass){
      pass = sum(!is.na(df2[, id_cols[2]])) != 0
    }
    if (!pass){
      tr[[e]] <- list(leaf_name = e)
      for (vc in var_cols){tr[[e]][[vc]] = A[e, vc]}
    } else {
      tr[[e]] <- df2tree(df = df[df[, id_cols[1]] == e, , drop = F], id_cols = id_cols[-1], var_cols = var_cols, func = func, name = e, ...)
      for (vc in var_cols){tr[[e]][[vc]] = A[e, vc]}
    }
  }

  return(tr)
}


# Given function is applied to all branches: trees and leafs
#
#' @export
treeApply = function(tr, func = function(tr) {tr}, ...){
  for(i in seq(tr)){
    if(inherits(tr[[i]], 'TREE')){
      tr[[i]] %<>% func(...)
      tr[[i]] %<>% treeApply(func, ...)}
    else if(inherits(tr[[i]], 'list')){
      tr[[i]] %<>% func(...)
    }}
  return(tr)
}



#ft = function(tr){tr$tree_name %<>% paste0('(', tr$Count, ')');return(tr)}
#fl = function(tr){tr$leaf_name %<>% paste0('(', tr$Count, ')');return(tr)}
#cwt = treeApply(cwt, func_tree = ft, func_leaf = fl)


