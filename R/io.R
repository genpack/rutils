

# Header
# Filename:      io.R
# Description:   This library, is a part of package viser and is used for data import and export to/from other environments
# Author:        Nicolas Berta
# Email :        nicolas.berta@gmail.com
# Start Date:    28 October 2016
# Last Revision: 17 July 2019
# Version:       2.0.1

# Version History:

# Version   Date               Action
# ----------------------------------
# 1.0.0     28 October 2016    Initial Issue: Started with function readODBC().
# 1.1.0     21 November 2016   Added filters to function readODBC().
# 1.2.0     01 December 2016   Function sqlFilter() generates script for filtering. Functions sqlFilterComponent.*() removed.
# 1.2.1     01 December 2016   Function sqlScript() generates complete sql query.
# 1.2.2     01 December 2016   Function readODBC() calls sqlScript() first to generate the query and then executes the query.
# 1.2.3     24 January 2017    Function readODBC() transfers additional arguments passed to function sqlQuery() in package RODBC.
# 1.3.0     17 May 2017        Function rScriptFilter() added, converts a filter into a R script of a condition
# 1.3.1     01 August 2017     Function rScriptFilter() exported.
# 1.3.2     15 May 2018        Function sqlScript() modified: changes column names with 'AS' if argument 'fields' is a named vector
# 1.3.3     12 September 2018  Function runSQL() added
# 2.0.0     26 June 2019       spark and athena tools and AGGERGATOR module added.
# 2.0.1     17 July 2019       Function sqlFilter() modified: small bug fixed.


################## ODBC TOOLS: ###############################

# This function reads a table by ODBC and generates an object of type TIME.SERIES
# Reads a table using ODBC
# similar functions will be like this:
# readODBC.TIME.SERIES
# readODBC.nibeTree
#
# Examples:
#
# D = readODBC('ASET_VALU_CURR', fields = c('ASET_I', 'ASET_VALU_D', 'ASET_VALU_A'), dbName = 'PVDATA')

# script.sql.select = function(fields = '*'){}
# script.sql.from   = function(dbName,tableName){}
# script.sql.where  =

valid.filter.elements = c('domain', 'min', 'max', 'type', 'na.rm', 'equal', 'query')
valid.field.types     = c('nominal', 'numeric', 'date', 'time')

#' @export
sqlFilter = function(colName, filter, vf = T){
  if(vf){
    verify(filter, 'list', names_domain = valid.filter.elements, varname = 'filter')
    filter$type  <- verify(filter$type, 'character', domain = valid.field.types, varname = 'filter$type', default = 'nominal')
    filter$na.rm <- verify(filter$na.rm, 'logical', domain = c(T,F), varname = 'filter$na.rm', default = F)
    filter$equal <- verify(filter$equal, 'logical', domain = c(T,F), varname = 'filter$equal', default = F)
    filter$query <- verify(filter$query, 'character', varname = 'filter$query')
  }
  if(!is.null(filter$query)){
    return(paste0(colName, " IN (", filter$query, ")"))
  }
  if (filter$equal){
    minopr = " >= "
    maxopr = " <= "
  } else {
    minopr = " > "
    maxopr = " < "
  }
  scr   = ""
  joint = ""
  switch(filter$type,
         'nominal' = {
           if(!is.null(filter$domain)){
             if (!inherits(filter$domain, 'character')){filter$domain %<>% as.character}
             if(length(filter$domain) == 1){
               scr = paste0(colName, " = '",  filter$domain, "'")
             } else {
               scr = paste0(colName, " IN ('",  paste(filter$domain, collapse = "','"), "')")
             }
             joint = " AND "
           }
           if(filter$na.rm){
             scr = paste0(scr, joint, colName, " IS NOT NULL")
           }
         },
         'numeric' = {
           if(!is.null(filter$domain)){
             if (!inherits(filter$domain, 'numeric')){filter$domain %<>% as.numeric}
             if(length(filter$domain) == 1){
               scr = paste0(colName, " = ",  filter$domain)
             } else {
               scr = paste0(colName, " IN (",  paste(filter$domain, collapse = ", "), ")")
             }
             joint = " AND "
           }
           if (!is.null(filter$min)){
             scr = paste0(colName, minopr, filter$min)
             joint = " AND "
           }
           if (!is.null(filter$max)){
             scr = paste0(scr, joint, colName, maxopr, filter$max)
             joint = " AND "
           }
           if(filter$na.rm){scr = paste0(scr, joint, colName, " IS NOT NULL")}
         },
         'date'    = {
           if (!is.null(filter$min)){
             if(vf){verify(filter$min, valid.time.classes, varname = 'filter$min')}
             if (inherits(filter$min, 'character')){minDate = char2Date(filter$min)} else {minDate = as.Date(filter$min)}
             minDate = paste0("CAST(", "'", as.character(minDate), "'", " AS DATE)")
             scr = paste0(colName, minopr, minDate)
             joint = " AND "
           }
           
           if (!is.null(filter$max)){
             verify(filter$max, valid.time.classes, varname = 'filter$max')
             if (inherits(filter$max, 'character')){maxDate = char2Date(filter$max)} else {maxDate = as.POSIXlt(filter$max)}
             maxDate = paste0("CAST(", "'", as.character(maxDate), "'", " AS DATE)")
             scr = paste0(scr, joint, colName, maxopr, maxDate)
             
             # maxDate = as.character(maxDate)
             # scr = paste0(scr, joint, colName, maxopr, "'", maxDate, "'")
             joint = " AND "
           }
           if(filter$na.rm){scr = paste0(scr, joint, colName, " IS NOT NULL")}
         },
         'time'    = {
           if (!is.null(filter$min)){
             verify(filter$min, valid.time.classes, varname = 'filter$min')
             # if (inherits(filter$min, 'character')){minTime = char2Time(filter$min)} else {minTime = as.POSIXlt(filter$min)}
             minTime = as.time(filter$min, target_class = 'POSIXlt')
             minTime = as.character(minTime, format = '%Y-%m-%d %H:%M:%S')
             minTime = paste0('CAST(', "'", minTime, "'", " AS TIMESTAMP)")
             scr = paste0(colName, minopr, minTime)
             joint = " AND "
           }
           
           if (!is.null(filter$max)){
             verify(filter$max, valid.time.classes, varname = 'filter$max')
             if (inherits(filter$max, 'character')){maxTime = char2Time(filter$max)} else {maxTime = as.POSIXlt(filter$max)}
             maxTime = as.character(maxTime, format = '%Y-%m-%d %H:%M:%S')
             maxTime = paste0('CAST(', "'", maxTime, "'", " AS TIMESTAMP)")
             scr = paste0(scr, joint, colName, maxopr, maxTime)
             joint = " AND "
           }
           if(filter$na.rm){scr = paste0(scr, joint, colName, " IS NOT NULL")}
         })
  
  return(scr)
}

#' @export
rScriptFilter = function(colName, filter, vf = T){
  if(vf){
    verify(filter, 'list', names_domain = valid.filter.elements, varname = 'filter')
    filter$type  <- verify(filter$type, 'character', domain = valid.field.types, varname = 'filter$type', default = 'nominal')
    filter$na.rm <- verify(filter$na.rm, 'logical', domain = c(T,F), varname = 'filter$na.rm', default = F)
    filter$equal <- verify(filter$equal, 'logical', domain = c(T,F), varname = 'filter$equal', default = F)
    filter$query <- verify(filter$query, 'character', varname = 'filter$query')
  }
  
  if (filter$equal){
    minopr = " >= "
    maxopr = " <= "
  } else {
    minopr = " > "
    maxopr = " < "
  }
  scr   = ""
  joint = ""
  switch(filter$type,
         'nominal' = {
           if(!is.null(filter$domain)){
             if (!inherits(filter$domain, 'character')){filter$domain %<>% as.character}
             scr = paste0(colName, " %in% c('",  paste(filter$domain, collapse = "','"), "')")
             joint = " & "
           }
           if(filter$na.rm){scr = paste0(scr, joint, "!is.na(", colName, ")")}
         },
         'numeric' = {
           if (!is.null(filter$min)){
             scr = paste0(colName, minopr, filter$min)
             joint = " & "
           }
           if (!is.null(filter$max)){
             scr = paste0(scr, joint, colName, maxopr, filter$max)
             joint = " & "
           }
           if(filter$na.rm){scr = paste0(scr, joint, "!is.na(", colName, ")")}
         },
         'date'    = {
           if (!is.null(filter$min)){
             if(vf){verify(filter$min, valid.time.classes, varname = 'filter$min')}
             if (inherits(filter$min, 'character')){minDate = char2Date(filter$min)} else {minDate = as.Date(filter$min)}
             minDate = as.character(minDate)
             scr = paste0(colName, minopr, "'", minDate, "'")
             joint = " & "
           }
           
           if (!is.null(filter$max)){
             verify(filter$max, valid.time.classes, varname = 'filter$max')
             if (inherits(filter$max, 'character')){maxDate = char2Date(filter$max)} else {maxDate = as.POSIXlt(filter$max)}
             maxDate = as.character(maxDate)
             scr = paste0(scr, joint, colName, maxopr, "('", maxDate, "' %>% as.Date)")
             joint = " & "
           }
           if(filter$na.rm){scr = paste0(scr, joint, "!is.na(", colName, ")")}
         },
         'time'    = {
           if (!is.null(filter$min)){
             verify(filter$min, valid.time.classes, varname = 'filter$min')
             # if (inherits(filter$min, 'character')){minTime = char2Time(filter$min)} else {minTime = as.POSIXlt(filter$min)}
             minTime = as.time(filter$min, target_class = 'POSIXlt')
             minTime = as.character(minTime, format = '%Y-%m-%d %H:%M:%S')
             scr = paste0(colName, minopr, "('", minTime, "' %>% as.POSIXlt)")
             joint = " & "
           }
           
           if (!is.null(filter$max)){
             verify(filter$max, valid.time.classes, varname = 'filter$max')
             if (inherits(filter$max, 'character')){maxTime = char2Time(filter$max)} else {maxTime = as.POSIXlt(filter$max)}
             maxTime = as.character(maxTime, format = '%Y-%m-%d %H:%M:%S')
             scr = paste0(scr, joint, colName, maxopr,"('", maxTime, "' %>% as.POSIXlt)")
             joint = " AND "
           }
           if(filter$na.rm){scr = paste0(scr, joint, "!is.na(", colName, ")")}
         })
  
  return(scr)
}


#' @export
sqlScript = function(tableName, fields = NULL, dbName = 'UDRBSCMS', filter = NULL, vf = T){
  if(vf){
    tableName  = verify(tableName, 'character', lengths = 1, varname = 'tableName')
    fields     = verify(fields, 'character', varname = 'fields', default = '*')
    filter     = verify(filter, 'list', varname = 'filter')
  }
  
  fldnms = names(fields)
  if(is.null(fldnms)){fields.str = paste(fields, collapse = " , ")} else {
    assert(length(fldnms) == length(fields))
    fldnms[fldnms == ''] <- fields[fldnms == '']
    fields.str = paste(paste(fields, 'AS', fldnms), collapse = " , ")
  }
  
  query      = paste0("SELECT ", fields.str, " from ", dbName, ".", tableName)
  
  if (!is.null(filter)){
    first = T
    query = query %++% " WHERE "
    for (fn in names(filter)){
      if (!first){query = query %++% " AND "}
      qfn   = sqlFilter(fn, filter[[fn]], vf = vf)
      query = query %++% qfn
      if (qfn != ""){first = F}
    }
  }
  return(query)
}


#' @export
readODBC <- function(tableName, fields = NULL, dbName = 'UDRBSCMS', filter = NULL, dsn = 'Teradata_Prod', vf = T){
  # Verifications:
  if(vf){assert(require(RODBC), "Package 'RODBC' is not installed!", err_src = match.call()[[1]])}
  
  channel    = odbcConnect(dsn = dsn, ...)
  D          = sqlQuery(channel = channel, query = sqlScript(tableName = tableName, fields = fields, dbName = dbName, filter = filter, vf = vf), ...)
  close(channel)
  return(D)
}

#' @export
runSQL <- function(query, ...){
  channel    = odbcConnect(...)
  D          = sqlQuery(channel = channel, query = query)
  close(channel)
  return(D)
}


################## AGGREGATOR MODULE: ###############################

AGGREGATOR = setRefClass(
  "AGGREGATOR", 
  
  fields = c(tables = "list"), 
  
  methods = c(
    initialize = function(data){
      tables <<- list(data = data)
    },
    
    # groups by given id columns and summarises sum of values but saves previous results in the tables list to do it faster for 
    # next runs
    # property 'tables' is a list of tables containing the original dataset and all pre-computed aggregated dataset
    # if aggregated data does not exist in the list, the function will compute it and adds the table to the list
    aggregate_by = function(id_cols, value_col){
      ds = tables$data
      mincols = ncol(ds) - 1
      for(tbl in tables){
        ns = colnames(tbl) %-% value_col
        if((ns %<% id_cols) & (id_cols %<% ns)){return(tbl)}
        nc = length(ns)
        if(id_cols %<% ns){
          if(mincols > nc){
            ds = tbl
            mincols = nc
          }
        }
      }
      scr = paste0("ds ", "%>% group_by(", id_cols %>% paste(collapse = ","), ") %>% summarise(", value_col, " = ", "sum", "(", value_col, ")) %>% ungroup()")
      out = parse(text = scr) %>% eval
      tables[[length(tables) + 1]] <<- out 
      return(out)
    },
    
    get.flowNetwork = function(id_cols, value_col, percentage = F){
      links = NULL
      
      for(i in sequence(length(id_cols) - 1)){
        aggregate_by(c(id_cols[i], id_cols[i + 1]), value_col) %>% 
          select_(source = id_cols[i], target = id_cols[i + 1], value = value_col) %>% 
          mutate(svname = id_cols[i], tvname = id_cols[i + 1]) %>% rbind(links) -> links
      }
      
      links %<>% mutate(hovertext = paste0(source, ' --> ', target, ': ', value))
      
      links$source = paste(links$svname, links$source, sep = "=")
      links$target = paste(links$tvname, links$target, sep = "=")
      
      links %<>% left_join(links %>% group_by(source) %>% summarise(sumval = sum(value)), by = 'source') %>% 
        mutate(ratio = round(100*value/sumval, digits = 2)) %>% 
        mutate(hovertext = hovertext %>% paste0(' (', ratio, '%)')) 
      #links$tooltip = paste()
      if(percentage){
        links %<>% left_join(links %>% group_by(target) %>% summarise(sumratio = sum(ratio)) %>% select(source = target, sumratio), by = 'source') %>% 
          mutate(sumratio = ifelse(is.na(sumratio), 100, sumratio)) %>% 
          mutate(pathratio = round(ratio*sumratio/100, digits = 2))
      }
      
      nodes = data.frame(id = c(links$source, links$target)) %>% 
        distinct(id, .keep_all = T) %>% mutate(label = id)
      
      list(nodes = nodes, links = links)
    },
    
    plot.sankey = function(id_cols, value_col, percentage = F, plotter = 'networkD3'){
      if(length(id_cols) == 0) return(NULL)
      
      get.flowNetwork(id_cols = id_cols, value_col = value_col, percentage = percentage) %>% 
        viser::viserPlot(key = 'id', linkTooltip = 'hovertext', label = 'label', linkWidth = ifelse(percentage, 'pathratio', 'value'), source = 'source', target = 'target', plotter = plotter, type = 'sankey')
    }
  )
)

################## SPARK TOOLS ###############################

#' @export
spark.read_s3 = function(path){
  sc <- sparklyr::spark_connect(master = 'local')
  el <- sparklyr::spark_read_csv(sc, name = 'el', path = path)
  return(el)
}


# bins numeric columns to categorical 
# binners is a list containing:
# input_col = colname to be binned, must have numeric values
# output_col = colname to be generated as a result of binning
# splits = named numeric vector containing splits and labels for them
#' @export
spark.bin = function(df, binners){
  if(inherits(df, c('data.frame', 'tibble'))){
    for (b in binners){
      df %>% pull(b$input_col) %>% cut(breaks = b$splits %>% unname, labels = names(b$splits)) -> df[, b$output_col]
    }
  } else if (inherits(df, 'tbl_spark')){
    
    df %<>% sparklyr::ft_bucketizer(
      input_cols   = binners %>% lapply(function(x) x$input_col) %>% unlist %>% as.character,
      output_cols  = binners %>% lapply(function(x) x$output_col) %>% unlist %>% as.character,
      splits_array = binners %>% lapply(function(x) x$splits %>% unname))
    
    for(b in binners){
      df %<>% rename(indexed = b$output_col) %>% sparklyr::ft_index_to_string(input_col = 'indexed', output_col = b$output_col, labels = b$splits %>% names) %>% 
        select(- indexed)
    }
    
  } else stop("Unsupported class for argument 'df'")
  return(df)
}


# 2btransferred to gener or niraspark (sparker)?
#' @export
spark.rename = function(tbl, ...){
  ns = c(...) %>% verify('character')
  lb = names(ns)
  if(is.null(lb)){return(tbl)}
  scr = paste0("tbl %>% dplyr::rename(", lb %>% paste(ns, sep = ' = ') %>% paste(collapse = ' , '), ")")
  parse(text = scr) %>% eval
}

#' @export
spark.mutate = function(tbl, ...){
  ns = c(...) %>% verify('character')
  lb = names(ns)
  if(is.null(lb)){return(tbl)}
  scr = paste0("tbl %>% dplyr::mutate(", lb %>% paste(ns, sep = ' = ') %>% paste(collapse = ' , '), ")")
  parse(text = scr) %>% eval
}

#' @export
spark.select = function(tbl, ...){
  ns = c(...) %>% verify('character') %>% intersect(colnames(tbl))
  if(length(ns) == 0){return(tbl)}
  scr = paste0("tbl %>% dplyr::select(", ns %>% paste(collapse = ' , '), ")")
  parse(text = scr) %>% eval
}

#' @export
spark.unselect = function(tbl, ...){
  ns = c(...) %>% verify('character') %>% intersect(colnames(tbl))
  if(length(ns) == 0){return(tbl)}
  scr = paste0("tbl %>% dplyr::select(", paste('-',ns) %>% paste(collapse = ' , '), ")")
  parse(text = scr) %>% eval
}

#' @export
spark.dummify = function(tbl, ..., keep_original = T, keep_index = F){
  cat_cols = c(...) %>% verify('character', domain = colnames(tbl), default = character()) %>% unique
  if(length(cat_cols) == 0) {return(tbl)}
  
  # scr     <- paste0('tbl %>% filter(', paste(paste(cat_cols, '>= 0'), collapse = ' ,'), ')')
  scr     <- paste0('tbl %>% mutate(', paste0(paste0(cat_cols, ' = as.character(abs(as.integer(', cat_cols, ')))'), collapse = ' ,'), ')')
  tbl     <- parse(text = scr) %>% eval
  uniques <- list()
  for(col in cat_cols){
    icn <- col %>% paste('Index', sep = '_') # icn: indexed column name
    tbl %<>% ft_string_indexer(input_col = col, output_col = icn, string_order_type = 'alphabetDesc')
    
    scr  <- paste0("tbl %>% distinct(", col, ", ", icn, ") %>% collect")
    parse(text = scr) %>% eval -> uniques[[col]]
    if(nrow(uniques[[col]]) > 2){
      labels <- paste(col, uniques[[col]] %>% pull(col), sep= '_')
      # uniques[[col]] <- uniques[[col]][- length(uniques[[col]])]
      tbl %<>% ft_one_hot_encoder(icn, 'outlist', drop_last = F) %>% 
        sdf_separate_column('outlist', into = labels) %>% select(-outlist)
      if(!keep_original){
        tbl <- parse(text = paste0("tbl %>% select(-", col, ")")) %>% eval
      }
      if(!keep_index){
        tbl <- parse(text = paste0("tbl %>% select(-", icn, ")")) %>% eval
      }
    }
  }
  return(tbl) 
}

################## ATHENA SQL TOOLS ###############################

#' @export
athena.buildConnection = function(bucket = "s3://aws-athena-query-results-192395310368-ap-southeast-2/", region = 'ap-southeast-2'){
  pyathena = reticulate::import('pyathena')
  pyathena$connect(s3_staging_dir = bucket, region_name = region)
}

#' @export
athena.read_s3   = function(con, query, max_rows = 100000){
  pandas   = reticulate::import('pandas')
  nrow     = pandas$read_sql(query %>% sql.nrow, con)[1,1]
  if(nrow > max_rows) stop('The query will result a table with' %>% paste(nrow, 'rows which is higher than maximum:', max_rows))
  else pandas$read_sql(query, con)
}

#' @export
athena.buildDataset = function(conn, dsName = 'dataset'){
  curser = conn$cursor()
  qry = paste("CREATE DATABASE", dsName)
  cursor$execute(qry)
}

# Example:
# acon = buildAthenaConnection()
# buildAthenaTable(acon, dsName = 'event', tblName = 'eventlogs', s3_data_path = datapath, drop_if_exists = T)

#' @export
sql.mlMapper = function(input, caseid_col, ts_col, et_col, var_col, value_col, tscat_col, variables){
  cvpr = paste0("select ", caseid_col, ", ", tscat_col, ", ", var_col, ", max_by(", value_col, ", ", ts_col, ") as latestValue from (", input, ") group by ", caseid_col, ", ", tscat_col, ", ", var_col)
  cvpr = sql.cast(cvpr, id_col = c(caseid_col, tscat_col), var_col = var_col, value_col = 'latestValue', variables = variables, aggregator = 'SUM')
  # cptime = paste0("select ", caseid_col, "\n")
  # cptime %<>% paste0(",MIN(", ts_col, ") AS caseStartTime", "\n")
  # cptime %<>% paste0(",MAX(", ts_col, ") AS latestEventTime", "\n")
  # cptime %<>% paste0(",MAX(IF(", et_col, " = 'LoanClosed', ", ts_col, ", cast('1900-01-01' as TIMESTAMP))) as closureTime", "\n")
  # cptime %<>% paste0(" from (", input, ") group by ", caseid_col)
  # sql.leftjoin(cvpr, cptime, by = caseid_col)
  cvpr
}

#' @export
sql.group_by = function(input, ...){
  by    = c(...)
  query = paste0("(", input, ") GROUP BY ", by %>% paste(collapse = ","))
  list(query = query, id_cols = by)
}

#' @export
sql.summarise = function(input, ...){
  if(inherits(input, 'list')){
    id_cols = input$id_cols
    input   = input$query
  } else {id_cols = ''}
  vars = list(...)
  if(length(vars) == 1 & inherits(vars[[1]], 'list')) vars = vars[[1]]
  vnms = names(vars)
  qry  = "SELECT "
  add  = id_cols %>% paste(collapse = ', ')
  if(add != ""){qry %<>% paste0(add, ", ")}
  
  for (vn in vnms){
    qry %<>% paste0(vars[[vn]]," AS ", vn)
    if(vn != vnms[length(vnms)]){
      qry  %<>% paste0(", \n")
    }
  }
  
  qry %>% paste0(" FROM ", input)
}

#' @export
sql.cast = function(input, id_col, var_col, value_col, variables, aggregator = 'SUM'){
  qry = paste0("SELECT ",  id_col %>% paste(collapse = ','),  ",", "\n")
  for (var in variables){
    qry %<>% paste0(aggregator, "(IF(variable = '", var, "',", value_col, ", NULL)) AS ", var, ", \n")
  }
  qry %>% paste0("count(", value_col, ") as Count from (", input, ")", " group by ", id_col %>% paste(collapse = ','))
}

# Example:
# athenaTableCast('event.eventlogs', 'caseid', 'variable', 'value', variables = c('OriginationChannel', 'Income'))

#' @export
sql.filter = function(input, ...){
  filter = list(...)
  if(length(filter) == 1 & inherits(filter[[1]], 'list')) filter = filter[[1]]
  fnms = names(filter)
  
  qry = paste0("SELECT * FROM (", input, ") WHERE ")
  for (fn in fnms){
    if(inherits(filter[[fn]], 'character')){
      qry %<>% paste0(fn," = '", filter[[fn]], "'")
    } else {
      qry %<>% paste0(fn," = ", filter[[fn]])
    }
    if(fn != fnms[length(fnms)]){qry %>% paste0(" AND ")}
  }
  return(qry)
}

# Returns the case profile containing the latest values of each variable
# input must be in eventlog format:
#' @export
sql.caseProfile = function(input, caseid_col, ts_col, et_col, var_col, value_col, variables, with_times = T){
  cvpr = paste0("select ", caseid_col, ", ", var_col, ", max_by(", value_col, ", ", ts_col, ") as latestValue from (", input, ") group by ", caseid_col, ", ", var_col)
  cprf = sql.cast(cvpr, id_col = caseid_col, var_col = var_col, value_col = 'latestValue', variables = variables, aggregator = 'SUM')
  if(with_times){
    cptime = paste0("select ", caseid_col, "\n")
    cptime %<>% paste0(",MIN(", ts_col, ") AS caseStartTime", "\n")
    cptime %<>% paste0(",MAX(", ts_col, ") AS latestEventTime", "\n")
    cptime %<>% paste0(",MAX(IF(", et_col, " = 'LoanClosed', ", ts_col, ", cast('1900-01-01' as TIMESTAMP))) as closureTime", "\n")
    # cptime %<>% paste0(",closureTime - caseStartTime AS LoanAge", "\n")
    cptime %<>% paste0(" from (", input, ") group by ", caseid_col)
    return(sql.leftjoin(cprf, cptime, by = caseid_col))
  } else {return(cprf)}
}


#' @export
sql.arrange = function(input, by){
  paste0(input, " ORDER BY ", by %>% paste(collapse = ","))
}

#' @export
sql.leftjoin = function(table1, table2, by){
  paste0("SELECT * FROM (", table1, ") a LEFT JOIN (", table2, ") b ON a.", by, " = b.", by)
}

#' @export
sql.select = function(input, ...){
  arg = list(...)
  if(length(arg) == 1 & inherits(arg[[1]], 'list')) arg = arg[[1]]
  anms = names(arg)
  
  qry = paste0("SELECT ")
  
  for (i in sequence(length(arg))){
    if(is.empty(anms[i])){
      qry %<>% paste0(arg[[i]])
    } else {
      qry %<>% paste0(arg[[i]], " AS ", anms[i])
    }
      
    if(i != length(arg)){qry %<>% paste0(", ")}
  }
  qry %>% paste0(" FROM (", input, ") ")
}

#' @export
sql.mutate = function(input, ...){
  arg = list(...)
  if(length(arg) == 1 & inherits(arg[[1]], 'list')) arg = arg[[1]]
  anms = names(arg)
  
  qry = paste0("SELECT *, ")
  
  for (a in anms){
    qry %<>% paste0(arg[[a]], " AS ", a)
    if(a != anms[length(anms)]){qry %<>% paste0(", ")}
  }
  qry %>% paste0(" FROM (", input, ") ")
}

# Example:
# qry = 'event.eventlogs' %>% sql.mutate(fvalue = "CAST(value as DOUBLE)")
# qry %<>% paste0(" where variable <> 'variable'")
# read_s3.athena(acon, query = qry %>% paste0(" limit 20")) %>% View

# qry = 'event.eventlogs' %>% sql.mutate(fvalue = "CAST(value as DOUBLE)")
# qry %<>% paste0(" where variable <> 'variable'")
# qry %<>% sql.mutate(status = "CASE WHEN fvalue > 1 THEN 'A' WHEN fvalue = 1 THEN 'C' ELSE 'B' END")
# qry = "SELECT * FROM event.eventlogs WHERE variable = 'ProductCode'"


# Example:
# sql.binColumn('LoanTenure', breaks = c('< 1 Yr' = 365, '1-2 Yrs' = 730, '2-3 Yrs' = 1095, '3-4 Yrs' = 1460, '4-5 Yrs' = 1825, '> 6 Yrs' = 2190))
#' @export
sql.binColumn = function(column, breaks){
  breaks %<>% sort
  mutscr = "CASE \n"
  nms = names(breaks)
  for (i in sequence(length(breaks))){
    if(i == 1){
      mutscr %<>% paste0("WHEN ", column, " < ", breaks[i], " THEN '", nms[i], "' \n")
    } else if (i < length(breaks)){
      mutscr %<>% paste0("WHEN ", column, " > ", breaks[i - 1], " AND ", column, " < ", breaks[i], " THEN '", nms[i], "' \n")
    } else {
      mutscr %<>% paste0("ELSE '", nms[i], "' END")
    }
  }  
  
  # CASE
  # WHEN Tenure < 0 THEN '0'
  # WHEN Tenure < 365 AND Tenure > 0 THEN '1'
  # WHEN Tenure < 730 AND Tenure > 365 THEN '2'
  # WHEN Tenure < 1095 AND Tenure > 730 THEN '3'
  # WHEN Tenure < 1460 AND Tenure > 1095 THEN '4'
  # WHEN Tenure < 1825 AND Tenure > 1460 THEN '5'
  # WHEN Tenure < 2190 AND Tenure > 1825 THEN '6'
  # ELSE '7'
  # END
  return(mutscr)
}

#' @export
sql.head = function(input, n = 6){
  'SELECT * from (' %>% paste(input, ') LIMIT', n)
}



#' @export
sql.nrow = function(input){
  "SELECT COUNT(*) FROM (" %>% paste0(input, ")")
}

#' @export
sql.cast = function(input, id_col, var_col, value_col, variables, aggregator = 'SUM'){
  qry = paste0("SELECT ",  id_col %>% paste(collapse = ','),  ",", "\n")
  for (var in variables){
    qry %<>% paste0(aggregator, "(IF(variable = '", var, "',", value_col, ", NULL)) AS ", var, ", \n")
  }
  qry %>% paste0("count(", value_col, ") as Count from (", input, ")", " group by ", id_col %>% paste(collapse = ','))
}
