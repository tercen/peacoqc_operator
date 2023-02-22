suppressPackageStartupMessages({
  library(tercen)
  library(dplyr, warn.conflicts = FALSE, quietly = TRUE)
  library(dtplyr)
  library(data.table)
  library(flowCore)
  library(PeacoQC)
})

matrix2flowFrame <- function(a_matrix){ 
  minRange <- matrixStats::colMins(a_matrix)
  maxRange <- matrixStats::colMaxs(a_matrix)
  rnge <- maxRange - minRange
  
  df_params <- data.frame(
    name = colnames(a_matrix),
    desc = colnames(a_matrix),
    range = rnge,
    minRange = minRange,
    maxRange = maxRange
  )
  params <- Biobase::AnnotatedDataFrame()
  Biobase::pData(params) <- df_params
  Biobase::varMetadata(params) <- data.frame(
    labelDescription = c("Name of Parameter",
                         "Description of Parameter",
                         "Range of Parameter",
                         "Minimum Parameter Value after Transformation",
                         "Maximum Parameter Value after Transformation")
  )
  flowFrame <- flowCore::flowFrame(a_matrix, params)
  
  return(flowFrame)
}

peacoqc_flowQC <- function(flowframe, input.pars) {
  QC <- PeacoQC(
    flowframe,
    channels = seq(length(colnames(flowframe))-1),
    determine_good_cells  = "all",
    plot = FALSE,
    save_fcs = FALSE,
    output_directory = NULL,
    report = FALSE,
    min_cells = 150,
    max_bins = 500,
    step = 500,
    MAD = input.pars$MAD,
    IT_limit = input.pars$IT_limit,
    consecutive_bins = 5,
    remove_zeros = input.pars$remove_zeros,
    force_IT = 150
  )
  return(QC$GoodCells)
}

ctx <- tercenCtx()

if(ctx$cnames[1] == "filename") {
  filename <- TRUE
  if(ctx$cnames[2] != "Time") {
    stop("Time not detected in the second column.")
  }
} else {
  filename <- FALSE
  if(ctx$cnames[1] != "Time") {
    stop("filename or Time not detected in the top column.")
  } 
}

input.pars <- list(
  MAD = ctx$op.value('MAD', as.double, 6),
  IT_limit = ctx$op.value('IT_limit', as.double, 0.55),
  remove_zeros = ctx$op.value('remove_zeros', as.logical, FALSE)
)

data <- ctx$as.matrix() %>% 
  t() %>% 
  cbind(ctx$cselect())

if(filename == FALSE){
  data$filename <- "singlefile"
}

df <- data.table::as.data.table(data)
df2 <- df[,{
  ff <- matrix2flowFrame(as.matrix(.SD))
  QC_vector <- peacoqc_flowQC(ff, input.pars)
  .(QC=QC_vector)
}, by = filename]

df2 %>% as_tibble() %>% 
  mutate(QC_flag = if_else(QC, "pass", "fail")) %>%
  mutate(.ci = 1:nrow(.) - 1L) %>%
  select(QC_flag, .ci) %>%
  ctx$addNamespace() %>%
  ctx$save()