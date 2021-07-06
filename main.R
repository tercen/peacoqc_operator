library(tercen)
library(dplyr, quietly = TRUE)
library(tibble)
library(flowCore)
library(magrittr)
library(PeacoQC)

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

peacoqc_flowQC <- function(flowframe, input.pars){
  QC <- try(PeacoQC(flowframe,
                    channels = seq(length(colnames(flowframe))-1),
                    determine_good_cells  = "all",
                    plot = FALSE,
                    save_fcs = FALSE,
                    output_directory = NULL,
                    #name_directory = ,
                    report = FALSE,
                    #events_per_bin = FindEventsPerBin(remove_zeros, ff, channels,min_cells, max_bins, step),
                    min_cells = 150,
                    max_bins = 500,
                    step = 500,
                    MAD = input.pars$MAD,
                    IT_limit = input.pars$IT_limit,
                    consecutive_bins = 5,
                    remove_zeros = input.pars$remove_zeros,
                    #suffix_fcs = "_QC",
                    force_IT = 150), silent = TRUE)
  return(QC$GoodCells)
}

ctx <- tercenCtx()

if(ctx$cnames[1] != "Time")stop("Time is not in the top column.")
celldf <- ctx %>% dplyr::select(.ri, .ci) 
if(nrow(celldf) != length(table(celldf)))stop("There are multiple values in one of the cells.")

input.pars <- list(
  MAD = ifelse(is.null(ctx$op.value('MAD')), 6, as.double(ctx$op.value('MAD'))),
  IT_limit = ifelse(is.null(ctx$op.value('IT_limit')),  0.55, as.double(ctx$op.value('IT_limit'))),
  remove_zeros = ifelse((ctx$op.value('remove_zeros') == "false"), FALSE, TRUE)
)

data <- ctx$as.matrix() %>% t() %>% cbind((ctx$cselect(ctx$cnames[[1]]))) %>% 
  as.matrix() %>% matrix2flowFrame()

qc_df <- data.frame(matrix(ncol=0, nrow=nrow(data)))
qc_df$QC_flag <- ifelse(peacoqc_flowQC(data, input.pars) == TRUE, "pass", "fail")
peacoqc_QC <- cbind(qc_df, .ci = (0:(nrow(qc_df)-1)))
ctx$addNamespace(peacoqc_QC) %>% ctx$save()


