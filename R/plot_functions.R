plot_data_sims <- function(x, decomposition, ...){

  if(decomposition){

    if(isFALSE(x@composite))
      stop(strwrap(prefix = " ", initial = "",
                   "Cannot plot decomposition for non-composite processes"))

    if(!("yG" %in% names(x@test)))
      stop(strwrap(prefix = " ", initial = "",
                   "Cannot plot decomposition for a composite process without
                   the global, local, and error processes decomposed."))

    toReturn <- suppressWarnings(plotDataSimsYGL(x, ...))
  }else{
    toReturn <- suppressWarnings(plotDataSimsY(x, ...))
  }
  return(toReturn)
}

plotDataSimsYGL <- function(x, ...){
  d <- ncol(x@training$x)

  titleStat <- ifelse(x@stationary, "Stationary", "Non-Stationary")
  titleComp <- ifelse(x@composite, "Composite", "Non-Composite")
  plotTitle <- paste(titleStat, titleComp, "BCGP", sep = " ")

  if(d == 1){
    dataToPlotPred <- data.frame(x = x@test$x, y = x@test$y, yG = x@test$yG,
                                 yL = x@test$yL) %>%
      tidyr::gather("process", "value", -x)

    dataToPlot <- data.frame(x = x@training$x, process = "data",
                             value = x@training$y) %>%
      dplyr::bind_rows(dataToPlotPred)

    dataPlot <- ggplot2::ggplot(mapping = ggplot2::aes(x = x, y = value,
                                                       color = process)) +
      ggplot2::ggtitle(plotTitle) +
      ggplot2::geom_point(data = dplyr::filter(dataToPlot, process == "data"),
                          ggplot2::aes(color = "red")) +
      ggplot2::geom_line(data = dplyr::filter(dataToPlot, process == "y"),
                         ggplot2::aes(color = "black")) +
      ggplot2::geom_line(data = dplyr::filter(dataToPlot, process == "yG"),
                         ggplot2::aes(color = "blue")) +
      ggplot2::geom_line(data = dplyr::filter(dataToPlot, process == "yL"),
                         ggplot2::aes(color = "green")) +
      ggplot2::theme_classic() +
      ggplot2::scale_color_manual(name = NULL,
                                  values = c("red" = "red", "black" = "black",
                                             "blue" = "blue",
                                             "green" = "green"),
                                  labels = c("Observed Data", "True Process",
                                             "Global Process", "Local Process"),
                                  breaks = c("red", "black", "blue", "green"),
                                  guide = ggplot2::guide_legend(
                                    override.aes = list(
                                      linetype = c("blank", "solid", "solid",
                                                   "solid"),
                                      shape = c(16, NA, NA, NA)))) +
      ggplot2::ylab("Y(x)") +
      ggplot2::theme(legend.position = "bottom",
                     plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                     axis.title = ggplot2::element_text(size = 16))
  }else if(d == 2){
    message(strwrap(prefix = " ", initial = "",
                    "The decomposition into global and local processes is not
                 plotted for 2-D data. The plot will be the overall process."))

    dataPlot <- plotDataSimsY(x, ...)

  }else{
    stop("Plotting is currently only supported for 1-D and 2-D data.")
  }
  return(dataPlot)
}

plotDataSimsY <- function(x, ...){
  d <- ncol(x@training$x)
  dots <- list(...)

  titleStat <- ifelse(x@stationary, "Stationary", "Non-Stationary")
  titleComp <- ifelse(x@composite, "Composite", "Non-Composite")
  plotTitle <- paste(titleStat, titleComp, "BCGP", sep = " ")

  if(d == 1){
    dataToPlot <- dplyr::bind_rows(data.frame(x = x@training$x,
                                              process = "data",
                                              value = x@training$y),
                                   data.frame(x = x@test$x, process = "y",
                                              value = x@test$y))

    dataPlot <- ggplot2::ggplot(mapping = ggplot2::aes(x = x, y = value,
                                                       color = process)) +
      ggplot2::ggtitle(plotTitle) +
      ggplot2::geom_point(data = dplyr::filter(dataToPlot, process == "data"),
                          ggplot2::aes(color = "red")) +
      ggplot2::geom_line(data = dplyr::filter(dataToPlot, process == "y"),
                         ggplot2::aes(color = "black")) +
      ggplot2::theme_classic() +
      ggplot2::scale_color_manual(name = NULL,
                                  values = c("red" = "red", "black" = "black"),
                                  labels = c("Observed Data", "True Process"),
                                  breaks = c("red", "black"),
                                  guide = ggplot2::guide_legend(
                                    override.aes = list(
                                      linetype = c("blank", "solid"),
                                      shape = c(16, NA)))) +
      ggplot2::ylab("Y(x)") +
      ggplot2::theme(legend.position = "bottom",
                     plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                     axis.title = ggplot2::element_text(size = 16))
  }else if(d == 2){

    dataToPlotPred <- data.frame(x1 = x@test$x[,1], x2 = x@test$x[,2],
                                 y = x@test$y, type = "test")
    dataToPlot <- data.frame(x1 = x@training$x[,1], x2 = x@training$x[,2],
                             y = x@training$y, type = "training") %>%
      dplyr::bind_rows(dataToPlotPred)

    if(isTRUE(x@test$grid) && isTRUE(dots$raster)){
      dataPlot <- ggplot2::ggplot(data = dplyr::filter(dataToPlot,
                                                       type == "test"),
                                  mapping = ggplot2::aes(x = x1, y = x2,
                                                         fill = y)) +
        ggplot2::ggtitle(plotTitle) +
        ggplot2::theme_classic() +
        ggplot2::geom_raster(interpolate = TRUE) +
        ggplot2::scale_fill_gradientn(name = "Y(x)",
                                      colors = c("red", "yellow", "blue")) +
        ggplot2::theme(legend.position = "bottom",
                       plot.title = ggplot2::element_text(size = 16,
                                                          hjust = 0.5),
                       axis.title = ggplot2::element_text(size = 16)) +
        ggplot2::labs(x = expression(x[1]), y = expression(x[2])) +
        ggplot2::geom_point(data = dplyr::filter(dataToPlot,
                                                 type == "training"),
                            mapping = ggplot2::aes(x = x1, y = x2, size = type),
                            inherit.aes = FALSE) +
        ggplot2::scale_size_manual(name = NULL,
                                   values = c("training" = 1.5),
                                   labels = c("Observed Data"),
                                   breaks = c("training"))

    }else{
      dataPlot <- ggplot2::ggplot(data = dataToPlot,
                                  mapping = ggplot2::aes(x = x1, y = x2,
                                                         color = y)) +
        ggplot2::ggtitle(plotTitle) +
        ggplot2::theme_classic() +
        ggplot2::geom_point(mapping = ggplot2::aes(size = type)) +
        ggplot2::scale_size_manual(name = NULL,
                                   values = c("test" = 1, "training" = 4),
                                   labels = c("True Process", "Observed Data"),
                                   breaks = c("test", "training"),
                                   guide = ggplot2::guide_legend(
                                     override.aes = list(
                                       size = c(1, 4),
                                       shape = c(16, 16)))) +
        ggplot2::scale_color_gradient2(name="Y(x)", mid = "yellow", low = "red",
                                       high = "blue", midpoint = 0) +
        ggplot2::theme(legend.position = "bottom",
                       plot.title = ggplot2::element_text(size = 16,
                                                          hjust = 0.5),
                       axis.title = ggplot2::element_text(size = 16)) +
        ggplot2::labs(x = expression(x[1]), y = expression(x[2]))
    }


  }else{
    stop("Plotting is currently only supported for 1-D and 2-D data.")
  }
  return(dataPlot)
}

plot_var_sims <- function(x, ...){

  dots <- list(...)
  d <- ncol(x@training$x)
  if(d == 1){
    predV <- data.frame(x = x@test$x, y = x@parameters$VTest)

    varPlot <- ggplot2::ggplot(mapping = ggplot2::aes(x, y)) +
      ggplot2::ggtitle("Non-Stationary BCGP:\nVariance Process") +
      ggplot2::geom_line(data = predV, color = "black") +
      ggplot2::theme_classic() +
      ggplot2::ylab(~ paste(sigma ^ 2, "(x)")) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
                     axis.title = ggplot2::element_text(size = 16))
  }else if(d == 2){
    predV <- data.frame(x1 = x@test$x[, 1], x2 = x@test$x[, 2],
                        y = x@parameters$VTest)

    if(isTRUE(x@test$grid) && isTRUE(dots$raster)){
      varPlot <- ggplot2::ggplot(data = predV,
                                 mapping = ggplot2::aes(x = x1, y = x2,
                                                        fill = y)) +
        ggplot2::ggtitle("Non-Stationary BCGP:\nVariance Process") +
        ggplot2::theme_classic() +
        ggplot2::geom_raster(interpolate = TRUE) +
        ggplot2::scale_fill_gradientn(name = expression(paste(sigma ^ 2, "(x)")),
                                      colors = c("red", "yellow", "blue")) +
        ggplot2::theme(legend.position = "bottom",
                       plot.title = ggplot2::element_text(size = 16,
                                                          hjust = 0.5),
                       axis.title = ggplot2::element_text(size = 16)) +
        ggplot2::labs(x = expression(x[1]), y = expression(x[2]))

    }else{
      varPlot <- ggplot2::ggplot(data = predV,
                                 mapping = ggplot2::aes(x = x1, y = x2,
                                                        color = y)) +
        ggplot2::ggtitle("Non-Stationary BCGP:\nVariance Process") +
        ggplot2::theme_classic() +
        ggplot2::geom_point() +
        ggplot2::scale_color_gradient2(name = expression(paste(sigma ^ 2, "(x)")),
                                       mid = "yellow", low = "red",
                                       high = "blue",
                                       midpoint = mean(predV$y)) +
        ggplot2::theme(legend.position = "bottom",
                       plot.title = ggplot2::element_text(size = 16,
                                                          hjust = 0.5),
                       axis.title = ggplot2::element_text(size = 16)) +
        ggplot2::labs(x = expression(x[1]), y = expression(x[2]))
    }

  }else{
    stop("Plotting is currently only supported for 1-D and 2-D data.")
  }
  return(varPlot)
}
