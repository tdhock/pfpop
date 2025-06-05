plot.pfpop_map_l1_verbose <- function(x, ...){
  if(requireNamespace("ggplot2")){
    size.values <- c(
      first=4.5,
      opt=3,
      last=1.5)
    color.values <- c(
      first="red",
      opt="blue",
      last="orange")
    ggplot2::ggplot()+
      ggplot2::theme_bw()+
      ggplot2::geom_vline(ggplot2::aes(
        xintercept=x),
        color="grey",
        data=x$vline)+
      ggplot2::geom_segment(ggplot2::aes(
        min_param, min_param*Linear+Constant,
        xend=max_param, yend=max_param*Linear+Constant),
        size=2,
        color="grey50",
        data=x$list_result$model)+
      ggplot2::facet_grid(data_i ~ step_i, scales="free", labeller=ggplot2::label_both)+
      ggplot2::geom_rect(ggplot2::aes(
        xmin=start, xmax=end,
        fill=opt,
        ymin=-Inf, ymax=Inf),
        data=x$rect,
        alpha=0.5,
        color="black")+
      ggplot2::geom_label(ggplot2::aes(
        param, ifelse(Linear_diff<0, -Inf, Inf),
        vjust=ifelse(Linear_diff<0, 0, 1),
        label=Linear_diff),
        data=x$breaks,
        alpha=0.5)+
      ggplot2::scale_fill_manual(values=c(
        before="violet",
        after="grey50"))+
      ggplot2::scale_size_manual(
        breaks=names(size.values),
        values=size.values)+
      ggplot2::scale_color_manual(
        breaks=names(color.values),
        values=color.values)+
      ggplot2::geom_abline(ggplot2::aes(
        slope=Linear, intercept=Constant,
        color=ptr),
        alpha=0.3,
        size=2,
        data=x$clusters.long)+
      ggplot2::geom_point(ggplot2::aes(
        param, param*Linear+Constant,
        color=ptr,
        size=ptr),
        data=x$clusters.long[order(-size.values[ptr])])+
      ggplot2::geom_abline(ggplot2::aes(
        slope=Linear, intercept=Constant),
        data=x$clusters.long)
  }
}
pfpop_map_l1_verbose <- function(degrees_vec, penalty, weight_vec=rep(1,length(degrees_vec))){
  fit <- pfpop_map_l1(degrees_vec, penalty, weight_vec, tempfile())
  fit$list_result <- pfpop_list_l1(degrees_vec, penalty, weight_vec, tempfile())
  fit$clusters.long <- melt(
    fit$clusters,
    measure.vars=measure(
      ptr, value.name,
      pattern="(first|opt|last)_(param|Constant|Linear)"))
  cldt <- function(opt, start, end)data.table(opt, start, end)
  fit$rect <- fit$clusters[, rbind(
    cldt("before", first_param, opt_param),
    cldt("after", opt_param, last_param)),
    by=.(data_i,step_i)]
  fit$vline <- unique(fit$list_result$model[
  , .SD[, .(x=unique(c(min_param,max_param)))]
  , by=data_i
  ])
  class(fit) <- c("pfpop_map_l1_verbose", class(fit))
  fit
}
pfpop_map_l1 <- function(degrees_vec, penalty, weight_vec=rep(1,length(degrees_vec)), verbose_file=""){
  fit <- pfpop_map_l1_interface(degrees_vec, penalty, weight_vec, verbose_file)
  if(verbose_file!=""){
    fit$clusters <- fread(verbose_file, colClasses=list(numeric=c(
      "first_param", "opt_param", "last_param",
      "first_Linear", "opt_Linear", "last_Linear",
      "first_Constant", "opt_Constant", "last_Constant")))
    breaks_file <- paste0(verbose_file,"_breaks")
    fit$breaks <- fread(breaks_file)
  }
  fit
}
pfpop_list_l1 <- function(degrees_vec, penalty, weight_vec=rep(1,length(degrees_vec)), verbose_file=""){
  fit <- pfpop_list_l1_interface(degrees_vec, penalty, weight_vec, verbose_file)
  if(verbose_file!=""){
    fit$model <- fread(verbose_file)
  }
  fit
}
