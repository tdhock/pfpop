plot.pfpop_map_l1_verbose <- function(x, ...){
  if(requireNamespace("ggplot2")){
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
      ggplot2::facet_grid(data_i ~ step_i)+
      ggplot2::scale_x_continuous(breaks=seq(0,360,by=90))+
      ggplot2::facet_grid(data_i ~ step_i, scales="free")+
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
        before="grey50",
        after="violet"))+
      ggplot2::geom_point(ggplot2::aes(
        param, param*Linear+Constant),
        size=4,
        shape=21,
        fill=NA,
        data=x$clusters.long)+
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
    fit$clusters <- fread(verbose_file)
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
