

#' Generate partial contrast matrices
#'
#' @param r2ob An r2 object from the r2beta function.
#' @param txtsize The text size of the axis labels.
#' @return A visual representation of the model and semi-partial R squared
#' from the r2 object provided.
#' @examples
#' library(nlme)
#' library(r2glmm)
#'
#' data(Orthodont)
#'
#' # Linear mixed model
#' lmemod = lme(distance ~ age*Sex, random = ~1|Subject, data = Orthodont)
#'
#' r2 = r2beta(model=lmemod,partial=TRUE,method='sgv')
#'
#' plot_r2(r2ob=r2)
#' @export plot_r2

plot_r2 <- function(r2ob, txtsize = 11){

  r2ob$Effect=with(r2ob, factor(Effect, ordered=T,levels = Effect))

  ggplot2::ggplot(data = r2ob, ggplot2::aes_string(x = 'Effect', y = 'Rsq'))+
    ggplot2::geom_point()+
    ggplot2::geom_errorbar(
      ggplot2::aes_string(ymin='lower.CL',ymax='upper.CL'),width=1/5)+
    ggplot2::labs(x = 'Fixed Predictor', y = 'R Squared') +
    ggplot2::theme_bw() +
    ggplot2::theme(
          panel.border = ggplot2::element_rect(colour = 'white'),
          axis.title   = ggplot2::element_text(face = "bold.italic",
                                    color ="black", size = txtsize),
          strip.background = ggplot2::element_rect(
                      fill = 'white', colour = 'white'),
          strip.text       = ggplot2::element_text(
                      face = 'bold.italic', color = 'black', size = txtsize)
          )

}



