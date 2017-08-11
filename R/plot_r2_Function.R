

#' Visualize standardized effect sizes and model R squared
#'
#' @param x An R2 object from the r2beta function.
#' @param y An R2 object from the r2beta function.
#' @param txtsize The text size of the axis labels.
#' @param maxcov Maximum number of covariates to include in the semi-partial plots.
#' @param r2labs a character vector containing labels for the models. The labels are printed as subscripts on a covariance model matrix.
#' @param cor An argument to be passed to the r2dt function. Only relevant if comparing two R2 objects.
#' @param r2mthd The method used to compute R2
#' @param ... Arguments to be passed to plot
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
#' plot(x=r2)
#' @export

plot.R2 <- function(x,y=NULL,txtsize=10, maxcov=3,
                    r2labs=NULL, r2mthd='sgv', cor=TRUE,
                    ...){

  maxcov= min(maxcov+1, nrow(x))
  if(is.null(r2labs)) r2labs = c('1','2')

  if(toupper(r2mthd)=='SGV'){
    c1.lab=bquote(paste(R[Sigma]^{2} ~~ (hat(Sigma)[.(r2labs[1])])))
    c2.lab=bquote(paste(R[Sigma]^{2} ~~ (hat(Sigma)[.(r2labs[2])])))
    diff.lab = bquote(paste(.(c2.lab)-.(c1.lab)))
  }

  if(toupper(r2mthd)=='KR'){
    c1.lab = bquote(paste(
      R[beta]^{2} * (hat(Sigma)[.(r2labs[1])] ~~ hat(nu)[.(toupper(r2mthd))])))
    c2.lab = bquote(paste(
      R[beta]^{2} * (hat(Sigma)[.(r2labs[2])] ~~ hat(nu)[.(toupper(r2mthd))])))
    diff.lab = bquote(paste(.(c2.lab)-.(c1.lab)))

  }

  if(toupper(r2mthd)=='NSJ'){
    c1.lab = bquote(paste(R[NSJ(m)]^{2} * ( Tr~( hat(Sigma)[.(r2labs[1])] ))))
    c2.lab = bquote(paste(R[NSJ(m)]^{2} * ( Tr~( hat(Sigma)[.(r2labs[2])] ))))
    diff.lab = bquote(paste(.(c2.lab)-.(c1.lab)))
  }


  p = function(xx, lab){

    nrows = min(nrow(xx), maxcov)

    xx$Effect=with(xx, factor(Effect, ordered=T,levels = Effect))
    ggplot2::ggplot(data = xx[1:maxcov,],
                    ggplot2::aes_string(x='Effect',y='Rsq'))+
      ggplot2::geom_point()+
      ggplot2::geom_errorbar(ggplot2::aes_string(
        ymin='lower.CL',ymax='upper.CL'),width=1/5)+
      ggplot2::labs(x = 'Fixed Predictor', y = '') +
      ggplot2::ggtitle(lab) + ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 1/2),
        panel.border = ggplot2::element_rect(colour = 'white'),
        strip.background = ggplot2::element_rect(
          fill = 'white', colour = 'white'),
        strip.text       = ggplot2::element_text(
          face = 'bold.italic', color = 'black', size = txtsize),
        axis.title.y=ggplot2::element_text(
          margin=ggplot2::margin(0,10,0,0)),
        axis.title.x=ggplot2::element_text(
          margin=ggplot2::margin(10,0,0,0))
      )
  }


  if(is.null(y)){

    p(xx=x, lab=c1.lab)

    } else {

    r2b1 = x ; r2b2 = y

    upper_y = round(max(r2b1[1,'upper.CL'], r2b2[1,'upper.CL'])+0.05,1)

    p1 = p(xx=r2b1, lab=c1.lab)+
      ggplot2::scale_y_continuous(limits = c(0, upper_y))
    p2 = p(xx=r2b2, lab=c2.lab)+
      ggplot2::scale_y_continuous(limits = c(0, upper_y))

   tst = r2dt(r2b2, r2b1, cor=cor)

   r2d=lcl=ucl=NULL

   df = data.frame(x=1, r2d = tst$d,
                   lcl = tst$ci[1],
                   ucl = tst$ci[2])

    p3 = ggplot2::ggplot(df, ggplot2::aes(x=x, y = r2d))+
      ggplot2::geom_point()+
      ggplot2::geom_hline(ggplot2::aes(yintercept=0))+
      ggplot2::geom_errorbar(ggplot2::aes(ymin=lcl, ymax=ucl), width=0.05)+
      ggplot2::scale_y_continuous(limits = c(-1,1))+
      ggplot2::scale_x_continuous(limits = c(0.95, 1.05))+
      ggplot2::coord_flip() +
      ggplot2::ylab(diff.lab)+
      ggplot2::xlab('') + ggplot2::theme_bw() + ggplot2::theme(
        axis.title.y= ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y= ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(colour = 'white'))

    grobs = list(p1, p2, p3)
    lay <- rbind(c(1,2), c(3,3))

    g=gridExtra::arrangeGrob(grobs=grobs,
                             layout_matrix=lay,
                             heights=c(.75,.25))

    grid::grid.draw(g)


  }

}




#' Print the contents of an R2 object
#' @param x an object of class R2
#' @param ... other arguments passed to the print function.
#' @export
print.R2 <- function(x, ...){

  # Clean up the output a bit
  nmrc = sapply(x,is.numeric)
  x[,nmrc] = apply(x[,nmrc], 2, round, 3)

  print.data.frame(x[,c('Effect', 'Rsq', 'upper.CL', 'lower.CL')])

}

