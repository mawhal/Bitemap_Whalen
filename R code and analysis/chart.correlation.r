function (R, histogram = TRUE, method = c("pearson", "kendall", 
                                          "spearman"),...){
  x = R
  # x = checkData(R, method = "matrix")
  if (missing(method)) 
    method = method[1]
  cormeth <- method
  panel.smooth <- function (x, y, col = "slateblue", bg = NA, pch = 19, 
                            cex = 3, col.smooth = "red", span = 2/3, iter = 3,...) 
  {
    smoothScatter(x, y, add=T, nrpoints=0, pch = pch, col = scales::alpha(col,0.25), bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
      lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth,lwd=2, ...)
  }
  panel.cor <- function(x, y, digits = 2, prefix = "", 
                        use = "pairwise.complete.obs", method = cormeth, 
                        cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = use, method = method)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor)) 
      cex <- 0.8/strwidth(txt)
    test <- cor.test(as.numeric(x), as.numeric(y), method = method)
    # Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
    #                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
    #                                                                           "**", "*", ".", " "))
    text(0.5, 0.5, txt, cex = cex )
    # text(0.8, 0.8, Signif, cex = cex, col = 2)
  }
  f <- function(t) {
    dnorm(t, mean = mean(x), sd = sd.xts(x))
  }
  dotargs <- list(...)
  dotargs$method <- NULL
  rm(method)
  hist.panel = function(x, ... = NULL) {
    par(new = TRUE)
    hist(x, col = "whitesmoke", probability = TRUE, 
         axes = FALSE, main = "", breaks = "FD")
    lines(density(x, na.rm = TRUE), col = "black")
    rug(x)
  }
  if (histogram) 
    pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, 
          diag.panel = hist.panel)
  else pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, 
             cex.labels  = 1)
}
