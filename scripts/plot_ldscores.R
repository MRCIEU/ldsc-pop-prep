setwd(dirname(rstudioapi::getSourceEditorContext()$path))

load_ldscores <- function(file, gz = TRUE) {
  if(gz) {
    file <- gzfile(file)
  }
  df <- read.table(file, header = T, sep = "\t")
}

for(POP in c("AFR", "AMR", "EAS", "EUR", "SAS", "ALL")) {
  objname <- paste0(POP, "_ldscores_new")
  for(CHR in 1:22) {
    filename <- paste0("../results/ldscores/", POP, "/", CHR, ".l2.ldscore.gz")
    df <- load_ldscores(filename)
    if(CHR == 1) {
      assign(objname, df)
    } else {
      assign(objname, rbind(get(objname), df))
    }
  }
}

for(CHR in 1:22) {
  filename <- paste0("../data/bluepebble/original/EUR/", CHR, ".l2.ldscore")
  df <- load_ldscores(filename, FALSE)
  if(CHR == 1) {
    EUR_ldscores_bp <- df 
  } else {
    EUR_ldscores_bp <- rbind(EUR_ldscores_bp, df)
  }  
}

rolling_windows_mean <- function(i, x, size, step) {
  starts <- seq(1, length(x), step)
  stops <- starts + size
  stops[length(stops)] <- length(x)
  windows <- cbind(starts, stops)
  means <- apply(windows, 1, FUN = function(rw) {
    mean(x[rw[1]:rw[2]], na.rm = T)
  })
  midpoints <- apply(windows, 1, FUN = function(rw) {
    i[round(rw[1] + ((rw[2] - rw[1] )/ 2))]
  })
  return(data.frame(midpoints, means))
}

plot_chr <- function(chr) {
  new <- rolling_windows_mean(
    EUR_ldscores_new[EUR_ldscores_new$CHR == chr,]$BP,
    EUR_ldscores_new[EUR_ldscores_new$CHR == chr,]$L2,
    1000,
    500
  )
  
  bp <- rolling_windows_mean(
    EUR_ldscores_bp[EUR_ldscores_bp$CHR == chr,]$BP,
    EUR_ldscores_bp[EUR_ldscores_bp$CHR == chr,]$L2,
    1000,
    500
  )
  
  plot(
    x = c(new$midpoints, bp$midpoints),
    y = c(new$means, bp$means),
    type = "n",
    xlab = "genomic position (bp)",
    ylab = "L2",
    main = paste("EUR chr", chr)
  )
  lines(
    x = new$midpoints,
    y = new$means,
    col = "black"
  )
  lines(
    x = bp$midpoints,
    y = bp$means,
    col = "red"
  )
  legend(
    "topright",
    legend = c("new", "old"),
    pch = "-",
    col = c("black", "red")
  )
}

plot_chr(1)


plot_chr_pop <- function(chr, pop) {
  new <- rolling_windows_mean(
    get(paste0(pop, "_ldscores_new"))[get(paste0(pop, "_ldscores_new"))$CHR == chr,]$BP,
    get(paste0(pop, "_ldscores_new"))[get(paste0(pop, "_ldscores_new"))$CHR == chr,]$L2,
    1000,
    500
  )
  
  bp <- rolling_windows_mean(
    EUR_ldscores_bp[EUR_ldscores_bp$CHR == chr,]$BP,
    EUR_ldscores_bp[EUR_ldscores_bp$CHR == chr,]$L2,
    1000,
    500
  )
  
  plot(
    x = c(new$midpoints, bp$midpoints),
    y = c(new$means, bp$means),
    type = "n",
    xlab = "genomic position (bp)",
    ylab = "L2",
    main = paste("EUR chr", chr)
  )
  lines(
    x = new$midpoints,
    y = new$means,
    col = "black"
  )
  lines(
    x = bp$midpoints,
    y = bp$means,
    col = "red"
  )
  legend(
    "topright",
    legend = c("new", "old"),
    pch = "-",
    col = c("black", "red")
  )
}

plot_chr_pop(1, "AFR")

