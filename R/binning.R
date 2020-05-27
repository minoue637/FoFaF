binning <- function(max_val, bin_num = 100) {
    start_deg <- 0
    if (bin_num <= max_val - start_deg + 1) {
      if (1 == bin_num) {
        base <- max_val - start_deg + 1
        interval_length <- max_val - start_deg + 1
      }
      else {
        is.warn <- options()$warn
        options(warn = -1)
        ff <- function(x) {
          max_val - start_deg + 1 - sum(floor(x^(0:(bin_num - 
                                                      1))))
        }
        base <- uniroot(ff, interval = c(1 + 1e-15, max_val - 
                                           start_deg + bin_num + 1.1), tol = .Machine$double.eps)$root
        options(warn = is.warn)
        interval_length <- floor(base^(0:(bin_num - 1)))
      }
    }
    else if ( (0 == bin_num) || (bin_num > max_val - start_deg + 1)) {
      bin_num         <- max_val - start_deg + 1
      interval_length <- rep(1, bin_num)
      base            <- 1
    }
    
    bin_vector <- rep(bin_num + start_deg - 1, max_val + 1)
    begin_deg <- c(start_deg, start_deg + cumsum(interval_length)[-bin_num])
    end_deg <- begin_deg + interval_length - 1
    if (start_deg > 0) 
      bin_vector[1:start_deg] <- 0:(start_deg - 1)
    for (i in 1:bin_num) bin_vector[(begin_deg[i]:end_deg[i]) + 1] <- i + start_deg - 1
    
    names(bin_vector) <- 0:(length(bin_vector) - 1) 
    if (start_deg  > 1) {
      center_bin        <- c(0:(start_deg - 1),sqrt(begin_deg * end_deg))
    } else center_bin <- sqrt(begin_deg * end_deg)
    bin_num <- max(bin_vector) + 1
    return(list(bin = bin_vector, center_bin = center_bin, 
                start = begin_deg, end = end_deg, bin_num = bin_num))
  }