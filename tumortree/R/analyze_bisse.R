#' Read log files
#'
#' Import BiSSE log files to get estimated rates 
#' 
#' @param file log file
#' @param percent.burnin fraction of mcmc burnin
#' @importFrom tidyr pivot_longer
#' @return data.frame
#' @export

read_log_file <- function(file, percent.burnin = 0.1) {
    
    #read in .log file and pivot to parameter and state variables
    if (file.exists(file)) {
        

        posteriors <- read.table(file, sep = "\t", header = TRUE, fill = TRUE) %>% 
            tidyr::pivot_longer(
                cols = diversification.1.:transition_rates.2.,
                names_to = c("parameter", "state"), 
                names_sep = "[.]",
                values_to = "value")
        #remove burnin iterations
        burnin <- round(0.1 * nrow(posteriors))
        posteriors <- posteriors[burnin:nrow(posteriors),]
                
        return(posteriors)

    } else {
        warning("File does not exist, returning NA...")
        return(NA)
    }
    
}