#analyse_mcmc.R
## Functions to read and format mcmc logs for analysi

#' Read MCMC logs
#'
#' Read mcmc logs into list of posteriors
#'
#' @param path directory in which log files are located
#' @param pattern string to match for files to include in list
#' @param burnin default 0.1
#' @param different_lengths if mcmc chains are different lengths in set to avoid running slowly
#' @importFrom beastio readLog
#' @return logs list of mcmc posteriors
#' @export
read_mcmc_logs <- function(path, pattern, burnin = 0.1, different_lengths = FALSE) {


    log_files <- paste(paste0(path,"/", sep=""),
                   list.files(path = path,
                              pattern = pattern), sep = "")

    logs <- NULL

    if (! different_lengths) {
        logs <- try(beastio::readLog(
        log_files,
        burnin = burnin,
        maxsamples = -1,
        as.mcmc = TRUE,
        burninAsSamples = FALSE
    ))
    }

    if (class(logs) == "try-error" | different_lengths) {
        if (class(logs) == "try-error") {
            print("Caught an error during readLog, trying to read individually....")
        }

        if (different_lengths) {

            print("Reading logs individually....")
        }

        logs <- list()
        for (i in 1:length(log_files)) {

            print(log_files[i])

            mcmc_logs <- beastio::readLog(log_files[i],
                burnin = burnin,
                maxsamples=-1,
                as.mcmc= TRUE,
                burninAsSamples=FALSE)

            logs[[i]] <-  mcmc_logs


        }
        # logs <- purrr::map(log_files, function(file) readLog(file, burnin = burnin,
        #                                                  maxsamples=-1,
        #                                                  as.mcmc= TRUE,
        #                                                  burninAsSamples=FALSE))
        names(logs) <- log_files
    }

    return(logs)
}

#' Calculate state-dependent growth rate difference
#'
#' Calculate mean and HPD intervals for difference between edge and center growth rates
#'
#' @param mcmc_obj  mcmc logs which include birthRateCanonical.1 and birthRateCanonical.0
#' @param dr true simulated death rate, if NULL will look for death rate in deathRateCanonical.1 and deathRateCanonical.0
#' @param fixed_dr  parameter given to BEAST inference, if NULL will look for death rate in deathRateCanonical.1 and deathRateCanonical.0
#' @importFrom HDInterval hdi
#' @return data.frame with mean and HPD intervals
#' @export
calc_state_dependent_growthRateDiff_means <- function(mcmc_obj,
                                     dr,
                                     fixed_dr = NULL) {

    if (is.null(fixed_dr)) {

        if("deathRateCanonical.1" %in% colnames(mcmc_obj)) { #check that death rates are logged if not given

            growthRate_0_posteriors <- mcmc_obj[,"birthRateCanonical.0"] - mcmc_obj[,"deathRateCanonical.0"]
            growthRate_1_posteriors <- mcmc_obj[,"birthRateCanonical.1"] - mcmc_obj[,"deathRateCanonical.1"]
            posteriors <- growthRate_1_posteriors - growthRate_0_posteriors

            } else {

            error("Must include fixed death rate if not logged in mcmc")
            }
    } else {

        growthRate_0_posteriors <- mcmc_obj[,"birthRateCanonical.0"] - fixed_dr
        growthRate_1_posteriors <- mcmc_obj[,"birthRateCanonical.1"] - fixed_dr

        posteriors <- growthRate_1_posteriors - growthRate_0_posteriors

    }


    return(data.frame("parameter" = "growthRateDiff",
                      "mean" = mean(posteriors),
                      "hdi95_lower" = HDInterval::hdi(posteriors,  credMass = 0.95)[1],
                      "hdi95_upper" = HDInterval::hdi(posteriors,  credMass = 0.95)[2],
                      "hdi85_lower" = HDInterval::hdi(posteriors,  credMass = 0.85)[1],
                      "hdi85_upper" = HDInterval::hdi(posteriors,  credMass = 0.85)[2],
                      "hdi75_lower" = HDInterval::hdi(posteriors,  credMass = 0.75)[1],
                      "hdi75_upper" = HDInterval::hdi(posteriors,  credMass = 0.75)[2],
                       "dr" = dr))

}


#' Extract death rate
#'
#' Extract simulated death rate from name of file
#'
#' @param filename
#'
#' @return dr numeric
#' @export
extract_death_rate <- function(filename){

    dr <- as.numeric(unlist(regmatches(filename,
                                       gregexpr("[[:digit:]]+\\.*[[:digit:]][[:digit:]]*",filename))))
    return(dr[-1])
}



