
#Get per day division prob based on lp and ld

get_division_sum <- function(lp, ld) {
    
    alpha = ld/24
    lambda = lp/24
    
    division_sum <- 0
    for (i in 0:11) {
        
        #(1-alpha)*(1-lambda) = P(cell didn't die and it didn't divide)
        add_term <- (((1-alpha)*(1-lambda))^i)*(1-alpha)*lambda
        
        division_sum <- division_sum + add_term 
    }
    
    return(division_sum)
}

get_division_sum(lp = 1, ld =0.0)
get_division_sum(lp = 1, ld =0.86)


get_death_sum <- function(ld) {
    
    division_sum <- 0
    for (i in 0:11) {
        add_term <- ((1 - ld/24)^i)*(ld/24)
        
        division_sum <- division_sum + add_term 
    }
    
    return(division_sum)
}


get_death_sum(ld=0.86)

b_p <- function(death_prob, birth_prob, num_timesteps){
    
    alpha = death_prob
    #(1 - alpha) = probability not dying 1 time step
    #dist = prob of dying at each time step
    d_dist <- (1 - alpha)^(0:(num_timesteps - 1)) * alpha
    
    
    #the complement tells us the probability that the cell enjoys num_timesteps opportunities at birth
    #death_dist <- c(d_dist, 1 - sum(d_dist))
    death_dist <- d_dist
    #death_dist[length(death_dist)] <- 1 - sum(death_dist[1:length(death_dist)-1])
    #death_dist <- 1 - d_dist
    #Now, we want to know - at the point of cell death, what proportion of cells fail to divide at each opportunity that they were able to
    #no_divisions <- sum(death_dist * (1 - birth_prob)^(0:num_timesteps))
    no_divisions_and_died_at_end <- sum(death_dist * (1 - birth_prob)^(0:(num_timesteps-1)))
    
    #but could also never die 
    prob_never_die <- 1 - sum(d_dist)
    
    prob_never_die_or_divide <- prob_never_die*(1-birth_prob)^(num_timesteps-1)
    no_divisions <- no_divisions_and_died_at_end + prob_never_die_or_divide
    #the remainder are cells that managed to divide at least one time before dying
    return(1 - no_divisions)
}


num_timesteps <- 12
num_hours <- 24
ld <- 0.84
lb <- 1 #lp

get_division_sum(lp = lb, ld = ld)
b_r(ld/num_hours, lb/num_hours, num_timesteps)
b_p(ld/num_hours, lb/num_hours, num_timesteps)

b_p <- function(death_prob, birth_prob, num_timesteps){
    
    d_dist <- (1 - death_prob)^(0:(num_timesteps - 1)) * death_prob
    
    #the complement tells us the probability that the cell enjoys num_timesteps opportunities at birth
    death_dist <- c(d_dist, 1 - sum(d_dist))
    
    #Now, we want to know - at the point of cell death, what proportion of cells fail to divide at each opportunity that they were able to
    no_divisions <- sum(death_dist * (1 - birth_prob)^(0:num_timesteps))
    
    #the remainder are cells that managed to divide at least one time before dying
    return(1 - no_divisions)
}


num_timesteps <- 12
num_hours <- 24
ld <- 0.84
lb <- 1 #lp

get_division_sum(lp = lb, ld = ld)
b_r(ld/num_hours, lb/num_hours, num_timesteps)
b_p(ld/num_hours, lb/num_hours, num_timesteps)
