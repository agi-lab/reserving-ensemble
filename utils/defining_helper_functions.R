
################################################################################
## This module defines helper functions used throughout the code.
################################################################################

# Helper Functions for Fitting Models
################################################################################
# Takes in array of N_i by occurence years and essentially performs a Vlookup on
# N_i depending on the occurence year in data
# i'th element of N_i corresponds to occurence year i
vlookup_N_i <- function(N, data){
    
    N_rep <- c()
    occurence_years <- as.numeric(data$origin)
    
    for (i in occurence_years) {
        N_rep <- c(N_rep, N[i])
    }
    
    return(N_rep)
}

################################################################################
# Writing the Minorization-Maximization Algorithm
MM_func<-function(w,dat){
    wd<-dat%*%w
    dens<-matrix(NA,nrow=nrow(dat),ncol=ncol(dat))
    for (i in 1:nrow(dat)){
        dens[i,]<-dat[i,]/wd[i]
    }
    mean<-apply(dens,MARGIN=2,FUN=mean)
    return(mean) 
}

MM_optim <- function(w_init,dat,testdat,nIters){
    w <- w_init
    i <- 0
    prev_loss <- -mean(log(dat%*%w_init))
    ## Terminate the algorithm if the iterations exceed the maximum number of iterations
    ##  specified or the loss difference between two subsequent iterations is less than a tolerable error: 
    while (i < nIters) {
        w <- w*MM_func(w, dat)
        current_loss <- -mean(log(dat%*%w))
        if (abs(current_loss - prev_loss) < 1e-16) {
            break
        }
        prev_loss <- current_loss
        i <- i + 1
    }
    return(list(finalparams=w, finalNLL = current_loss))
}

################################################################################
## Transforms the origin and development columns to numerics and factors respectively
################################################################################
df.to.factor <- function(df) {
    df$origin<-as.factor(df$origin)
    df$dev<-as.factor(df$dev)
    df$calendar<-as.numeric(df$calendar)
    
    return(df)
}

df.to.numeric <- function(df) {
    df$origin=as.numeric(as.character(df$origin))
    df$dev=as.numeric(as.character(df$dev))
    df$calendar=as.numeric(as.character(df$calendar))
    
    return(df)
}

################################################################################
# Creates 'object' that has the weights for EW model in the same format as other ensembles
weights_for_EW_all_sims <- function(n.sims, weights){
    w <- c()
    
    for (i in 1:n.sims) {
        w <- c(w, list(list(weights)))
    }
    
    return (w)
}

# Creates 'object' that has the weights for BMV model in the same format as other ensembles
weights_for_BMV_all_sims <- function(n.models, BMV_ind) {
    w <- c()
    
    for (i in BMV_ind) {
        w_i <- rep(0, n.models)
        w_i[i] <- 1
        
        w <- c(w, list(list(w_i)))
    }
    
    return (w)
}

################################################################################
average_by_index <- function(index, lst_of_lst, avg_median=F) {
    # Assume lst_of_lst of unknown length has sublists of the same length as index
    # This function then averages over each index over each list
    
    out_index <- sort(unique(index))
    all_vals <- unlist(lst_of_lst)
    val_index <- sapply(names(all_vals), FUN = function(x) as.numeric(unlist(strsplit(x, "-"))[1]))
    
    out_val <- c()
    for (i in out_index) {
        avg <- ifelse(avg_median, median(all_vals[val_index == i]), mean(all_vals[val_index == i]))
        out_val <- c(out_val, avg)
    }
    
    return (out_val)
}

average_by_index_sim <- function(index, lst_of_mat, avg_median=F) {
    # Assume lst_of_mat of unknown length has matrices of nrow = length(index)
    # and unknown length
    # This function then averages over each index over each matrix
    
    out_index <- sort(unique(index))
    model_index <- rep(rep(1:ncol(lst_of_mat[[1]]), each = length(index)), length(lst_of_mat))
    val_index <- rep(index, times = (ncol(lst_of_mat[[1]]) * length(lst_of_mat)))
    all_vals <- unlist(lst_of_mat)
    
    model_vals <- list()
    for (j in 1:ncol(lst_of_mat[[1]])) {
        out_val <- c()
        for (i in out_index) {
            avg <- ifelse(avg_median, median(all_vals[(val_index == i)&(model_index == j)]), mean(all_vals[(val_index == i)&(model_index == j)]))
            out_val <- c(out_val, avg)
        }
        model_vals[[j]] <- out_val
    }
    
    return (model_vals)
}