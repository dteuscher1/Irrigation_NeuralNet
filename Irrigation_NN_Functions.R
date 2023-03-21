library(mvtnorm)
library(truncnorm)
library(latex2exp)
source("AMCMC/AMCMCUpdate.R")
source("ClustFunctions/ClusteringRFunctions.R")
# Functions for updating parameters

update_beta <- function(xs, smat, zs, m, n) {
    beta_cov_mat <- solve((t(xs) %*% xs)/1 + solve(smat)) 
    beta_mean <- beta_cov_mat %*% ((t(xs)%*%zs)/1 + solve(smat) %*% m)
    betas <- beta_mean + t(chol(beta_cov_mat))%*%rnorm(length(beta_mean))
    return(betas)
}

activation <- function(input, fx="ident"){
    if(fx%in%c("relu", "tanh", "sigmoid", "ident")){
        output <- switch(fx,
                         "relu" = pmax(matrix(0, nrow=nrow(input), ncol=ncol(input)), input),
                         "tanh" = tanh(input),
                         "sigmoid" = plogis(input),
                         "ident" = input)
        output <- 
            return(output)
    } else {
        stop(paste(fx, 'not a valid activation function'))
    }
}
pdp_plot <- function(var_name, data){
    ag <- readxl::read_excel("../Data/MultiVariable_withZones.xlsx") %>% 
        mutate(Aspect.x = cos(AspectDeg_),
               Aspect.y = sin(AspectDeg_)) %>% data.frame()
    var_vector <- seq(min(ag[,var_name]), max(ag[,var_name]), length.out = 100)
    explan_vars <- c("ELEVATION", "kryield19", "SlopeDEg_m", "Aspect.x",
                     "Aspect.y", "Kr19NDVI_m", "krslopedeg", "TWI", "Kr18NDVI_m")
    value <- which(var_name %in% explan_vars)
    lower <- apply(avg_zs[,,value], 2, quantile, .025)
    avg <- apply(avg_zs[,,value], 2, median)
    upper <- apply(avg_zs[,,value], 2, quantile, .975)
    graph <- data.frame(values= var_vector, lower = lower, est = avg, upper = upper)
    p <- ggplot(graph, aes(x = values)) + 
        geom_line(aes(y = est)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), fill = "dodgerblue", alpha = .3) +
        geom_line(aes(y = lower), lty = 2, color = "red") + 
        geom_line(aes(y = upper), lty = 2, color = "red") + 
        labs(x = var_name, y = TeX("PDP(X_p)")) + 
        theme_minimal()
    p
}

# Function to calculate the mode; Used on the chains of predictions to 
# determine the classified zone
getmode <- function(x) {
    uniqv <- unique(x)
    uniqv[which.max(tabulate(match(x, uniqv)))]
}

# Function that calculates the adjusted RAND index
adj.rand.index <- function(c1,c2) {
    n <- length(c1)
    if ( length(c2) != n ) stop("Clustering must be the same length.")
    t1 <- table(c1)
    t2 <- table(c2)
    t12 <- table(c1,c2)
    expected <- sum(choose(t1,2)) * sum(choose(t2,2)) / choose(n,2)
    numerator <- sum(choose(t12,2)) - expected
    denominator <- 0.5*(sum(choose(t1,2))+sum(choose(t2,2))) - expected
    numerator / denominator
} 
# Function to fit the model
bayes_nn_model <- function(nzones, nneurons, act.fns, spatial_basis = 0, ndraws = 1000, nburn = 100, nthin = 1, feature_importance = FALSE, permutations = 50, pdp = FALSE){  
    ## Model settings
    nzones <- nzones #integer in {2, 3, 4, ...}
    nneurons <- nneurons #Both nlayers and nneurons NOT including input and output layer
    act.fns <- act.fns #any of c("relu", "tanh", "sigmoid", "ident")
    ndraws <- ndraws
    nburn <- nburn
    amcmc.it <- 250
    thin <- nthin
    amcmc.eps <- 1e-7
    
    ## Inferred model settings
    tot.it <- nburn+ndraws*thin
    kp.seq <- seq(nburn+thin, tot.it, by=thin)
    kp <- 0
    nlayers <- length(nneurons)
    spatial_basis <- spatial_basis
    
    ###############################
    ## Data Reading and Cleaning ##
    ###############################
    
    ## Read in the data
    ag <- readxl::read_excel("../Data/MultiVariable_withZones.xlsx")
    if(spatial_basis == 0){
        spatial <- read_csv("../BasisFunctions/SpatialBasisFunctions.csv", col_types = cols())[,0]
    } else{
        spatial <- read_csv("../BasisFunctions/SpatialBasisFunctions.csv", col_types = cols())[,1:spatial_basis]
    }
    # Select the necessary variables and create an average VWC variable
    # Create an x and y part of the the aspect degree
    ag <- ag %>% dplyr::select(ObjectId, X, Y, ELEVATION, kryield19, SlopeDEg_m, 
                               AspectDeg_, AspectClas, krVWC8_29, krVWC5_31, krVWC6_25a, 
                               krVWC4_29, krslopedeg, Kr19NDVI_m, Kr18NDVI_m, TWI, VWCchange) %>%
        mutate(AvgVWC = (krVWC8_29 + krVWC5_31 + krVWC6_25a + krVWC4_29)/4,
               Aspect.x = cos(AspectDeg_),
               Aspect.y = sin(AspectDeg_))
    
    # Split the AvgVWC variable into three groups: Low, Medium, and High
    ag <- ag %>% mutate(Water = cut(AvgVWC, 
                                    breaks = quantile(AvgVWC, probs=seq(0, 1, length=nzones+1)),
                                    labels = 1:nzones,
                                    include.lowest=TRUE)) %>% 
        dplyr::select(-AvgVWC)
    
    prediction_matrix <- model.matrix(Water~.-1,
                                      data=ag %>% 
                                          dplyr::select(ELEVATION, kryield19, SlopeDEg_m, Aspect.x, Aspect.y,
                                                        Kr19NDVI_m, krslopedeg, TWI, Water, Kr18NDVI_m))
    prediction_matrix <- apply(prediction_matrix, 2, scales::rescale, to=c(-1, 1))
    preds_spatial <- spatial
    
    
    ag <- read.csv("../Data/ag_train.csv")
    if(spatial_basis == 0){
        spatial <- read_csv("../Data/spatial_train.csv", col_types = cols())[,0]
    } else{
        spatial <- read_csv("../Data/spatial_train.csv", col_types = cols())[,1:spatial_basis]
    }
    
    ag <- ag %>% dplyr::select(ObjectId, X, Y, ELEVATION, kryield19, SlopeDEg_m, 
                               AspectDeg_, AspectClas, krVWC8_29, krVWC5_31, krVWC6_25a, 
                               krVWC4_29, krslopedeg, Kr19NDVI_m, Kr18NDVI_m, TWI, VWCchange) %>%
        mutate(AvgVWC = (krVWC8_29 + krVWC5_31 + krVWC6_25a + krVWC4_29)/4,
               Aspect.x = cos(AspectDeg_),
               Aspect.y = sin(AspectDeg_))
    
    # Split the AvgVWC variable into three groups: Low, Medium, and High
    ag <- ag %>% mutate(Water = cut(AvgVWC, 
                                    breaks = quantile(AvgVWC, probs=seq(0, 1, length=nzones+1)),
                                    labels = 1:nzones,
                                    include.lowest=TRUE)) %>% 
        dplyr::select(-AvgVWC)
    
    model_matrix <- model.matrix(Water~.-1,
                                 data=ag %>% 
                                     dplyr::select(ELEVATION, kryield19, SlopeDEg_m, Aspect.x, Aspect.y,
                                                   Kr19NDVI_m, krslopedeg, TWI, Water, Kr18NDVI_m))
    input_scaled <-apply(model_matrix, 2, scales::rescale, to=c(-1, 1))
    
    ####################################
    ## Functions to update parameters ##
    ####################################
    
    update_beta <- function(xs, smat, zs, m, n) {
        beta_cov_mat <- solve((t(xs) %*% xs)/1 + solve(smat)) 
        beta_mean <- beta_cov_mat %*% ((t(xs)%*%zs)/1 + solve(smat) %*% m)
        betas <- beta_mean + t(chol(beta_cov_mat))%*%rnorm(length(beta_mean))
        return(betas)
    }
    
    activation <- function(input, fx="ident"){
        if(fx%in%c("relu", "tanh", "sigmoid", "ident")){
            output <- switch(fx,
                             "relu" = pmax(matrix(0, nrow=nrow(input), ncol=ncol(input)), input),
                             "tanh" = tanh(input),
                             "sigmoid" = plogis(input),
                             "ident" = input)
            output <- 
                return(output)
        } else {
            stop(paste(fx, 'not a valid activation function'))
        }
    }
    
    ##################################
    ## Initial Values of Parameters ##
    ##################################
    
    ## Neural Network Parameters
    nn <- vector("list", length=nlayers)
    for(l in 1:nlayers){
        if(l==1){
            ## Define the input matrix w/no intercept
            nn[[l]]$input <-input_scaled
            
            ## Set Random weights
            nn[[l]]$wgts <- matrix(rnorm(ncol(nn[[l]]$input)*nneurons[l], 0, 0.1), 
                                   nrow=ncol(nn[[l]]$input), ncol=nneurons[l])
        } else {
            nn[[l]]$input <- cbind(1,activation(nn[[l-1]]$input%*%nn[[l-1]]$wgts, fx=act.fns[l-1]))
            nn[[l]]$wgts <- matrix(rnorm(ncol(nn[[l]]$input)*nneurons[l], 0, 0.1), 
                                   nrow=ncol(nn[[l]]$input), ncol=nneurons[l])
        }
        
    }
    
    X_nn <- cbind(1,activation(nn[[nlayers]]$input%*%nn[[nlayers]]$wgts, fx=act.fns[nlayers]))
    X <- as.matrix(cbind(X_nn, spatial))
    beta <- matrix(c(-qnorm(mean(ag$Water==1)), rep(0, ncol(X)-1)), ncol=1)
    
    # Cut points
    if(nzones==2){
        cuts <- c(-Inf, 0, Inf)
    } else {
        cuts <- c(-Inf, 0, cumsum(prop.table(table(ag$Water)))[-1] %>% qnorm(., mean=beta[1], sd=1))
        est.cuts <- 3:(length(cuts)-1)
        nc <- length(est.cuts)
        if(nc==1){ 
            trans.cuts <- log(cuts[est.cuts])
        } else {
            trans.cuts <- log(c(cuts[est.cuts[1]], diff(cuts[est.cuts])))
        }
        cut.amcmc <- list(mn=matrix(0, nrow=nc, ncol=1),
                          var=matrix(0, nrow=nc, ncol=nc))
    }
    
    
    #####################################
    ## Matrices for Holding MCMC draws ##
    #####################################
    
    beta.draws <- matrix(0, nrow=ndraws, ncol=length(beta))
    if(nzones>2){
        cut.draws <- matrix(0, nrow=ndraws, ncol=nc)
    }
    weight.draws <- vector("list", length = nlayers)
    weight.amcmc <- vector("list", length = nlayers)
    
    for(i in 1:nlayers){
        weight.draws[[i]] <- array(0,c(nrow(nn[[i]]$wgts), ncol(nn[[i]]$wgts), tot.it)) # Dimensions of this 
        weight.amcmc[[i]] <- vector("list", length = ncol(nn[[i]]$wgts))
        weight.amcmc[[i]] <- list(mn=matrix(0,nrow=ncol(nn[[i]]$input) * nneurons[i],ncol=1),
                                  var=matrix(0,nrow=ncol(nn[[i]]$input)*nneurons[i],ncol=ncol(nn[[i]]$input)*nneurons[i]))
    }
    
    predictions <- matrix(0, nrow = ndraws, ncol = nrow(prediction_matrix))
    ############################
    ## Run the MCMC Algorithm ##
    ############################
    pb <- txtProgressBar(min = 0, max = tot.it, style = 3)
    for(it in 1:tot.it){
        
        ## Sample the latent Gaussian variables
        z <- sapply(1:nrow(ag), FUN=function(ind){
            lvl <- as.numeric(ag$Water[ind])
            lwr <- cuts[lvl]
            upr <- cuts[lvl+1]
            mn <- sum(X[ind,]*beta)
            return(rtruncnorm(1, mn, sd=1, a=lwr, b=upr))
        })
        
        ## Sample the beta coefficients
        beta <- update_beta(X, smat = 100* diag(ncol(X)), 
                            zs = z, m = rep(0, ncol(X)), n = nrow(X))
        beta_nn <- beta[1:ncol(X_nn),]
        beta_spatial <- beta[-(1:ncol(X_nn)),]
        ## Sample cut points
        if(nzones>2) {
            
            ## Propose new cuts points by proposing new transformed cut points
            if(it > amcmc.it){
                prop.var <- (2.4^2/nc)*(cut.amcmc$var+amcmc.eps*diag(nc))
            } else {
                prop.var <- amcmc.eps*diag(nzones-2)
            }
            prop.cuts <- cuts
            prop.trans <- trans.cuts + t(chol(prop.var))%*%rnorm(nc) #MVN Draw
            prop.cuts[est.cuts] <- cumsum(exp(prop.trans))
            
            ## Evaluate the likelihood under proposed cutpoints
            llike.diff <- sapply(1:nrow(ag), function(ind){
                lvl <- as.numeric(ag$Water[ind])
                mn <- X[ind,]%*%beta
                
                #Proposed
                lwr <- prop.cuts[lvl]
                upr <- prop.cuts[lvl+1]
                prop.prob <- log(pnorm(upr, mn, 1)-pnorm(lwr, mn, 1))
                
                #Current
                lwr <- cuts[lvl]
                upr <- cuts[lvl+1]
                cur.prob <- log(pnorm(upr, mn, 1)-pnorm(lwr, mn, 1))
                
                return(prop.prob-cur.prob)
            }) %>% sum()
            
            ## Evaluate the Metropolis-Hasting Ratio
            MH <- llike.diff + sum(dnorm(prop.trans, 0, 10, log=TRUE) -
                                       dnorm(trans.cuts, 0, 10, log=TRUE))
            
            ## Determine whether to accept or reject
            if(log(runif(1)) < MH){
                cuts <- prop.cuts
                trans.cuts <- prop.trans
            }
            
            ## Update the AMCMC matrix
            cut.amcmc <- AMCMC.update(trans.cuts, cut.amcmc$mn, cut.amcmc$var, it)
            
        } #End if(nzones>2)
        
        ## Update NN Weights
        for(l in 1:nlayers){
            nn_prop <- nn
            if(it > amcmc.it){
                prop.var <- (2.4^2/nrow(weight.amcmc[[l]]$var))*(weight.amcmc[[l]]$var+amcmc.eps*diag(nrow(weight.amcmc[[l]]$var))) # Is the is number of parameters??
            } else {
                prop.var <- .00001 * diag(nrow(nn[[l]]$wgts) * ncol(nn[[l]]$wgts)) # Not sure here
            }
            nn_prop[[l]]$wgts <- matrix(rmvnorm(1, as.vector(nn[[l]]$wgts), prop.var), nrow = nrow(nn[[l]]$wgts), ncol = ncol(nn[[l]]$wgts)) #sigma of some sort using AMCMC
            for(L in l:nlayers){
                if(L==1){
                    ## Define the input matrix w/no intercept
                    nn_prop[[L]]$input <- input_scaled
                } else {
                    nn_prop[[L]]$input <- cbind(1,activation(nn_prop[[L-1]]$input%*%nn_prop[[L-1]]$wgts, fx=act.fns[L-1]))
                }
            }
            X_prop <- cbind(1,activation(nn_prop[[nlayers]]$input%*%nn_prop[[nlayers]]$wgts, fx=act.fns[nlayers]))
            # Add spatial part in here
            MH <- sum(dnorm(z, cbind(X_prop, X[,-(1:ncol(X_nn))])  %*% beta, 1, log = TRUE) - dnorm(z, X %*% beta, 1, log = TRUE)) + 
                sum(dnorm(nn_prop[[l]]$wgts, 0, 1, log = TRUE) - dnorm(nn[[l]]$wgts, 0, 1, log = TRUE))
            if(log(runif(1)) < MH){
                X <- cbind(X_prop, X[,-(1:ncol(X_nn))])
                nn <- nn_prop
            }
            weight.amcmc[[l]] <- AMCMC.update(as.vector(nn[[l]]$wgts), weight.amcmc[[l]]$mn, weight.amcmc[[l]]$var, it)
            weight.draws[[l]][,,it] <- nn[[l]]$wgts
        }
        
        X_part <- cbind(1,activation(nn[[nlayers]]$input%*%nn[[nlayers]]$wgts, fx=act.fns[nlayers]))
        X <- as.matrix(cbind(X_part, spatial))
        
        for(l in 1:nlayers){
            if(l == 1){
                new_input <- prediction_matrix
            } else{
                new_input <- cbind(1, activation(new_input %*% nn[[l-1]]$wgts, fx = act.fns[l-1]))
            }
        }
        X_pred <- cbind(1,activation(new_input %*%nn[[nlayers]]$wgts, fx=act.fns[nlayers]))
        X_preds <- as.matrix(cbind(X_pred, preds_spatial))
        preds <- rnorm(nrow(X_preds), X_preds %*% beta, 1)
        
        ## Save parameters if necessary
        if(it %in% kp.seq){
            kp <- kp+1
            beta.draws[kp,] <- beta
            if(nzones>2){
                cut.draws[kp,] <- cuts[est.cuts]
                predictions[kp, ] <-  cut(preds, breaks=c(-Inf, 0, cut.draws[kp,], Inf), labels=1:nzones)
            } else {
                predictions[kp, ] <-  cut(preds, breaks=c(-Inf, 0, Inf), labels=1:nzones)
            }
            
            for(i in 1:nlayers){
                weight.draws[[i]][,,kp] <- nn[[i]]$wgts
            }
        }
        
        ## Increment progress bar
        setTxtProgressBar(pb, it)
    }
    close(pb)

    if(feature_importance == TRUE){
        n_beta <- nrow(beta.draws)
        if(spatial_basis == 0){
            spatial <- read_csv("../Data/spatial_train.csv", col_types = cols())[,0]
        } else{
            spatial <- read_csv("../Data/spatial_train.csv", col_types = cols())[,1:spatial_basis]
        }
        
        ag <- read.csv("../Data/ag_train.csv")
        ag <- ag %>% dplyr::select(ObjectId, X, Y, ELEVATION, kryield19, SlopeDEg_m, 
                                   AspectDeg_, AspectClas, krVWC8_29, krVWC5_31, krVWC6_25a, 
                                   krVWC4_29, krslopedeg, Kr19NDVI_m, Kr18NDVI_m, TWI, VWCchange) %>%
            mutate(AvgVWC = (krVWC8_29 + krVWC5_31 + krVWC6_25a + krVWC4_29)/4,
                   Aspect.x = cos(AspectDeg_),
                   Aspect.y = sin(AspectDeg_))
        
        # Split the AvgVWC variable into three groups: Low, Medium, and High
        ag <- ag %>% mutate(Water = cut(AvgVWC, 
                                        breaks = quantile(AvgVWC, probs=seq(0, 1, length=nzones+1)),
                                        labels = 1:nzones,
                                        include.lowest=TRUE)) %>% 
            dplyr::select(-AvgVWC)
        model_matrix <- model.matrix(Water~.-1,
                                     data=ag %>% 
                                         dplyr::select(ELEVATION, kryield19, SlopeDEg_m, Aspect.x, Aspect.y,
                                                       Kr19NDVI_m, krslopedeg, TWI, Water, Kr18NDVI_m))
        input_scaled <-apply(model_matrix, 2, scales::rescale, to=c(-1, 1))
        preds <- matrix(0, nrow = n_beta, ncol = nrow(model_matrix))
        for(i in 1:n_beta){
            for(l in 1:nlayers){
                if(l == 1){
                    new_input <- input_scaled
                } else{
                    new_input <- cbind(1, activation(new_input %*% weight.draws[[l-1]][,,i], fx = act.fns[l-1]))
                }
            }
            X <- cbind(1,activation(new_input %*% weight.draws[[nlayers]][,,i], fx=act.fns[nlayers]))
            X <- cbind(X, spatial)
            z_s <- as.matrix(X) %*% as.matrix(beta.draws[i,], ncol = 1)
            preds[i, ] <-  cut(z_s, breaks=c(-Inf, 0, cut.draws[i,], Inf), labels=1:nzones)
        }
        
        pred_zone <- apply(preds, 2, getmode)
        model_rand <- adj.rand.index(pred_zone, ag$Water)
        explan_vars <- c("ELEVATION","kryield19", "SlopeDEg_m", "Aspect.x",
                         "Aspect.y", "Kr19NDVI_m", "krslopedeg", "TWI", "Kr18NDVI_m")
        
        permutations <- permutations
        feature_importance_matrix <- matrix(0, nrow = permutations, ncol = length(explan_vars))
        rand_index_vars <- numeric(length(explan_vars))
        for(permute in 1:permutations){
            preds2 <- matrix(0, nrow = n_beta, ncol = nrow(model_matrix))
            pb <- txtProgressBar(min = 0, max = 10, style = 3)
            for(k in 1:length(explan_vars)){
                var <- explan_vars[k]
                input_scaled_permuted <- input_scaled
                input_scaled_permuted[,var] <- input_scaled_permuted[,var][sample(1:nrow(ag), nrow(ag), replace = FALSE)]
                for(i in 1:n_beta){
                    for(l in 1:nlayers){
                        if(l == 1){
                            new_input <- input_scaled_permuted
                        } else{
                            new_input <- cbind(1, activation(new_input %*% weight.draws[[l-1]][,,i], fx = act.fns[l-1]))
                        }
                    }
                    X <- cbind(1,activation(new_input %*% weight.draws[[nlayers]][,,i], fx=act.fns[nlayers]))
                    X <- cbind(X, spatial)
                    z_s <- as.matrix(X) %*% as.matrix(beta.draws[i,], ncol = 1)
                    preds2[i, ] <-  cut(z_s, breaks=c(-Inf, 0, cut.draws[i,], Inf), labels=1:nzones)
                }
                pred_zone <- apply(preds2, 2, getmode)
                rand_index_vars[k] <- adj.rand.index(pred_zone, ag$Water)
                setTxtProgressBar(pb, k)
            }
            
            feature_importance_matrix[permute,] <- model_rand - rand_index_vars
        }
        
        feature_importance_all <- apply(feature_importance_matrix, 2, mean)
        graph <- data.frame(vars = factor(explan_vars), importance = feature_importance_all)
        p1 <- graph %>% 
            mutate(vars = fct_reorder(vars, importance)) %>% 
            ggplot(aes(x = vars, y = importance)) +
            geom_bar(stat = "identity", fill = "wheat") +
            scale_x_discrete(labels = rev(c("Elevation", "Yield 2019", "NDVI 2018", 
                                            "NDVI 2019", "Aspect Y",
                                            "TWI", "Slope Measure 1", "Slope Measure 2", "Aspect X"))) +
            coord_flip() +
            theme_minimal() +
            labs(y = "Feature Importance", x = "Variable")
    }
    else{
        feature_importance_matrix <- "None"
        p1 <- "No plot"
    }
    if(pdp == TRUE){
        if(spatial_basis == 0){
            spatial <- read_csv("../Data/spatial_train.csv", col_types = cols())[,0]
        } else{
            spatial <- read_csv("../Data/spatial_train.csv", col_types = cols())[,1:spatial_basis]
        }
        explan_vars <- c("ELEVATION", "kryield19", "SlopeDEg_m", "Aspect.x",
                         "Aspect.y", "Kr19NDVI_m", "krslopedeg", "TWI", "Kr18NDVI_m")
        nvals <- 100
        avg_zs <- array(0, c(n_beta, nvals, length(explan_vars)))
        pb <- txtProgressBar(min = 0, max = 10, style = 3)
        for(k in 1:length(explan_vars)){
            var <- explan_vars[k]
            ag <- suppressMessages(read_csv("../Data/ag_train.csv"))
            ag <- ag %>% dplyr::select(ObjectId, X, Y, ELEVATION, kryield19, SlopeDEg_m, 
                                       AspectDeg_, AspectClas, krVWC8_29, krVWC5_31, krVWC6_25a, 
                                       krVWC4_29, krslopedeg, Kr19NDVI_m, Kr18NDVI_m, TWI, VWCchange) %>%
                mutate(AvgVWC = (krVWC8_29 + krVWC5_31 + krVWC6_25a + krVWC4_29)/4,
                       Aspect.x = cos(AspectDeg_),
                       Aspect.y = sin(AspectDeg_))
            
            # Split the AvgVWC variable into three groups: Low, Medium, and High
            ag <- ag %>% mutate(Water = cut(AvgVWC, 
                                            breaks = quantile(AvgVWC, probs=seq(0, 1, length=nzones+1)),
                                            labels = 1:nzones,
                                            include.lowest=TRUE)) %>% 
                dplyr::select(-AvgVWC)
            
            input_vector <- seq(min(ag[,var]), max(ag[,var]), length.out = 100)
            vector_scaled <- scales::rescale(input_vector, to=c(-1,1))
            model_matrix <- model.matrix(Water~.-1,
                                         data=ag %>% 
                                             dplyr::select(ELEVATION, kryield19, SlopeDEg_m, Aspect.x, Aspect.y,
                                                           Kr19NDVI_m, krslopedeg, TWI, Kr18NDVI_m, Water))
            input_scaled <-apply(model_matrix, 2, scales::rescale, to=c(-1, 1))
            for(i in 1:n_beta){
                for(j in 1:length(vector_scaled)){
                    input_scaled[,var] <- vector_scaled[j]
                    for(l in 1:nlayers){
                        if(l == 1){
                            new_input <- input_scaled
                        } else{
                            new_input <- cbind(1, activation(new_input %*% weight.draws[[l-1]][,,i], fx = act.fns[l-1]))
                        }
                    }
                    X <- cbind(1,activation(new_input %*% weight.draws[[nlayers]][,,i], fx=act.fns[nlayers]))
                    X <- cbind(X, spatial)
                    z_s <- as.matrix(X) %*% as.matrix(beta.draws[i,], ncol = 1)
                    avg_zs[i, j, k] <- mean(z_s)
                }
            }
            setTxtProgressBar(pb, k)
        }
        plots <- list()
        for(var in explan_vars){
            plots[[length(plots) + 1]] <- pdp_plot(var, avg_zs)
        }
    }
    else {
        avg_zs <- "Not calculated"
        plots <- "None"
    }
    
    prob_1_s <- apply(predictions, 2, function(x){mean(x == 1)})
    prob_2_s <- apply(predictions, 2, function(x){mean(x == 2)})
    prob_3_s <- apply(predictions, 2, function(x){mean(x == 3)})
    ag <- readxl::read_excel("../Data/MultiVariable_withZones.xlsx")
    nzones <- 3
    ag <- ag %>% mutate(AvgVWC = (krVWC8_29 + krVWC5_31 + krVWC6_25a + krVWC4_29)/4,
                        Water = cut(AvgVWC, breaks = quantile(AvgVWC, probs=seq(0, 1, length=nzones+1)),
                                    labels = 1:nzones,
                                    include.lowest=TRUE)) %>% 
        dplyr::select(-AvgVWC)
    
    plot_prob1 <- ggplot(ag, aes(x = X, y = Y, fill = prob_1_s)) +   geom_tile() + 
        scale_fill_distiller(palette = "Spectral", na.value = NA) +
        labs(fill = "VWC", x = "Longitude", y = "Latitude", title = "Probability: Zone 1") + theme_minimal()
    plot_prob2 <- ggplot(ag, aes(x = X, y = Y, fill = prob_2_s)) +   geom_tile() + 
        scale_fill_distiller(palette = "Spectral", na.value = NA) +
        labs(fill = "VWC", x = "Longitude", y = "Latitude", title = "Probability: Zone 2") + theme_minimal()
    plot_prob3 <- ggplot(ag, aes(x = X, y = Y, fill = prob_3_s)) +   geom_tile() + 
        scale_fill_distiller(palette = "Spectral", na.value = NA) +
        labs(fill = "VWC", x = "Longitude", y = "Latitude", title = "Probability: Zone 3") + theme_minimal()
    preds <- apply(predictions, 2, getmode)
    plot_pred <- ggplot(ag, aes(x = X, y = Y, fill = factor(preds))) +   geom_tile() + 
        scale_fill_brewer(palette = "Blues", na.value = NA) +
        labs(fill = "VWC", x = "Longitude", y = "Latitude", title = "Predictions") + theme_minimal()
    
    prediction_plots <- list(plot_prob1, plot_prob2, plot_prob3, plot_pred) 
    
    y <- 1*prob_1_s + 2*prob_2_s + 3*prob_3_s
    s <- cbind(ag$X, ag$Y)
    boundary_pts <- c(min(ag$X), max(ag$X), min(ag$Y), max(ag$Y))
    clust <- 10
    clust.agg <- agg.gradclust(y,s,clust,boundary_pts ,num.gpts=45)
    
    map_data <- data.frame(Long = ag$X, Lat = ag$Y, clust_agg = clust.agg)
    preds <- numeric(length = nrow(data))
    for(i in 1:length(prob_1_s)){
        preds[i] <- which.max(c(prob_1_s[i], prob_2_s[i], prob_3_s[i]))
    }
    
    clusters <- map_data %>% mutate(cluster = clust.agg) %>% group_by(cluster) %>%
        summarize(avg_prob = median(y), 
                  dist_1 = abs(avg_prob - 1.40),
                  dist_2 = abs(avg_prob - 1.94),
                  dist_3 = abs(avg_prob - 2.77)
        )
    zones <- numeric(length = nrow(clusters))
    for(i in 1:nrow(clusters)){
        zones[i] <- which.min(c(clusters$dist_1[i], clusters$dist_2[i], clusters$dist_3[i]))
    }
    clusters <- clusters %>% mutate(zone = zones)
    overall <- ag %>% mutate(cluster = clust.agg) %>% inner_join(clusters, by = "cluster")
    s_clust <- ggplot(overall, aes(x = X, y = Y, fill = factor(zone))) +   geom_tile() + 
        scale_fill_brewer(palette = "Blues", na.value = NA) +
        labs(fill = "Zone", x = "Longitude", y = "Latitude", title = "Spatial Clustering") + theme_minimal()
    
    return(list(betas = beta.draws, cut_points = cut.draws, weights = weight.draws, 
                preds = predictions, 
                feature_importance = list(fi = feature_importance_matrix, plot = p1),
                pdp = list(avg_z = avg_zs,
                           plots = plots),
                prediction_plots = prediction_plots,
                final_zones = s_clust))

}
