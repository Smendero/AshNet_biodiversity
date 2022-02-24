#### Gain Loss Analysis Functions ####
#### Written by Emily Smenderovac ####

GLglm <- function(lr, conditions, verbose = FALSE, mc.samples = 128){
       
        ## Run a GLM on a Gain-Loss Matrix, represented as %of replicates where a zOTU is gained or lost 
        
        lr2glm <- function(lr, conditions){
                ## Adapted from aldex.glm - Montecarlo permutations removed
                if (nrow(lr) != nrow(conditions)) {
                        stop("Input data and 'model.matrix' should have same number of rows.")
                }
                model. <- conditions
                glms <- apply(lr, 2, function(x) {
                        glm(x ~ model.)
                })
                extract <- function(model) {
                        x <- coef(summary(model))
                        coefs <- lapply(1:nrow(x), function(i) {
                                y <- x[i, , drop = FALSE]
                                colnames(y) <- paste(rownames(y), colnames(y))
                                y
                        })
                        do.call("cbind", coefs)
                }
                extracts <- lapply(glms, extract)
                df <- do.call("rbind", extracts)
                rownames(df) <- colnames(lr)
                df <- as.data.frame(df)
                pvals <- colnames(df)[grepl("Pr\\(>", colnames(df))]
                df.bh <- df[, pvals]
                colnames(df.bh) <- paste0(colnames(df.bh), ".BH")
                for (j in 1:ncol(df.bh)) {
                        df.bh[, j] <- p.adjust(df.bh[, j], method="BH")
                }
                data.frame(cbind(df, df.bh))
        }
        rdirichlet<-function(n,a)
                ## pick n random deviates from the Dirichlet function with shape parameters a
        {
                if(length(n) > 1 || length(n) < 1 || n < 1) stop("n must be a single positive integer value")
                if(length(a) < 2) stop("a must be a vector of numeric value")
                n <- floor(n)
                l<-length(a);
                x<-matrix(rgamma(l*n,a),ncol=l,byrow=TRUE);
                # sm<-x%*%rep(1,l);
                return((x-100)/100);# Back transform to Gain/Loss
        }
        ## Monte-Carlo permutations of Dirichlet distributions of Gains/Losses
        #  remove all rows with reads less than the minimum set by minsum
        minsum <- 0
        
        # remove any row in which the sum of the row is 0
        z <- as.numeric(apply(lr, 1, sum))
        lr <- as.data.frame( lr[(which(z > minsum)),]  )
        
        rn <- rownames( lr )
        
        mc <- lapply( lr , 
                      function(col) {
                              q <- t( rdirichlet( mc.samples, col ) ) ;
                              rownames(q) <- rn ; q } )
        r=0
        
        for(i in 1:mc.samples){
                mci_lr <- t(sapply(mc, function(x) x[, i]))
                r <- r + lr2glm(mci_lr, conditions)
        }
        r/mc.samples
}



## Change Log ##

## 2021-04-14
## Created initial draft of GLglm function from adapted aldex.glm code