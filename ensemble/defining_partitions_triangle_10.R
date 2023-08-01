###########################
##Note: the numbering of the ADLP ensemble is not the same as the numbering of ADLP ensemble in the paper due to better presentation in paper
##########################
####please refer to the following index section for the actual partition of subsets

################################################################################
claims_df_train_val_test_10 <- function(df) {
    
    df <- df.to.numeric(df)
    # Subsets claims dataframe into 6 training, validation and test sets
    # Construction of Training Set:
    train_1 <- df[df$calendar<=9, ]
    train_2 <- df[df$origin==1 & 9<=df$dev, ]
    train_3 <- df[9<=df$origin & df$dev==1, ]
    train<-rbind(train_1,train_2,train_3)
    train <- df.to.factor(train)
    train[order(as.numeric(as.character(train$origin))),]
    
    # Construction of Validation set:
    valid <- df[10<=df$calendar & df$calendar<=11 & 2<=df$origin & 2<=df$dev, ]
    valid <- df.to.factor(valid)
    valid[order(as.numeric(as.character(valid$origin))),]
    
    # Construction of Test Data
    test <- df[12<=df$calendar, ]
    test <- df.to.factor(test)
    test<-test[order(as.numeric(as.character(test$origin))),]
    
    return(list(train=train, valid=valid, test=test))
}

claims_df_in_out_10 <- function(df) {
    return(list(
        train = df[as.numeric(df$calendar) <= 11, ], 
        test = df[as.numeric(df$calendar) > 11,]
    ))
}
################################################################################
### Size 10 triangle partitions
################################################################################

# Two Subset Partitions

split_points_10 <- c(3, 4, 5, 6, 7, 8)

par0 <- function(df) {
    return (list(
        subset1 = df
    ))
}

par1_2_10 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_10[1]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_10[1] + 1) & (as.numeric(as.character(df$origin)) <= 10), ]
    ))
}

par2_2_10 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_10[2]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_10[2] + 1) & (as.numeric(as.character(df$origin)) <= 10), ]
    ))
}

par3_2_10 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_10[3]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_10[3] + 1) & (as.numeric(as.character(df$origin)) <= 10), ]
    ))
}

par4_2_10 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_10[4]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_10[4] + 1) & (as.numeric(as.character(df$origin)) <= 10), ]
    ))
}

par5_2_10 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_10[5]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_10[5] + 1) & (as.numeric(as.character(df$origin)) <= 10), ]
    ))
}

par6_2_10 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_10[6]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_10[6] + 1) & (as.numeric(as.character(df$origin)) <= 10), ]
    ))
}