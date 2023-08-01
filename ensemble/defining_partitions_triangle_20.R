###########################
##Note: the numbering of the ADLP ensemble is not the same as the numbering of ADLP ensemble in the paper due to better presentation in paper
##########################
####please refer to the following index section for the actual partition of subsets

################################################################################
claims_df_train_val_test_20 <- function(df) {
    
    df <- df.to.numeric(df)
    # Subsets claims dataframe into 6 training, validation and test sets
    # Construction of Training Set:
    train_1 <- df[df$calendar<=17, ]
    train_2 <- df[df$origin<=1 & 17<=df$dev, ]
    train_3 <- df[df$dev<=1 & 17<=df$origin, ]
    train_4 <- df[df$origin==16 & df$dev==2, ]
    train<-rbind(train_1,train_2,train_3,train_4)
    train <- df.to.factor(train)
    train[order(as.numeric(as.character(train$origin))),]
    
    # Construction of Validation set:
    valid_1 <- df[18<=df$calendar & df$calendar<=21 & 2<=df$origin & 3<=df$dev, ]
    valid_2 <- df[df$dev==2 & 17<=df$origin & df$origin<=19, ]
    valid<-rbind(valid_1,valid_2)
    valid <- df.to.factor(valid)
    valid[order(as.numeric(as.character(valid$origin))),]
    
    # Construction of Test Data
    test <- df[22<=df$calendar, ]
    test <- df.to.factor(test)
    test<-test[order(as.numeric(as.character(test$origin))),]
    
    return(list(train=train, valid=valid, test=test))
}

claims_df_in_out_20 <- function(df) {
    return(list(
        train = df[as.numeric(df$calendar) <= 21, ], 
        test = df[as.numeric(df$calendar) > 21,]
    ))
}
################################################################################
### Size 20 triangle partitions
################################################################################

# Two Subset Partitions

split_points_20 <- c(4, 6, 8, 10, 12, 14, 16, 18)

par0 <- function(df) {
    return (list(
        subset1 = df
    ))
}

par1_2_20 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_20[1]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_20[1] + 1) & (as.numeric(as.character(df$origin)) <= 20), ]
    ))
}

par2_2_20 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_20[2]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_20[2] + 1) & (as.numeric(as.character(df$origin)) <= 20), ]
    ))
}

par3_2_20 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_20[3]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_20[3] + 1) & (as.numeric(as.character(df$origin)) <= 20), ]
    ))
}

par4_2_20 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_20[4]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_20[4] + 1) & (as.numeric(as.character(df$origin)) <= 20), ]
    ))
}

par5_2_20 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_20[5]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_20[5] + 1) & (as.numeric(as.character(df$origin)) <= 20), ]
    ))
}

par6_2_20 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_20[6]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_20[6] + 1) & (as.numeric(as.character(df$origin)) <= 20), ]
    ))
}

par7_2_20 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_20[7]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_20[7] + 1) & (as.numeric(as.character(df$origin)) <= 20), ]
    ))
}

par8_2_20 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_20[8]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_20[8] + 1) & (as.numeric(as.character(df$origin)) <= 20), ]
    ))
}
