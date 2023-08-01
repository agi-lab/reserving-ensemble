###########################
##Note: the numbering of the ADLP ensemble is not the same as the numbering of ADLP ensemble in the paper due to better presentation in paper
##########################
####please refer to the following index section for the actual partition of subsets

#full_data <- read.csv('simulation/triangle_40-data/sim1-full-data.csv')
#parition_in_idx <- full_data$calendar <= 41
#parition_out_idx <- full_data$calendar > 41

################################################################################
claims_df_train_val_test_40 <- function(df) {
    
    df <- df.to.numeric(df)
    # Subsets claims dataframe into 6 training, validation and test sets
    # Construction of Training Set:
    train_1 <- df[df$calendar<=34, ]
    train_2 <- df[df$origin==1 & 34<=df$dev & df$dev<=40, ]
    train_3 <- df[35<=df$calendar & df$calendar<=36 & df$dev <=2, ]
    train_4 <- df[36<=df$origin & df$dev==1, ]
    train<-rbind(train_1,train_2,train_3,train_4)
    train <- df.to.factor(train)
    train[order(as.numeric(as.character(train$origin))),]
    
    # Construction of Validation set:
    valid_1 <- df[35<=df$calendar & df$calendar<=41 & 2<=df$origin & 3<=df$dev, ]
    valid_2 <- df[35<=df$origin & df$origin<=39 & df$dev == 2, ]
    valid<-rbind(valid_1,valid_2)
    valid <- df.to.factor(valid)
    valid[order(as.numeric(as.character(valid$origin))),]
    
    # Construction of Test Data
    test<-df[42<=df$calendar, ]
    test <- df.to.factor(test)
    test<-test[order(as.numeric(as.character(test$origin))),]
    
    return(list(train=train, valid=valid, test=test))
}

claims_df_in_out_40 <- function(df) {
    return(list(
        train = df[as.numeric(df$calendar) <= 41, ], 
        test = df[as.numeric(df$calendar) > 41,]
    ))
}
################################################################################
### Size 40 triangle partitions
################################################################################

# Two Subset Partitions

split_points_40 <- c(3, 4, 5, 7, 9, 11, 13, 14, 15, 16, 17, 18, 19, 23, 26, 28, 31, 33)

par0 <- function(df) {
    return (list(
        subset1 = df
    ))
}

par1_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[1]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[1] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par2_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[2]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[2] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par3_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[3]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[3] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par4_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[4]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[4] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par5_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[5]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[5] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par6_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[6]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[6] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par7_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[7]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[7] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par8_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[8]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[8] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par9_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[9]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[9] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par10_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[10]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[10] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par11_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[11]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[11] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par12_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[12]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[12] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par13_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[13]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[13] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par14_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[14]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[14] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par15_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[15]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[15] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par16_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[16]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[16] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par17_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[17]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[17] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par18_2_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40[18]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40[18] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}


# Three Subset Partitions

split_points_40_3 <- list(c(5, 15), c(15, 29), c(17, 31), c(23, 33))

par1_3_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40_3[[1]][1]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40_3[[1]][1] + 1) & (as.numeric(as.character(df$origin)) <= split_points_40_3[[1]][2]), ],
        subset3 = df[(as.numeric(as.character(df$origin)) >= split_points_40_3[[1]][2] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par2_3_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40_3[[2]][1]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40_3[[2]][1] + 1) & (as.numeric(as.character(df$origin)) <= split_points_40_3[[2]][2]), ],
        subset3 = df[(as.numeric(as.character(df$origin)) >= split_points_40_3[[2]][2] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par3_3_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40_3[[3]][1]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40_3[[3]][1] + 1) & (as.numeric(as.character(df$origin)) <= split_points_40_3[[3]][2]), ],
        subset3 = df[(as.numeric(as.character(df$origin)) >= split_points_40_3[[3]][2] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}

par4_3_40 <- function(df) {
    return (list(
        subset1 = df[(as.numeric(as.character(df$origin)) >= 2) & (as.numeric(as.character(df$origin)) <= split_points_40_3[[4]][1]), ],
        subset2 = df[(as.numeric(as.character(df$origin)) >= split_points_40_3[[4]][1] + 1) & (as.numeric(as.character(df$origin)) <= split_points_40_3[[4]][2]), ],
        subset3 = df[(as.numeric(as.character(df$origin)) >= split_points_40_3[[4]][2] + 1) & (as.numeric(as.character(df$origin)) <= 40), ]
    ))
}
