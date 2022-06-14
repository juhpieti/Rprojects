library(mlegp)
library(ggplot2)
load("toJuho.Rdata")

##### Fitting GP, predicting ######

df <- as.data.frame(toJuho)
View(df)
colnames(df)[6] <- "likelihood"
plot(df$likelihood)
dim(df) # 750 x 6 

shuffled_df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
plot(shuffled_df$likelihood, ylab = "likelihood")

training_set <- shuffled_df[1:500,]
test_set <- shuffled_df[501:nrow(shuffled_df), ]

start_clock <- proc.time()

GPmodel = mlegp(X = training_set[, -ncol(df), drop = FALSE], Z = training_set[, ncol(df), drop = FALSE], 
                nugget = 0, nugget.known = 1, verbose = 0) #sets nugget to 0 --> no variance on training points? test verbose 1, 2

stop_clock <- proc.time()
time <- stop_clock - start_clock
time[3] #384.934 seconds to fit the model

start_clock <- proc.time()
GP_pred = predict(GPmodel, shuffled_df[,-ncol(df)], se.fit = TRUE)
stop_clock <- proc.time()
time <- stop_clock - start_clock
time[3] #906.408 seconds to predict the test set

shuffled_df <- cbind(shuffled_df, GP_pred$fit, GP_pred$se.fit) #see how the predictions match with the training data, differs in test data

colnames(shuffled_df)[7:ncol(shuffled_df)] <- c("pred", "se")

test_set <- shuffled_df[501:nrow(shuffled_df), ]

sqrt(mean((test_set$likelihood - test_set$pred)^2)) #156.3078 (RMSE root-mean-square-error)

ggplot(test_set, mapping = aes(x = as.integer(row.names(test_set)), y = likelihood)) +
  geom_point(color = "black") +
  geom_point(aes(y = pred), color = "red") + 
  geom_segment(mapping = aes(xend = as.integer(row.names(test_set)), yend = pred), alpha = 0.1) +
  labs(x = "index") #+
  #scale_colour_manual(values = c("black" = "black", "red" = "red"), labels = c("likelihood", "prediction"))


### marginal plots of predicted likelihoods by different predictors
par(mfrow = c(2,3))
for (i in 1:5) {
  print(i)
  plot(shuffled_df[,i], shuffled_df$pred, ylab = "predicted likelihood", xlab = paste(colnames(shuffled_df)[i]))
} 

