#' @title Incorporates historical control into the meta analysis
#' @description This function takes the data into N-n form and run meta analysis using Begg & Pilote Method and DarSimonian & Laird method
#' @param data A data frame having 5 columns 1. Study.Name, 2. Trt1_N, 3. Trt1_n, 4. Trt2_N and 5. Trt2_n. The column order is not important but the column name must be same as specified above
#' @param method which method user wants to use specify one among 1. beggPilote or 2. DerLiard. Default is "beggPilote
#' @param alpha alpha value for confidence interval. Default is 0.95 results 95% confidence Interval
#' @author Shubhram Pandey
#' @export
metaHistoric = function(data,method = "beggPilote", alpha = 0.95, roundDigits = 4){
  if (!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
      #To ignore the warnings during usage
      options(warn = -1)
      options("getSymbols.warning4.0" = FALSE)

      data = as.data.frame(data)
      colnames = c("Study.Name","Trt1_N","Trt1_n","Trt2_N","Trt2_n")
      defaultMethod = c("beggPilote","DerLiard")
      if (!identical(colnames(data),colnames)) {
        print("Column names are not in standard format")
      } else {
        if (!(method %in% defaultMethod)) {
          print("Method is not correctlty specified")
        } else if (method == "beggPilote") {
          totStudy = nrow(data)
          data = data %>%
            mutate(propTrt1 = Trt1_n/Trt1_N,
                   propTrt2 = Trt2_n/Trt2_N,
                   seTrt1 = sqrt((propTrt1*(1 - propTrt1))/Trt1_n),
                   seTrt2 = sqrt((propTrt2*(1 - propTrt2))/Trt2_n),
                   weightInterTrt1 = 1/(seTrt1^2),
                   weightInterTrt2 = 1/(seTrt2^2),
                   weightsTrt1 = weightInterTrt1/sum(weightInterTrt1,na.rm = T),
                   weightsTrt2 = weightInterTrt2/sum(weightInterTrt2,na.rm = T)
            )
          estimateTrt1 = sum(data$propTrt1 * data$weightsTrt1, na.rm = T)
          estimateTrt2 = sum(data$propTrt2 * data$weightsTrt2, na.rm = T)
          diffEstimate = estimateTrt2 - estimateTrt1
          seEstimate1 = sqrt(1/sum(data$weightInterTrt1,na.rm = T))
          seEstimate2 = sqrt(1/sum(data$weightInterTrt2,na.rm = T))
          seCombined = sqrt(seEstimate1^2 + seEstimate2^2 )
          lowerBound = diffEstimate - (qnorm(alpha + (1 - alpha)/2) * seCombined)
          upperBound = diffEstimate + (qnorm(alpha + (1 - alpha)/2) * seCombined)
          results = data.frame("Treatment1" = c(round(estimateTrt1,roundDigits),round(seEstimate1,roundDigits),"","",""),
                               "Treatment2" = c(round(estimateTrt2,roundDigits),
                                                round(seEstimate2,roundDigits),
                                                round(diffEstimate,roundDigits),
                                                round(seCombined,roundDigits),
                                                paste0(alpha*100,"% CI: (",
                                                       round(lowerBound,roundDigits),
                                                       ",",
                                                       round(upperBound,roundDigits),
                                                       ")"
                                                )
                               )

          )
          rownames(results) = c("Pooled Estimate",
                                "Pooled SE",
                                "Diffrence in estimate",
                                "Combined SE",
                                "Confidence interval"
          )
        } else {
          data = data %>%
            mutate(propTrt1 = Trt1_n/Trt1_N,
                   propTrt2 = Trt2_n/Trt2_N,
                   seTrt1 = sqrt((propTrt1*(1 - propTrt1))/Trt1_n),
                   seTrt2 = sqrt((propTrt2*(1 - propTrt2))/Trt2_n),
                   maxSE = pmax(seTrt1,seTrt2,na.rm = T),
                   weightInter = 1/(maxSE^2),
                   weights = weightInter/sum(weightInter,na.rm = T),
                   wmd = propTrt2 - propTrt1
            )
          pooledEstimate = sum(data$weights * data$wmd, na.rm = T)
          pooledSE = sqrt(1/sum(data$weightInter,na.rm = T))
          lowerBound = pooledEstimate - (qnorm(alpha + (1 - alpha)/2) * pooledSE)
          upperBound = pooledEstimate + (qnorm(alpha + (1 - alpha)/2) * pooledSE)
          results = data.frame("Values" = c(round(pooledEstimate,roundDigits),
                                            round(pooledSE,roundDigits),
                                            paste0(alpha*100,"% CI: (",
                                                   round(lowerBound,roundDigits),
                                                   ",",
                                                   round(upperBound,roundDigits),
                                                   ")"
                                            )
          )

          )
          rownames(results) = c("Pooled Estimate",
                                "Pooled SE",
                                "Confidence interval"
          )
        }

      }
  return(list(data = data,results = results))
}
