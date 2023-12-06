## ----setup, include=FALSE-----------------------------------------------------
## Global options
knitr::opts_chunk$set(cache = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("rmdformats")

## -----------------------------------------------------------------------------
grade = data.frame(name = c("小明","小红","小方"), chinese = c(90,89,93), math = c(88,95,89), english = c(94,95,91))

## -----------------------------------------------------------------------------
knitr::kable(grade)

## -----------------------------------------------------------------------------
grade$total = grade$chinese + grade$math + grade$english
knitr::kable(grade)

## ----message = FALSE----------------------------------------------------------
library(tidyverse)

## -----------------------------------------------------------------------------
mpg %>% 
  ggplot(aes(x=class, y=hwy, fill=class)) + 
  geom_boxplot()

