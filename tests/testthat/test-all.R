library(GPSCDF)
library(testthat)
context("tests GPSCDF")


x1<- c(.2,.3,.5)
x2<- c(.1,.6,.3)
x3<- c(.6,.3,.1)
x4<- c(.7,.1,.2)
x5<- c(.33,.33,.34)
x6<- c(.35,.25,.4)

set.seed(18201)
trt<- sample(1:3, 6, replace=TRUE)
id<- sort(sample(1:6,6))

#Data inputs
X<- cbind(id, trt, as.data.frame(rbind(x1,x2,x3,x4,x5,x6)))

scores<- rbind(x1,x2,x3,x4,x5,x6)


test_that("tests sum of GPS", {
  scores[6,3]<- 0.5
  expect_message(GPSCDF(pscores=scores, data=X), "Ensure pscores sum to 1")


  scores[6,3]<- 0.3
  expect_message(GPSCDF(pscores=scores, data=X), "Ensure pscores sum to 1")
})


test_that("tests GPSCDF ppar", {
  fit<- GPSCDF(pscores=scores, data=X)

  expect_equal(fit$ppar[1], 0.46667569, tolerance = 0.00001)
  expect_that(length(fit$ppar), is_equivalent_to(6))
})


test_that("tests GPSCDF stratification", {
  fit1<- GPSCDF(pscores=scores, stratify=TRUE)
  fit2<- GPSCDF(pscores=scores, stratify=TRUE, nstrat=3)

  expect_identical(fit1$nstrat, 5)
  expect_identical(fit2$nstrat, 3)
})


test_that("tests GPSCDF optimal matching", {

  fit3<- GPSCDF(pscores=scores, data=X, trt=X$trt, optimal=TRUE, ordinal=TRUE)
  fit4<- GPSCDF(pscores=scores, data=X, trt=X$trt, optimal=TRUE, multinomial=TRUE)

  expect_error(GPSCDF(pscores=scores, optimal=TRUE), "Specify a dataframe to attach matches")
  expect_error(GPSCDF(pscores=scores, data=X, optimal=TRUE), "Specify a treatment variable to proceed with matching")
  expect_error(GPSCDF(pscores=scores, data=X, trt=X$trt, optimal=TRUE), "Specify Ordinal or Multinomial treatments")

  expect_equal(fit3$optdistance, 0.1147744, tolerance = 0.00001)
  expect_equal(fit4$optdistance, 0.1147744, tolerance = 0.00001)
})


test_that("tests GPSCDF greedy matching", {

  fit5<- GPSCDF(pscores=scores, data=X, trt=X$trt, greedy=TRUE, ordinal=TRUE)
  fit6<- GPSCDF(pscores=scores, data=X, trt=X$trt, greedy=TRUE, multinomial=TRUE)

  expect_error(GPSCDF(pscores=scores, greedy=TRUE), "Specify a dataframe to attach matches")
  expect_error(GPSCDF(pscores=scores, data=X, greedy=TRUE), "Specify a treatment variable to proceed with matching")
  expect_error(GPSCDF(pscores=scores, data=X, trt=X$trt, greedy=TRUE), "Specify Ordinal or Multinomial treatments")

})

