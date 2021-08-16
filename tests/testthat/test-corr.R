test_that("multiplication works", {
  expect_equal(corr.test(data.frame(a=1:3,b=1:3),
                         data.frame(c=2:4,d=2:4))$r,
               matrix(1,nrow=2,ncol=2,
                      dimnames = list(c('a','b'),
                                     c('c','d'))))
})
