#require(testthat)
context("fSampleFunc")

test_that("two vectors are added correctly",{
		a = 1:10
		b = 1
		c = a+b
		expect_that( fSampleFunc(a,b), equals(c) )
		expect_that( fSampleFunc(a), throws_error() )
	})

