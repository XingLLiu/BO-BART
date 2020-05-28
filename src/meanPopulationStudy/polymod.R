library("socialmixr")
data(polymod)
uk_polymod <- polymod$participants[polymod$participants$country == "United Kingdom"]
dim(uk_polymod)
uk_polymod
contact_matrix(polymod, age.limits = c(0, 1,5,15))
