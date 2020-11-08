library( rfishbase )
library( ggplot2 )


fam <- "Sparidae"
fam2 <- "Cottidae"

spec <- species_list( Family=fam2 )

#ecosystem(spec)
#ecol.spec <- ecology(spec)
morp.spec <- morphometrics(spec)
morp.spec$SnoutTipX
morp.spec$SnoutTipY
matu.spec <- maturity(spec)
matu.spec$LengthMatMin
matu.spec$Lm
with(matu.spec, cbind(LengthMatMin, LengthMatMin2, Lm))


morp.list <- morp.spec


# copied from fishbase.org webpages
test <- read.csv( "test.csv", stringsAsFactors = FALSE, header=TRUE )
head(test)
str(test)
# only keep rows with family names
test <- test[ test$family != "", ]

# compare some of the traits
ggplot( data=test, aes(x=family,y=generation) ) + geom_violin() 
ggplot( data=test, aes(x=family,y=Linf) ) + geom_violin()
ggplot( data=test, aes(x=family,y=life) ) + geom_violin()

