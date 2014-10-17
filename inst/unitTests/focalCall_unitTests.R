library("focalCall")

# Load data
data(BierkensCNA)
print(CGHset)
stopifnot(inherits(CGHset, "cghCall"))

print(CNVset)
stopifnot(inherits(CNVset, "data.frame"))

ExampleRun<-focalCall(CGHset, CNVset, focalSize=3, minFreq=2)
print(ExampleRun)
stopifnot(inherits(ExampleRun, "cghCall"))

