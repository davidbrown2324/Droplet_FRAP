for(input in inputs) {
    tryCatch(print(paste("log of", input, "=", log(input))),
                 warning = function(w) {print(paste("negative argument", input)); 
                   log(-input)},
                 error = function(e) {print(paste("non-numeric argument", input));
                   NaN})
}


for(input in inputs) {
  tryCatch(print(paste("log of", input, "=", log(input))),
           warning = function(w) {print(paste("negative argument", input)); 
             log(-input)},
           error = function(e) {print(paste("non-numeric argument", input));
             NaN})
}

biexp_fit<-tryCatch(BiExpFit(time, frap),
                    warning = function(w) {print("watchout"); NaN},
                    error = function(e) {print("can't do it"); NaN})

f=1
testFit<-function(f){
  frap=Irel_frame[26:325,f]
  plot(time, frap)
  biexp_fit<-tryCatch(BiExpFit(time, frap),
                      warning = function(w) {print(paste0("warning for curve ", f));},
                      error = function(e) {print(paste0("can't process curve ", f)); NaN})
  biexp_fit<-try(fitted(biexp_fit))
  try(lines(time, biexp_fit, col='red'))
}

for (f in 1:46){
  testFit(f)
}
