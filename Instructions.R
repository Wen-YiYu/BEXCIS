#download the latest version of R
  install.packages("installr")
  library(installr)
  updateR()
  getwd()
  Sys.getenv("BINPREF")

#download the latest version of Rtools


# create a .Renviron file in the folder path of "document", writing the folder path of the "bin" in Rtools in it

  writeLines('PATH="${RTOOLS40_HOME}//usr//bin;${PATH}"', con = "~/.Renviron")

# succeed if the following command returns the path of "make.exe"
  Sys.which("make")


  install.packages("jsonlite", type = "source")

# run the following commands
  dotR <- file.path(Sys.getenv("HOME"), ".R")
  if (!file.exists(dotR)) dir.create(dotR)
  M <- file.path(dotR, "Makevars.win")
  if (!file.exists(M)) file.create(M)
  cat("\n CXX14FLAGS += -mtune=native -O3 -mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2",
    file = M, sep = "\n", append = FALSE)

# delete the former version of rstan

  remove.packages("rstan")
  if (file.exists(".RData")) file.remove(".RData")


# install the latest version of rstan

#Sys.setenv(MAKEFLAGS = "-j4") # four cores used
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)

  install.packages("rstan")

  library(rstan)

# the following warning is normal and can be ignored
# Warning message:
# In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#  'C:/rtools40/usr/mingw_/bin/g++' not found


####################################################################
# install the package
# setwd to the root of the folder path of "BEXCIS_1.0.zip" and run the following command
  install.packages("BEXCIS_1.0.zip",repos = NULL)

# or setwd to the root of the folder path of "BEXCIS_1.0.tar.gz" and run the following command
  install.packages("BEXCIS_1.0.tar.gz",type='source')

# load the package
  library(BEXCIS)

# check the help files
  ?BEXCIS

# setwd to the path of the example data to run the examples
  setwd(".../example data")

###########################################################################
##################    examples for Bayes_XCI()    #########################
###########################################################################
# The following warning is normal and can be ignored:

 Warning message:
 In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
  'C:/rtools40/usr/mingw_/bin/g++' not found

###########################################################################
## If the following error is returned:

 [1] "Error in sampler$call_sampler(args_list[[i]]) : Initialization failed."

 [1] "error occurred during calling the sampler; sampling not done".

Then, the value limits set for the coefficients of the covariates may be inappropriate.

  
  
# If the following error is returned: 

  Error in h(simpleError(msg, call)) : 

  error in evaluating the argument 'object' in selecting a method for function 'extract': object 'Fit_H1' not found.

Then, the prior distributions set for the coefficients of the covariates may be inappropriate.

###########################################################################
##example 1:
#the ped file without header
#set "ped='ped_qualitative.txt'" and "ped.header=F"

#the trait is qualitative
#set "trait_type='qualitative'"

#no covariate
#set "covariate=NULL"

#using uniform distribution as the prior distribution for gamma
#set "gamma_prior='uniform'"

#set seed to get the fixed result
  set.seed(123)
  example1<-Bayes_XCI(ped='ped_qualitative.txt', trait_missing=9, allele_missing=0,
					ped.header=F, trait_type='qualitative', covariate=NULL, 
					covariate_prior=NULL, covariate_missing=NA, covariate.header=F, 
					covariate_prior_limit=NULL, gamma_prior='uniform', chains_num=2, 
					total_sample=12000, warmup_sample=2000, adapt_delta_value=0.85)

#Although the number of Markov chains, the number of iterations for each
#chain, the number of warmup iterations and the target acceptance rate
#are respectively set to 2, 12000, 2000 and 0.85 in this example, we recommend
#respectively using 8, 20000, 10000 and 0.99.

#Result
  write.table(example1,"example1.txt",row.names=F,quote=F)
  print(example1)

#      Point_Estimate HPDI_Lower HPDI_Upper
#snp_1       1.715434  0.1687426   1.999994



##example 2:
#the ped file with header
#set "ped='ped_qualitative_header.txt'" and "ped.header=T"

#the trait is qualitative
#set "trait_type='qualitative'"

#the covariate file without header
#set "covariate='covariate.txt'" and "covariate.header=F"

#using normal distribution as the prior distribution for gamma
#set "gamma_prior='normal'"

#define the prior distributions of the coefficients of the covariates 
  covariate_prior<-rbind(c('uniform','-2,2'),c('normal','0,10'))

#define the value limits of the coefficients of the covariates
  covariate_prior_limit <- rbind(c(NA,'2'),c(NA,NA))

#set seed to get the fixed result  set.seed(123)
  example2<-Bayes_XCI(ped='ped_qualitative_header.txt', trait_missing=9, allele_missing=0,
			ped.header=T, trait_type='qualitative', covariate='covariate.txt',
			covariate_prior=covariate_prior, covariate_missing=NA, covariate.header=F, 
			covariate_prior_limit=covariate_prior_limit, gamma_prior='normal',
			chains_num=2, total_sample=12000, warmup_sample=2000, adapt_delta_value=0.85)

#Although the number of Markov chains, the number of iterations for each
#chain, the number of warmup iterations and the target acceptance rate
#are respectively set to 2, 12000, 2000 and 0.85 in this example, we recommend
#respectively using 8, 20000, 10000 and 0.99.

#Result
  write.table(example2,"example2.txt",row.names=F,quote=F)
  print(example2)

#        Point_Estimate HPDI_Lower HPDI_Upper
#rs001_1       1.601871  0.2058925   1.999482




##example 3:
#the ped file with header
#set "ped='ped_quantitative_header.txt'" and "ped.header=T"

#the trait is quantitative
#set "trait_type='quantitative'"

#the covariate file with header
#set "covariate='covariate_header.txt'" and "covariate.header=T"

#using uniform distribution as the prior distribution for gamma
#set "gamma_prior='uniform' "

#define the prior distributions of the coefficients of the covariates 
  covariate_prior<-rbind(c('uniform','-2,2'),c('normal','0,10'))

#define the value limits of the coefficients of the covariates
  covariate_prior_limit <- rbind(c(NA,'2'),c(NA,NA))

#set seed to get the fixed result
  set.seed(123)
  example3<-Bayes_XCI(ped='ped_quantitative_header.txt', trait_missing=NA, allele_missing=0, 
					ped.header=T, trait_type='quantitative', covariate='covariate_header.txt', 
					covariate_prior=covariate_prior, covariate_missing=NA, 
					covariate.header=T, covariate_prior_limit=covariate_prior_limit,
					gamma_prior='uniform', chains_num=2, total_sample=12000,  
					warmup_sample=2000, adapt_delta_value=0.85)

#Although the number of Markov chains, the number of iterations for each
#chain, the number of warmup iterations and the target acceptance rate
#are respectively set to 2, 12000, 2000 and 0.85 in this example, we recommend
#respectively using 8, 20000, 10000 and 0.99.

#Result
  write.table(example3,"example3.txt",row.names=F,quote=F)
  print(example3)

#        Point_Estimate HPDI_Lower HPDI_Upper
#rs001_1       1.417922  0.3589927   1.999366




##example 4:
#the ped dataframe with colnames
#set "ped=ped"
  ped<-read.table('ped_quantitative_header.txt',header=T)

#the trait is quantitative
#set "trait_type='quantitative'"

#the covariate dataframe with colnames
#set "covariate=covariate"
  covariate<-read.table('covariate_header.txt',header=T)

#using uniform distribution as the prior distribution for gamma
#set "gamma_prior='uniform'"

#define the prior distributions of the coefficients of the covariates
  covariate_prior<-rbind(c('uniform','-4,4'),c('normal','0,10'))

#if no constraint for the coefficients of the two covariates
  covariate_prior_limit <- rbind(c(NA,NA),c(NA,NA))

#set seed to get the fixed result
  set.seed(123)
  example4<-Bayes_XCI(ped=ped, trait_missing=NA, allele_missing=0, ped.header=F, trait_type='quantitative',
			covariate=covariate, covariate_prior=covariate_prior, covariate_missing=NA,
			covariate.header=F, covariate_prior_limit=covariate_prior_limit, gamma_prior='uniform',
			chains_num=2, total_sample=12000, warmup_sample=2000, adapt_delta_value=0.85)

#Although the number of Markov chains, the number of iterations for each
#chain, the number of warmup iterations and the target acceptance rate
#are respectively set to 2, 12000, 2000 and 0.85 in this example, we recommend
#respectively using 8, 20000, 10000 and 0.99.

#Result
  write.table(example4,"example4.txt",row.names=F,quote=F)
  print(example4)

#setting different prior distributions and different limits
#for the coefficients of the covariates affects the results
#        Point_Estimate HPDI_Lower HPDI_Upper
#rs001_1       1.458581  0.3624141   1.999038



###########################################################################
##################    examples for Frequen_XCI()    #######################
###########################################################################
##example 5:
#the ped file without header 
#set "ped='ped_qualitative.txt'" and "ped.header=F"

#no covariate 
#set "covariate=NULL"

#the trait is qualitative 
#set "trait_type='qualitative'"

  example5<-Frequen_XCI(ped='ped_qualitative.txt', trait_missing=9, allele_missing=0,
			ped.header=F, trait_type='qualitative',
			covariate=NULL, covariate.header=F, covariate_missing=NA)


#Result
  write.table(example5,"example5.txt",row.names=F,quote=F)
  print(example5)

#      F_Point_Estimate F_Lower F_Upper F_discontinuous PF_Point_Estimate
#snp_1                2       0       2               0          1.756293
#      PF_Lower PF_Upper
#snp_1        0        2



##example 6:
#the ped file without header 
#set "ped='ped_qualitative.txt'" and "ped.header=F"

#the covariate file with header 
#set "covariate='covariate_header.txt'" and "covariate.header=T"

#the trait is qualitative 
#set "trait_type='qualitative'"

  example6<-Frequen_XCI(ped='ped_qualitative.txt', trait_missing=9, allele_missing=0,
			ped.header=F, trait_type='qualitative', covariate='covariate_header.txt', 
			covariate.header=T, covariate_missing=NA)


#Result
  write.table(example6,"example6.txt",row.names=F,quote=F)
  print(example6)

#      F_Point_Estimate F_Lower F_Upper F_discontinuous PF_Point_Estimate
#snp_1                2       0       2               0                 2
#      PF_Lower PF_Upper
#snp_1        0        2




##example 7:
#the ped dataframe without colnames
  ped<-read.table('ped_quantitative.txt',header=F)

#the trait is quantitative 
#set "trait_type='quantitative'"

#the covariate dataframe without colnames
  covariate<-read.table('covariate.txt',header=F)

  example7<-Frequen_XCI(ped=ped, trait_missing=NA, allele_missing=0,
			ped.header=F, trait_type='quantitative',
			covariate=covariate, covariate.header=F, covariate_missing=NA)

#Result
  write.table(example7,"example7.txt",row.names=F,quote=F)
  print(example7)

#      F_Point_Estimate F_Lower F_Upper F_discontinuous PF_Point_Estimate
#snp_1         1.271208       0       2               0          1.239847
#      PF_Lower PF_Upper
#snp_1        0        2








