# wCI
Modified Concordance Index for pharmacological readout

How to install: At present there are 2 ways to install:
First, you will need devtools. Install devtools in R as:
```R
install.packages("devtools")
library(devtools) #load library
```

1. Download zip file from mCI github page. Unzip it. and install 
```R
devtools::install("mCI-mCI_Chris")
```
 
2. Using install_github :  
```R
devtools::install_github("bhklab/mCI", ref="mCI_Chris")
```

3. This is private repository, you will require a personal access token on GitHub. Go to the url https://github.com/settings/tokens and login to your account. Type a token description and click on Generate Token. Copy the token text. Now in R:
```R
devtools::install_github("bhklab/mCI", ref="mCI_Chris", auth_token = "your access token")
```
