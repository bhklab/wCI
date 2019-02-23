# mCI
Modified Concordance Index for pharmacological readout

How to install: At present there are 2 ways to install:

* Download zip file from mCI github page. Unzip it. and install 
    ```R
    library(devtools) 
    devtools::install("mCI-master")
    ```
 
 * Using install_github : To do this you will required devtools. Install devtools in R as: 
 	```R
      install.packages("devtools")
      library(devtools) #load library
      ```
      This is privet reposeatry, you will required a a personal access tokens on github. Go to the url https://github.com/settings/tokens and login to your account. Type a Token description and click on Generate Token. Copy the token text. Now in R:
      ```R
      devtools::install_github("bhklab/mCI", auth_token = "your access token")
      ```
