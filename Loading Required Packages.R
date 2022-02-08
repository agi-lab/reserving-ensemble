

###Installing packages if required
if (!require("PACLasso")) install.packages("PACLasso")
if (!require("groupdata2")) install.packages("groupdata2")
if (!require("SynthETIC")) install.packages("SynthETIC")
if (!require("gamlss.inf")) install.packages("gamlss.inf")
if (!require("gamlss")) install.packages("gamlss")
if (!require("ChainLadder")) install.packages("ChainLadder")
#install.packages("ChainLadder")
if (!require("ggplot2")) install.packages("ggplot2")




###Loading the required packages
library("PACLasso")
library("groupdata2")
library("SynthETIC")
library("gamlss.inf")
library("gamlss")
library("ChainLadder")
library("ggplot2")