###################################################
#                poiCalc.R                        #
#                                                 #
# This script will calculate four things related  #
# to the probability of identity. First, it will  #
# calculate the POI for each specific genotype    #
# (function poi.ind). Second, it will calculate   #
# the POI-Sib for each genotype (function         #
# poi.ind.sib). Third, it will calculate the POI  #
# for the population as a whole (function         #
# poi.all). And fourth, it will calculate the     #
# POI-sib for the whole population (function      #
# poi.all.sib).                                   #
#-------------------------------------------------#
# Tested by hand, and against results from        #
# Gimlet, on 04-Ocr-2016                          #
#                                                 #
#                    by                           #
#               Tim Frasier                       #
#                                                 #
#       Last updated: 04-Oct-2016                 #
###################################################




###################################################
#                   poi.all                       #
#                                                 #
# This function calculates the POI for the        #
# population, based on the equation described in  #
# Paetkau & Strobeck (1994) Molecular Ecology 3:  #
# 489-495.                                        #
###################################################

poi.all <- function(freqs) {
  
  #----------------------------------------#
  # Get the number of loci and their names #
  #----------------------------------------#
  nLoci = length(unique(freqs$Locus))
  locusNames = unique(freqs$Locus)
  
  #-----------------------------------------#
  # Get the number of alleles at each locus #
  #-----------------------------------------#
  nAlleles = rep(0, times = nLoci)
  
  for (i in 1:nLoci) {
    nAlleles[i] = sum(freqs$Locus == locusNames[i])
  }
  
  #-----------------------------------------#
  #  Calculate pi^4 for each locus          #
  #-----------------------------------------#
  locusPOI1 = rep(0, times = nLoci)
  rowCount = 1
  
  for (i in 1:nLoci) {
      tempPOI = 0
      for (j in 1:nAlleles[i]) {
        tempPOI = tempPOI + (freqs[rowCount, 3]^4)
        rowCount = rowCount + 1
      }
      locusPOI1[i] = tempPOI
  }
  
  #-----------------------------------------#
  #  Calculate (2PiPj)^2 for each locus     #
  #-----------------------------------------#
  locusPOI2 = rep(0, times = nLoci)
  countSum = 0
  
  for (i in 1:nLoci) {
    tempPOI = 0
    for (j in 1:(nAlleles[i] - 1)) {
      for (k in (j + 1):nAlleles[i]) {
        tempPOI = tempPOI + (2 * (freqs[(countSum + j), 3] * freqs[(countSum + k), 3]))^2
      }
    }
    locusPOI2[i] = tempPOI
    countSum = countSum + nAlleles[i]
  }
  
  #----------------------------------------#
  #  Combine values for total POI          #
  #----------------------------------------#
  poiTemp = rep(0, times = nLoci)
  poi = 1
  
  for (i in 1:nLoci) {
    poiTemp[i] = locusPOI1[i] + locusPOI2[i]
  }
  
  for (i in 1:nLoci) {
    poi = poi * poiTemp[i]
  }
  
  return(poi)
}

###################################################
#                   poi.all.sib                   #
#                                                 #
# This function calculates the POI among siblings #
# for the population, based on the equation       #
# described in Evett & Weir (1998) Interpreting   #
# DNA Evidence: Statistical Genetics for Forensic #
# Scientists, and Waits et al. (2001) Molecular   #
# Ecology 10: 249-256.                            #
###################################################
poi.all.sib <- function(freqs) {
  
  #----------------------------------------#
  # Get the number of loci and their names #
  #----------------------------------------#
  nLoci = length(unique(freqs$Locus))
  locusNames = unique(freqs$Locus)
  
  #-----------------------------------------#
  # Get the number of alleles at each locus #
  #-----------------------------------------#
  nAlleles = rep(0, times = nLoci)
  
  for (i in 1:nLoci) {
    nAlleles[i] = sum(freqs$Locus == locusNames[i])
  }
  
  #-----------------------------------------#
  #  Calculate pi^4 for each locus          #
  #-----------------------------------------#
  locusPOI1 = rep(0, times = nLoci)
  rowCount = 1
  
  for (i in 1:nLoci) {
    tempPOI = 0
    for (j in 1:nAlleles[i]) {
      tempPOI = tempPOI + (freqs[rowCount, 3]^4)
      rowCount = rowCount + 1
    }
    locusPOI1[i] = tempPOI
  }
  
  #-----------------------------------------#
  #  Calculate pi^2 for each locus          #
  #-----------------------------------------#
  locusPOI2 = rep(0, times = nLoci)
  rowCount = 1
  
  for (i in 1:nLoci) {
    tempPOI = 0
    for (j in 1:nAlleles[i]) {
      tempPOI = tempPOI + (freqs[rowCount, 3]^2)
      rowCount = rowCount + 1
    }
    locusPOI2[i] = tempPOI
  }
  
  #----------------------------------------#
  #  Combine values for total POI          #
  #----------------------------------------#
  poiTemp = rep(0, times = nLoci)
  poi = 1
  
  for (i in 1:nLoci) {
    poiTemp[i] = 0.25 + (0.5 * locusPOI2[i]) + (0.5 * locusPOI2[i]^2) - (0.25 * locusPOI1[i])
  }
  
  for (i in 1:nLoci) {
    poi = poi * poiTemp[i]
  }
  
  return(poi)
}


#############################################
#               poi.ind                     #
#                                           #
# This function calculates the probability  #
# of identity for individual genotypes      #
# within a population based on the          #
# calculations described in Paetkau &       #
# Strobeck (1994) Molecular Ecology 3: 489- #
# 495                                       #
#############################################

poi.ind <- function(microdata, freqs) {

  #-----------------------------------------#
  # Transform the data to its multivariate  #
  # counterpart                             #
  #-----------------------------------------#
  data2 = to_mv(microdata)
  
  #-----------------------------------------#
  # Get the number of individuals in the    #
  # genotype file, and the # of loci        #
  #-----------------------------------------#
  nInds = length(microdata[, 1])
  nLoci = length(unique(freqs$Locus))
  
  #-----------------------------------------#
  # Get the number of alleles at each locus #
  #-----------------------------------------#
  locusNames = unique(freqs$Locus)
  nAlleles = rep(0, times = nLoci)
  
  for (i in 1:nLoci) {
      nAlleles[i] = sum(freqs$Locus == locusNames[i])
  }
  

  #-------------------------------------------#
  #        Calculate the POI                  #
  #-------------------------------------------#

  
  #-------------------------------------------#
  # This first matrix goes through the        #
  # genotypes and conducts the appropriate    #
  # calculation for homozygotes, but replaces #
  # all other values with 1                   #
  #-------------------------------------------#
  poiHomo = matrix(NA, nrow = nInds, ncol = length(data2[1, ]))
  
  for (i in 1:nInds) {
      for (j in 1:length(data2[1, ])) {
          if (data2[i, j] == 1.0) {
              poiHomo[i, j] = (freqs[j, 3] * freqs[j, 3])
          } else {
              poiHomo[i, j] = 999
          }
      }
  }
  
  #----------------------------------#
  #  Multiply Across Homozygous Loci #
  #----------------------------------#
  homo_product <- rep(NA, times = nInds)
  for (i in 1:nInds) {
    holder <- 1
    for(j in 1:length(poiHomo[1, ])) {
      if (poiHomo[i, j] < 999) {
        holder <- holder * poiHomo[i, j]
      } else {
        holder <- holder * 1
      }
    }
    homo_product[i] <- holder
  }
  
  
  #------------------------------------------#
  # Get heterozygote frequencies             #
  #------------------------------------------#
  poiHetero = matrix(NA, nrow = nInds, ncol = length(data2[1, ]))
  
  for (i in 1:nInds) {
      for (j in 1:length(data2[1, ])) {
          if (data2[i, j] == 0.5) {
              poiHetero[i, j] = freqs[j, 3]
          } else {
              poiHetero[i, j] = 999
          }
      }
  }
  
  #--------------------------------------------#
  # Organize heterozygote information by locus #
  #--------------------------------------------#
  poiHetero2 <- matrix(999, nrow = nInds, ncol = (2 * nLoci))
  
  for (i in 1:nInds) {
    counter <- 1
    for (j in 1:length(poiHetero[1, ])) {
      if(poiHetero[i, j] < 999) {
        poiHetero2[i, counter] <- poiHetero[i, j]
        counter <- counter + 1
      }
    }
  }
  

    #-------------------------------------------#
    # Calculate appropriate heterozygote values #
    #-------------------------------------------#
    hetero_product <- rep(NA, times = nInds)
    
    for (i in 1:nInds) {
      holder <- 1
      for (j in seq(1, nLoci, 2)) {
        if (poiHetero2[i, j] < 999) {
          holder <- holder * (2 * poiHetero2[i, j] * poiHetero2[i, j + 1])
        }
      }
    hetero_product[i] <- holder
    }
    
    #-------------------------------------#
    # Now, multiply together to get final #
    # POI value                           #
    #-------------------------------------#
    poiFinal = homo_product * hetero_product

    #------------------------------------#
    # Combine the data to return to the  #
    # user                               #
    #------------------------------------#
    names = microdata[, 1]
    results = data.frame(names, poiFinal)
    return(results)
}  


#############################################
#               poi.ind.sib                 #
#                                           #
# This function calculates the probability  #
# of identity for individual genotypes      #
# within a population based on the          #
# sibling calculations described in Woods   #
# et al. (1999) Widlife Society Bulletin    #
# 27(3): 616-627                            #
#############################################

poi.ind.sib <- function(microdata, freqs) {
    
    #-----------------------------------------#
    # Transform the data to its multivariate  #
    # counterpart                             #
    #-----------------------------------------#
    data2 = to_mv(microdata)
    
    #-----------------------------------------#
    # Get the number of individuals in the    #
    # genotype file, and the # of loci        #
    #-----------------------------------------#
    nInds = length(microdata[, 1])
    nLoci = length(unique(freqs$Locus))
    
    #-----------------------------------------#
    # Get the number of alleles at each locus #
    #-----------------------------------------#
    locusNames = unique(freqs$Locus)
    nAlleles = rep(0, times = nLoci)
    
    for (i in 1:nLoci) {
        nAlleles[i] = sum(freqs$Locus == locusNames[i])
    }
    
    
    #-------------------------------------------#
    #        Calculate the POI                  #
    #-------------------------------------------#
    poiHomo = matrix(NA, nrow = nInds, ncol = length(data2[1, ]))
    
    #-------------------------------------------#
    # This first matrix goes through the        #
    # genotypes and conducts the appropriate    #
    # calculation for homozygotes, but replaces #
    # all other values with 1                   #
    #-------------------------------------------#
    for (i in 1:nInds) {
        for (j in 1:length(data2[1, ])) {
            if (data2[i, j] == 1.0) {
                poiHomo[i, j] = ((1 + (2 * freqs[j, 3]) + (freqs[j, 3] * freqs[j, 3])) / 4)
            } else {
                poiHomo[i, j] = 999
            }
        }
    }
    
    #----------------------------------#
    #  Multiply Across Homozygous Loci #
    #----------------------------------#
    homo_product <- rep(NA, times = nInds)
    for (i in 1:nInds) {
      holder <- 1
      for(j in 1:length(poiHomo[1, ])) {
        if (poiHomo[i, j] < 999) {
          holder <- holder * poiHomo[i, j]
        } else {
          holder <- holder * 1
        }
      }
      homo_product[i] <- holder
    }
    
    #------------------------------------------#
    # Get heterozygote frequencies             #
    #------------------------------------------#
    poiHetero = matrix(NA, nrow = nInds, ncol = length(data2[1, ]))
    
    for (i in 1:nInds) {
        for (j in 1:length(data2[1, ])) {
            if (data2[i, j] == 0.5) {
                poiHetero[i, j] = freqs[j, 3]
            } else {
                poiHetero[i, j] = 999
            }
        }
    }
    

    #--------------------------------------------#
    # Organize heterozygote information by locus #
    #--------------------------------------------#
    poiHetero2 <- matrix(999, nrow = nInds, ncol = (2 * nLoci))
    
    for (i in 1:nInds) {
      counter <- 1
      for (j in 1:length(poiHetero[1, ])) {
        if(poiHetero[i, j] < 999) {
          poiHetero2[i, counter] <- poiHetero[i, j]
          counter <- counter + 1
        }
      }
    }
    
    
    #-------------------------------------------#
    # Calculate appropriate heterozygote values #
    #-------------------------------------------#
    hetero_product <- rep(NA, times = nInds)
    
    for (i in 1:nInds) {
      holder <- 1
      for (j in seq(1, nLoci, 2)) {
        if (poiHetero2[i, j] < 999) {
          holder <- holder * ((1 + poiHetero2[i, j] + poiHetero2[i, j + 1] + (2 * poiHetero2[i, j] * poiHetero2[i, j + 1])) / 4)
        }
      }
      hetero_product[i] <- holder
    }
    
    
    #-------------------------------------#
    # Now, multiply together to get final #
    # POI-sib value                       #
    #-------------------------------------#
    poiFinal = homo_product * hetero_product
    
    #------------------------------------#
    # Combine the data to return to the  #
    # user                               #
    #------------------------------------#
    names = microdata[, 1]
    results = data.frame(names, poiFinal)
    return(results)
}


#############################################
#               poi.obs                     #
#                                           #
# This function calculates the observed     #
# probability of identity for a given data  #
# set based on incrementing the number of   #
# loci examined. Specifically, it starts    #
# with one locus and calculates the         #
# proportion of pairs with the same         #
# genotype, then adds a second locus and    #
# does the same thing, for a user-defined   #
# number of loci.                           #
#############################################

poi.obs <- function(microdata, freqs) {

    #-----------------------------------------#
    # Get the number of individuals in the    #
    # genotype file, and the # of loci        #
    #-----------------------------------------#
    nInds = length(microdata[, 1])
    nLoci = length(unique(freqs$Locus))
    
    
    #-----------------------------------------#
    # Calculate the total number of pairwise  #
    # comparisons.                            #
    #-----------------------------------------#
    nComparisons = (nInds * (nInds - 1)) / 2
    
    
    #-----------------------------------------#
    # Create array to hold results (one       #
    # matrix for each locus)                  #
    #-----------------------------------------#
    poiCompare = array(0, dim = c(nInds, nInds, nLoci))
    
    for (i in 1:(nInds - 1)) {
        
        for (j in (i + 1):nInds) {
            
            for (k in 2:(nLoci + 1)) {

                if (microdata[i, k] == microdata[j, k]) {
                    poiCompare[j, i, (k - 1)] = 1
                } else {
                    poiCompare[j, i, (k - 1)] = 0
                }
            }

        }
    }
    
    
    #----------------------------------------#
    # Calculate individuals that have the    #
    # same profile across all chosen loci    #
    #----------------------------------------#
    
    #--- Create one matrix that sums values across arrays ---#
    poiCompare2 = matrix(0, nrow = nInds, ncol = nInds)
    
    for (i in 1:nInds) {
        
        for (j in 1:nInds) {
            
            poiCompare2[i, j] = sum(poiCompare[i, j, ])
        }
    }
    
    #--- Create another matrix that converts ---#
    #--- these to 1s or 0s indicating whether --#
    #--- or not they match at all loci       ---#
    poiCompare3 = matrix(0, nrow = nInds, ncol = nInds)
    
    for (i in 1:nInds) {
        
        for (j in 1:nInds) {
            
            if (poiCompare2[i, j] == nLoci) {
                poiCompare3[i, j] = 1
            } else {
                poiCompare2[i, j] = 0
            }
        }
    }
    
    
    #---------------------------------------#
    # Calculate Observed POI                #
    #---------------------------------------#
    poiObs <- (sum(poiCompare3)) / nComparisons
    return(poiObs)
}


#############################################
#               poi.compare                 #
#                                           #
# This function compares poi values         #
# ("regular, "sibs", "observed") across     #
# different number of loci. It does this    #
# iiteratively by starting with one locus   #
# and then sequentially adding on           #
# additional loci. The output is a file of  #
# the values obtained, as well as a plot.   #
#############################################

poi.compare <- function(microdata, loci) {
    
    #------------------------------------#
    # Create vectors to hold results     #
    #------------------------------------#
    poi_norm = rep(0, times = loci)
    poi_sib = rep(0, times = loci)
    poi_obs = rep(0, times = loci)
    
    
    #------------------------------------#
    # Sequentially add loci, and conduct #
    # calculations.                      #
    #------------------------------------#
    for (l in 1:loci) {
        
        
        #--------------------------------#
        # Generated subsetted data       #
        #--------------------------------#
        subData <- microdata[, 1:(l + 1)]
        subFreqs <- frequencies(subData)
        
        #--------------------------------#
        # Conduct calculations           #
        #--------------------------------#
        poi_norm[l] = poi.all(subFreqs)
        poi_sib[l] = poi.all.sib(subFreqs)
        poi_obs[l] = poi.obs(subData, subFreqs)
    }
    
    #----------------------------------#
    # Join data and plot               #
    #----------------------------------#
    numLoci = 1:loci
    results = cbind(numLoci, poi_norm, poi_obs, poi_sib)
    
    plot(numLoci, poi_norm, ylim = c(0, max(poi_sib)), type = "b", pch = 16, xlab = "Number of Loci", ylab = "POI")
    lines(numLoci, poi_obs, ylim = c(0, max(poi_sib)), type = "b", pch = 16, col = "red")
    lines(numLoci, poi_sib, ylim = c(0, max(poi_sib)), type = "b", pch = 16, col = "blue")
    
    legend("top", legend = c("POI", "POI-Obs", "POI-Sib"), pch = c(16, 16, 16), col = c("black", "red", "blue"), bty = "n")
    
    return(results)
    
}