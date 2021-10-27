## Code by Femke Smit
# Final Project Code
# Handed in on 27-10-2021
# Course: MSB1015 - Scientific Programming 
# Teacher: Agnieszka Smolinska

##Package installation and loading
packages <- c('readr', 'ggplot2', 'vegan', 'ape', 'ecodist', 'stringr', 'sparsepca','devtools', 'BiocManager', 'robustbase', 'mrfDepth', 'splus2R', 'ggfortify', 'data.table' ,'RColorBrewer')
biomart_packages <- c('biomaRt','pcaMethods','PCAtools', 'ComplexHeatmap') 
all_packages <- c(packages, biomart_packages)

#Run the below code to install all the packages that are not yet installed on your computer
'%ni%' <- Negate('%in%')
packages_install <- packages[packages %ni% rownames(installed.packages())]
biomart_install <- biomart_packages[biomart_packages %ni% rownames(installed.packages())]
install.packages(packages_install)
if (!requireNamespace("pcaMethods", quietly = TRUE))
  BiocManager::install("pcaMethods", dependencies = TRUE)
BiocManager::install(biomart_install)

#Load all the packages
lapply(all_packages, require, character.only = TRUE)
#these packages gave me trouble later on in the file so I wanted to lead them again here just to be sure
require(PCAtools)
require(ComplexHeatmap)
require(sparsepca)
require(data.table)


##Loading the main files
file_location <- "~/Google Drive/systeem vakken/network/"     #replace with your own file location!!!

bee_flower_data <- read_csv(paste(file_location, "network.csv", sep = ""))
plant_data <- read_csv(paste(file_location, "floral_survey.csv", sep = ""))
site_meta <- read_csv(paste(file_location, "site_data.csv", sep = ""))
weather_data <- read_delim(paste(file_location, "weather_data", sep = ""), delim = ";", escape_double = FALSE, trim_ws = TRUE)
bee_meta <- read_csv(paste(file_location, "bb_traits.csv", sep = ""))
bee_meta <- bee_meta[-14,]
#I ended up not really using the plant data but the data adaptation still took quite some work so I left it in the code


##Adapting the files in various ways

#streamlined files
bee_sl <- bee_flower_data[,c("date", "site", "bb.sp", "plant.sp")]
plant_sl <- plant_data[,c("date", "site", "plant.sp")] 
weather_sl <- weather_data[,c("MESS_DATUM", "V_TE010M")]
colnames(weather_sl) <- c("date", "temp")

bee_sl$date <- str_replace_all(bee_sl$date, "-", "")
plant_sl$date <- str_replace_all(plant_sl$date, "-", "")

bee_sl$date_site <- paste(bee_sl$date, bee_sl$site)

#Finding the relevant dates for the weather data
min_date <- min(bee_sl$date)
max_date <- max(bee_sl$date)

weather_sl <- weather_sl[weather_sl$date <= max_date,]
weather_sl <- weather_sl[weather_sl$date >= min_date,]

#Creating the count files
#finding all unique bee and plant species
unique_bee_species <- unique(bee_sl$bb.sp)
unique_plant_species <- unique(plant_sl$plant.sp)
colnames_bee <- c("date", "site", unique_bee_species)
colnames_plant <- c("date", "site", unique_plant_species)

#preparing the shape of the count files
nrow_bee_count <- nrow(unique(bee_sl[,1:2]))
nrow_plant_count <- nrow(unique(plant_sl[,1:2]))

#I attempted another form for a count file for the kernelPCA but could for some reason not get this file to work
nrow_plant_count2 <- nrow_bee_count

bee_count <- data.frame(matrix(ncol = length(colnames_bee) + 2, nrow = nrow_bee_count))
colnames(bee_count) <- colnames_bee

plant_count <- data.frame(matrix(ncol = length(colnames_plant) + 2, nrow = nrow_plant_count))
colnames(plant_count) <- colnames_plant

plant_count2 <- data.frame(matrix(ncol = length(colnames_plant) + 2, nrow = nrow_plant_count2))
colnames(plant_count2) <- colnames_plant


#filling in the count files
bee_count[,1:2] <- unique(bee_sl[,1:2])
rownames(bee_count) <- paste(bee_count$date, bee_count$site)
plant_count[,1:2] <- unique(plant_sl[,1:2])
rownames(plant_count) <- paste(plant_count$date, plant_count$site)
plant_count2[,1:2] <- unique(bee_sl[,1:2])
rownames(plant_count2) <- paste(bee_count$date, bee_count$site)

#filling the count files with zeros
bee_zeros <- data.frame(matrix(0, ncol = (ncol(bee_count) - 2), nrow = nrow_bee_count))
bee_count[,3:ncol(bee_count)] <- bee_zeros
plant_zeros <- data.frame(matrix(0, (ncol = ncol(plant_count) - 2), nrow = nrow_plant_count))
plant_count[,3:ncol(plant_count)] <- plant_zeros
plant_zeros2 <- data.frame(matrix(0, (ncol = ncol(plant_count2) - 2), nrow = nrow_bee_count))
plant_count2[,3:ncol(plant_count)] <- plant_zeros2


#counting the bees and the plants
for (row in 1:nrow(bee_sl)){
  counted_bee <- as.character(bee_sl[row, 3])
  date_site <- paste(as.character(bee_sl[row, 1]), as.character(bee_sl[row, 2]))
  bee_count[date_site, counted_bee] <- bee_count[date_site, counted_bee] + 1
}

for (row in 1:nrow(plant_sl)){
  counted_plant <- as.character(plant_sl[row, 3])
  date_site <- paste(as.character(plant_sl[row, 1]), as.character(plant_sl[row, 2]))
  plant_count[date_site, counted_plant] <- plant_count[date_site, counted_plant] + 1
}

#this second plant file caused an error for me here, so I couldn't continue using it for the Kernel PCA
for (row in 1:nrow(bee_sl)){
  counted_plant <- as.character(bee_sl[row, 4])
  date_site <- paste(as.character(bee_sl[row, 1]), as.character(bee_sl[row, 2]))
  plant_count2[date_site, counted_plant] <- plant_count2[date_site, counted_plant] + 1
}


#inspecting the created bee and plant count files
View(bee_count)
View(plant_count)

#somehow the last two columns of the files are entirely empty, so lets remove those
bee_count <- bee_count[,1:(ncol(bee_count)-2)]
plant_count <- plant_count[,1:(ncol(plant_count)-2)]


#saving the dataframes so we can load them in directly later
write.csv(bee_count, paste(file_location, "bee_count.csv"), row.names = FALSE)
write.csv(plant_count, paste(file_location, "plant_count.csv"), row.names = FALSE)

#load the count files if you've created them before
bee_count <- read.csv(paste(file_location, "bee_count.csv", sep = ""))
plant_count <- read.csv(paste(file_location, "plant_count.csv", sep = ""))

##Completing the metadata and creating files that include both metadata and counts
#creating the temperature difference column for the sites
mean_temps <- mean(site_meta$temp.mean)
site_meta$temp.diff <- site_meta$temp.mean - mean_temps

#matching weather data to site-day measurements
rownames(weather_sl) <- weather_sl$date
bee_sl_backup <- bee_sl

bee_count_temps <- merge(bee_count, weather_sl, by = "date")
bee_count_temps<- merge(bee_count_temps, site_meta[,c("site","temp.diff"), by = "site"])
bee_count_temps$local_temp <- bee_count_temps$temp + bee_count_temps$temp.diff

plant_count_temps <- merge(plant_count, weather_sl, by = "date")
plant_count_temps<- merge(plant_count_temps, site_meta[,c("site","temp.diff"), by = "site"])
plant_count_temps$local_temp <- plant_count_temps$temp + plant_count_temps$temp.diff

#making bee and count files that contain all relevant info
bee_total <- bee_count_temps
plant_total <- plant_count_temps

bee_total <- merge(bee_total, site_meta[,c("site", "elev.mean", "management", "transect", "elev.class","elev.class2")])
plant_total <- merge(plant_total, site_meta[,c("site", "elev.mean", "management", "transect", "elev.class","elev.class2")])

date_site <- paste(as.character(bee_total$date), as.character(bee_total$site))
rownames(bee_total) <- date_site
bee_total$date_site <- date_site

bee_total_with_outliers <- bee_total[,-18] #removing humi since only counted 3 times
bee_total <- bee_total[-649,-18]  #removing 649 since outlier in PCA

bee_count_only <- bee_total[,3:17] 
bee_meta_only <- bee_total[,c(1,2,18:25)]
bee_meta_only$date_site <- rownames(bee_meta_only)


#creating a summary file for each bee species
bee_summary <- data.frame(matrix(0, ncol = length(unique(bee_total$site)), nrow = ncol(bee_count_only)))
colnames(bee_summary) <- c(unique(bee_total$site))
rownames(bee_summary) <- colnames(bee_count_only)

for (row in 1:nrow(bee_summary)){
  for (col in 1:nrow(bee_summary)){
    site <- colnames(bee_summary)[col]
    bee <- rownames(bee_summary)[row]
    bee_summary[row,col] <- sum(bee_total[bee_total$site == site, bee])
  }
}

bee_summary$total <- rowSums(bee_summary)

#adding elevation to the summary as that came out as important in the PCA (found later in code, but most clear to add here)
bee_summary_elev <- data.frame(t(bee_summary))
bee_summary_elev$site <- rownames(bee_summary_elev)
bee_summary_elev <-  merge(bee_summary_elev, site_meta[c("elev.mean","site")], by = "site")
#removing sites not used
bee_summary_elev <- bee_summary_elev[1:16,]


#Exploratory plots:
bee_count_only_nas <- bee_count_only
bee_count_only_nas[bee_count_only == 0] <- NA
boxplot(bee_count_only_nas, ylab = "Amount counted per sample", xlab = "Bee species") #this will only show the sample values > 0 to give an indication of the value distribution 

#making a violin plot of bee counts in both different mean elevation levels and different local temperatures
bee_sl_ev <- merge(bee_sl, site_meta[c("elev.mean","site", "temp.mean")], by = "site")
bee_sl_ev$date_site <- paste(as.character(bee_sl_ev$date), as.character(bee_sl_ev$site))
bee_sl_ev <- merge(bee_sl_ev, bee_meta_only[c("date_site","temp", "temp.diff", "local_temp", "management", "transect", "elev.class", "elev.class2")], by = "date_site")

ggplot(bee_sl_ev[!bee_sl_ev$bb.sp=="humi",], aes(x=bb.sp, y=elev.mean)) +     #not using humi cause only 2 counts in total
  geom_violin() + geom_boxplot(width=0.1) + theme_minimal() + labs(x = "Bee species", y = "Mean elevation")
ggplot(bee_sl_ev[!bee_sl_ev$bb.sp=="humi",], aes(x=bb.sp, y=local_temp)) +     #not using humi cause only 2 counts in total
  geom_violin() + geom_boxplot(width=0.1) + theme_minimal() + labs(x ="Bee species", y =  "Local temperature")

##Transforming the data using an inverse hyperbolic sine transformation, to make further analysis possible
bee_ihs <- apply(bee_count_only, 2, sinh)
#plant_ihs <- apply(plants_bees[,29:382], 2, sinh)
ihs_with_outliers <- apply(bee_total_with_outliers[3:17], 2, sinh)

#making ihs transformed files with the metadata attached (for future analysis)
ihs_outliers_total <- bee_total_with_outliers
ihs_outliers_total[,3:17] <- ihs_with_outliers

bee_total_ihs <- bee_total
bee_total_ihs[,3:17] <- bee_ihs



##PCA
pca_ihs <- PCAtools::pca(t(bee_ihs), metadata = bee_meta_only, rank = 10, scale = TRUE)

#Performing a PCA on the file with the outliers for completeness as I of course did this to determine what the outliers are
pca_ihs_outliers <- PCAtools::pca(t(ihs_outliers_total[,3:17]), metadata = ihs_outliers_total[,c(1,2,18:ncol(ihs_outliers_total))], rank = 10, scale = TRUE)  
PCAtools::biplot(pca_ihs_outliers, x = "PC1", y= "PC2", showLoadings = FALSE, colby = "elev.mean",lab = rownames(pca_ihs_outliers$metadata), legendPosition = "right")

#Inspecting the proper PCA output
summary(pca_ihs)
PCAtools::screeplot(pca_ihs)

#I made many more pairsplots colored by many covariates but I didn't want to bog down the code file with all of them. The colored by elevation turned out to be most interesting.
pca_pairsplot_elev <- pairsplot(pca_ihs,
                                colby = 'elev.mean',
                                title = 'Pairs plot, coloured by mean elevation',
                                titleLabSize = '20')
pca_pairsplot_elev

#some biplots/loading plots to more clearly show the relationship between the covariates and the PCs
biplot(pca_ihs, x = "PC1", y= "PC2", showLoadings = TRUE, colby = "elev.mean",lab = rownames(pca_ihs$metadata), legendPosition = "right")
biplot(pca_ihs, x = "PC1", y= "PC2", showLoadings = FALSE, colby = "elev.class",lab = rownames(pca_ihs$metadata), legendPosition = "right")

#correlation plot PCs
eigencorplot(pca_ihs,
             components = getComponents(pca_ihs, 1:10),
             metavars = colnames(bee_meta_only)[1:9],
             col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
             cexCorval = 0.7,
             colCorval = 'white',
             fontCorval = 2,
             posLab = 'bottomleft',
             rotLabX = 45,
             posColKey = 'top',
             cexLabColKey = 1.5,
             scale = TRUE,
             main = 'PC1-10 correlations',
             colFrame = 'white',
             plotRsquared = FALSE)


#function for autoscaling
autoscaling <- function(dat){
  m <- colMeans(dat)
  s <- apply(dat,2,sd)
  res <- vector('list',3)
  
  dat2 <- sweep(dat,2,m,'-')
  dat2 <- sweep(dat2,2,s,'/')
  res[[1]] <- dat2
  res[[2]] <- m
  res[[3]] <- s
  names(res) <- c('matrix','mean','sd')
  return(res)
}

bee_ihs_auto <- autoscaling(bee_ihs)[[1]]
rownames(bee_ihs_auto) <- rownames(bee_ihs)
colnames(bee_ihs_auto) <- colnames(bee_ihs)


#Making a heatmap including dendograms
# Create the heatmap annotation
col = list(clade = c("SF" = "darkgreen", "K" = "blue", "LF" = "lightgreen", "M" = "lightblue"), tongue_length = c("short" = "pink", "med" = "red", "long" = "purple"))
ha <- HeatmapAnnotation(
  clade = bee_meta$clade,
  tongue_length = bee_meta$pbl.w.class2,
  col = col
)

col_samples <- list(elevation_class = c("mitte" = "orange", "oben" = "yellow", "unten" = "red"), management = c("grasing" = "white", "mowing" = "gray", "none" = "black"))
ha_samples <- rowAnnotation(
  elevation_class = bee_total$elev.class,
  management = bee_total$management,
  col = col_samples
)

Heatmap(bee_ihs_auto, 
        name = "bee_count", #title of legend
        column_title = "Bee species", row_title = "Samples",
        row_names_gp = gpar(fontsize = 1), # Text size for row names
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "ward.D",
        clustering_method_columns = "ward.D",
        #top_annotation = ha,
        show_heatmap_legend = TRUE,
        top_annotation = ha,
        left_annotation = ha_samples
)


#robust sparce pca
#using autoscaled data
robspca <- robspca(bee_ihs_auto, center = FALSE)
summary(robspca)
scores <- data.frame(robspca$scores)
scores_meta <- cbind(scores, bee_meta_only)
ggplot(scores_meta, aes(X1,X2, color = elev.mean)) + 
  geom_point(shape = 16, size = 3, show.legend = TRUE, alpha = .5) +
  theme_minimal() +
  scale_color_gradient(low = "#0091ff", high = "#f0650e") +
  labs(x="PC1",y="PC2")

#using non-autoscaled data
robspca2 <- robspca(bee_ihs, scale = TRUE)
scores2 <- data.frame(robspca2$scores)
scores_meta2 <- cbind(scores2, bee_meta_only)
ggplot(scores_meta2, aes(X1,X2, color = elev.mean)) + 
  geom_point(shape = 16, size = 3, show.legend = TRUE, alpha = .5) +
  theme_minimal() +
  scale_color_gradient(low = "#0091ff", high = "#f0650e") +
  labs(x="PC1",y="PC2")



##Kernel PCA (failed)
#I wanted to include my attempt at a kernel PCA combining plant and bee data, even though it failed
plant_count$date_site <- paste(plant_count$date, plant_count$site)

#Not all the dates in the bee_count data matched the dates in the plant_count data, therefore this needed to be matched up
#finding which dates in the plant data don't match the dates in the bee_data
bee_not_plant <- setdiff(bee_count$date, plant_count$date)
a <- data.table(Value = plant_count$date)
b <- data.table(Value = bee_not_plant)
a[,merge:=Value]
b[,merge:=Value]
setkeyv(a,c('merge'))
setkeyv(b,c('merge'))
Merge_a_b <- a[b,roll='nearest']

#Giving the non-matching plant sample dates the dates of the closest matching bee sample dates
plant_count$date_matched <- plant_count$date
for (row in 1:nrow(Merge_a_b)){
  old_date <- as.numeric(Merge_a_b[row,1])
  new_date <- as.numeric(Merge_a_b[row,2])
  plant_count$date_matched[as.numeric(plant_count$date_matched) == old_date] <- new_date
}

plant_count$date_site <- paste(plant_count$date_matched, plant_count$site)
plants_bees <- merge(plant_count, bee_total, by = "date_site")
#the problem is that I couldn't get the combined file to have the same amount of samples as the bee file, it remained smaller

#this was the code I wrote for the final kernel pca, though I couldn't use it in the end
bee_kernel <- bee_ihs%*%t(bee_ihs)
plant_kernel <- plant_ihs%*%t(plant_ihs)
list_both <- list(bee_kernel, plant_kernel)
both_kernel <- Reduce("+", list_both)/2
#then, a normal PCA could've been performed on the "both_kernel" dataframe

