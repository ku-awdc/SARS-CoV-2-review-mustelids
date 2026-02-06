#####Load all packages#########################################################
library(readxl)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(stringr)
library(nlme)
library(viridis)
library(stringr)
library(writexl)
library(purrr)

#####READ IN FILES AND CREATE DATASETS##########################################
## Read in the SLR file
readings <- read_xlsx(
  'Boklund 2026 Suppl data 1 SLR Data extracted from 47 articles on SARS-CoV-2 in mink and ferrets.xlsx',
  sheet = 'Disease development',
  range = 'A1:BL168',
  col_types = c(rep("guess", 58), "date", "date", "date",rep("guess", 3))  # BG = 59, BH = 60, BI=61
)

## Transforming data
## Give some variables easier names and make some numeric
readings$logNasalPCR <- as.numeric(readings$`Nasal, MaxVirusExcreted_PCR_transformeret`)

## ND is put to 0
readings$`Oral/Throat, MaxVirusExcreted_PCR_transformeret` <- ifelse(readings$`Oral/Throat, MaxVirusExcreted_PCR_transformeret`=="ND",0,
                                                                     readings$`Oral/Throat, MaxVirusExcreted_PCR_transformeret`)
readings$logOralPCR <- as.numeric(readings$`Oral/Throat, MaxVirusExcreted_PCR_transformeret`)

readings$logNasalVirus <- as.numeric(readings$`Nasal, MaxVirusExcreted_INFvirus_transformeret`)
readings$logOralVirus  <- as.numeric(readings$`Oral/Throat, MaxVirusExcreted_INFvirus_transformeret`)


## Renaming some variables
readings$Viral_load <- readings$`Viral load`
readings$durSympt   <- readings$`Duration of symptoms`
readings$firstRNA   <- readings$dpi_RNAdetection
readings$peakRNA    <- readings$dpi_PeakRNA
readings$lastRNA    <- readings$`end_of_Shedding_RNA (last day pos)`
readings$firstVirus <- readings$dpi_infectious_virus
readings$peakVirus  <- readings$dpi_PeakVirus
readings$lastVirus  <- readings$end_of_Shedding_virus

## Catching one extra omicron
readings$Virus_cat[which(readings$Title=='Experimental Infection of Mink with SARS-COV-2 Omicron Variant and Subsequent Clinical Disease')] <- 'Omicron'

## Defining routes as either inoculation or transmission, including all types
readings$route <- NA

DirContact <- c("Direct","1st contact","2nd contact","Direct contact")

IndContact <- c("Indirect","Indirect >1m","Indirect contact (10 cm)","Indirect contact")

readings$route <- ifelse(readings$`TransmissionRoute (direct, neighbour, airborne)`%in%c("Young","Aged","Inoculated"),
                             'Inoculated',
                      ifelse(readings$`TransmissionRoute (direct, neighbour, airborne)`%in%DirContact,'Direct contact',
                      ifelse(readings$`TransmissionRoute (direct, neighbour, airborne)`%in%IndContact,'Indirect contact',
                         ifelse(is.na(readings$`TransmissionRoute (direct, neighbour, airborne)`),NA,'Direct infection'))))

## Removes those with "Direct infection"
readings <- readings[which(readings$route!='Direct infection'),]

## Cleaining and making all numerical
readings$dpi_firstSymptoms[which(readings$dpi_firstSymptoms!="ND")] <- as.numeric(readings$dpi_firstSymptoms[which(readings$dpi_firstSymptoms!="ND")])
readings$durSympt[which(readings$durSympt!="ND")]   <- as.numeric(readings$durSympt[which(readings$durSympt!="ND")]) 
readings$firstRNA   <- as.numeric(readings$firstRNA)
readings$peakRNA[which(readings$peakRNA=="5 or 8")] <-"5"
readings$peakRNA    <- as.numeric(readings$peakRNA)
readings$lastRNA[which(readings$lastRNA==">7")] <-"7"
readings$lastRNA    <- as.numeric(readings$lastRNA)
readings$firstVirus <- as.numeric(readings$firstVirus)
readings$peakVirus[which(readings$peakVirus=="5 or 1")] <-"5"
readings$peakVirus  <- as.numeric(readings$peakVirus)
readings$lastVirus  <- as.numeric(readings$lastVirus)

## Categorizing viral load inoculated
readings$ViralLoad1 <- NA
readings$ViralLoad1[which(readings$route=="Inoculated")] <- str_extract(readings$Viral_load[which(readings$route=="Inoculated")],
                                                                    "\\^(\\d+(\\.\\d+)?)") %>%
  str_remove("\\^")  # fjern ^

readings$ViralLoad1[which(readings$route=="Inoculated"&readings$Viral_load%in%c("6.0 log10 PFU","2.1×106 TCID50"))]   <- 6
readings$ViralLoad1[which(readings$route=="Inoculated"&readings$Viral_load%in%c("5.0 log10 TCID50/mL","10⁵ TCID50","1.5×105 TCID50"))] <- 5

readings$ViralLoad2 <- as.numeric(readings$ViralLoad1)
m <- max(readings$ViralLoad2[which(!is.na(readings$ViralLoad2))])
# Split inoculated viral load in groups
readings$ViralLoad3 <- NA
readings$ViralLoad3[which(!is.na(readings$ViralLoad2))] <- cut(readings$ViralLoad2[which(!is.na(readings$ViralLoad2))], breaks = seq(0, max(m), by = 1), include.lowest = TRUE, right = TRUE)

readings$ViralLoad <- NA
readings$ViralLoad <- ifelse(readings$ViralLoad3<5,"<10^5",
                           ifelse(readings$ViralLoad3==5,"10^5",">10^5"))

## Removed 8 observations, where animals were re-infected, 
## and 1 observation which was contact animals to animals infected by droplets. 
readings <- readings[which(!readings$Group%in%c("h","j","k","l","m")),]

## Creates smaller datasets; 

## Observations meassuring RNA by PCR
unitList1 <- c("Log10 copies/mL","Log10TCID50","Log10TCID50/mL")
logPCR <- readings[which(readings$PCR_unit_transformeret%in%unitList1),]

## Observations meassuring infectious virus by Log10TCID50
unitList2 <- c("Log10TCID50/mL")
logVirus <- readings[which(readings$INFvirus_unit_t%in%unitList2),]

## Creating a small dataset for comparison of RNA in nasal wash, nasal swab and oral swab

compRNA <- logPCR %>%
  filter(!is.na(logNasalPCR)) %>%
  select(Title,route,Species,logNasalPCR,logOralPCR,Nasal_excr_measured_in)

compRNA$type <- ifelse(compRNA$Nasal_excr_measured_in=="Swabs","NasalSwabs","NasalWash")
compRNA$logPCR <- compRNA$logNasalPCR

compRNA2 <- logPCR %>%
  filter(!is.na(logOralPCR)) %>%
  select(Title,route,Species,logNasalPCR,logOralPCR,Nasal_excr_measured_in)

compRNA2$type <- "OralSwabs"
compRNA2$logPCR <- compRNA2$logOralPCR

compRNA3 <- rbind(compRNA,compRNA2)

PCR_nasal_oral <- logPCR %>%
  filter(!is.na(logNasalPCR), !is.na(logOralPCR))

## A small dataset for comparison of inf virus in nasal wash, nasal swab and oral swab
compVirus <- logVirus %>%
  filter(!is.na(logNasalVirus)) %>%
  select(Title,route,Species,logNasalVirus,logOralVirus,Nasal_excr_measured_in)

compVirus$type <- ifelse(compVirus$Nasal_excr_measured_in=="Swabs","NasalSwabs","NasalWash")
compVirus$logVirus <- compVirus$logNasalVirus

compVirus2 <- logVirus %>%
  filter(!is.na(logOralVirus)) %>%
  select(Title,route,Species,logNasalVirus,logOralVirus,Nasal_excr_measured_in)

compVirus2$type <- "OralSwabs"
compVirus2$logVirus <- compVirus2$logOralVirus

compVirus3 <- rbind(compVirus,compVirus2)

Virus_nasal_oral <- logVirus %>%
  filter(!is.na(logNasalVirus), !is.na(logOralVirus))


## For days of first, peak and last, we do not need to differentiate in 
## which unit is used for excretions, and we therefor use another subsample of data
## not excluding based on unit

## Making dys to symptoms numeric, after we have looked at how many registered NO sympt.
readings$dpi_firstSymptoms <- as.numeric(readings$dpi_firstSymptoms)
readings$durSympt          <- as.numeric(readings$durSympt) 

sympt2 <- readings %>%
  filter(`Experimental/Field`=="Experimentally") %>%
            select(Title,Virus_cat,ViralLoad,route,
                   dpi_firstSymptoms, durSympt,
                   firstRNA, peakRNA,lastRNA,
                   firstVirus, peakVirus, lastVirus,Species)

## Creating dataset with only observations measuring by PCR AND inf virus in same animals
duplex <- readings %>%
  filter((!is.na(logOralPCR)&!is.na(logOralVirus))|
           (!is.na(logNasalPCR)&!is.na(logNasalVirus)))%>%
  select(Title,logNasalPCR,logOralPCR,logNasalVirus,logOralVirus)

unique(duplex$Title)
duplex$title <- NA
duplex$title[which(duplex$Title=="Severe acute respiratory disease in American mink experimentally infected with SARS-CoV-2")] <- "A"
duplex$title[which(duplex$Title=="Intranasal fusion inhibitory lipopeptide prevents direct-contact SARS-CoV-2 transmission in ferrets")]  <- "B"
duplex$title[which(duplex$Title=="Anti-SARS-CoV-2 antibodies in a nasal spray efficiently block viral transmission between ferrets")] <- "C"
duplex$title[which(duplex$Title=="Critical role of neutralizing antibody for SARS-CoV-2 reinfection and transmission")] <- "D"
duplex$title[which(duplex$Title=="Infection and rapid transmission of SARS-CoV-2 in ferrets")] <-"E"
duplex$title[which(duplex$Title=="Transmission and Protection against Reinfection in the Ferret Model with the SARS-CoV-2 USA-WA1/2020 Reference Isolate")] <-"F"
duplex$title[which(duplex$Title=="Development of an Inactivated Vaccine against SARS CoV-2")] <- "G"
duplex$title[which(duplex$Title=="Hamster and ferret experimental infection with intranasal low dose of a single strain of SARS- CoV-2")] <- "H"
duplex$title[which(duplex$Title=="SARS-CoV-2 is transmitted via contact and via the air between ferrets")] <- "I"
duplex$title[which(duplex$Title=="Therapeutic effect of CT-P59 against SARS-CoV-2 South African variant")] <- "J"


#####DESCRIPTION OF DATASETS###################################################

## Then we start plotting and analysing
## First, lets make an overview of what we want on which datasets

## logPCR         - overview of how RNA is measured in oral and nasal samples,
##                - comparing virus type in nasal by RNA
## PCR_nasal_oral - comparing nasal versus oral samples by RNA when both are sampled
## compRNA3       - comparing RNA in nasal swabs, nasal wash and oral swabs

## logVirus       - overview of how virus is measured in oral and nasal samples,
##                - comparing virus type in nasal by inf virus
## Virus_nasal_oral - comparing nasal versus oral samples by RNA when both are sampled
## compVirus3       - comparing RNA in nasal swabs, nasal wash and oral swabs

## duplex         - Comparing infectious virus versus PCR from articles including both, 24 obs

## sympt2         - comparing dpi of first, peak and last by RNA and in virus
##                - comparing grouped by inoculated/transmission and virus type

## sympt3         - comparing dpi of first, peak and last by RNA and in virus
##                - comparing grouped by inoculation load

#####Generating a plot-list#####################################################
## To save all plots in one list
all_plot_lists <- list()
#####DESCR ANALYSES OF DAYS (DPI) OF EXCRETION##################################

## Days of excretion

list_comparisons  <- c(2,4,13)

y_max_vec <- c(rep(max(sympt2[,c(5)],sympt2[,c(6)],na.rm=T),times=2),
               rep(max(sympt2[,c(7:12)],na.rm=T),times=6))

y_axis_titles <- c("First Symptoms (days)", "Duration of Symptoms (days)",
                   "First day of RNA pos. (days)", "Peak day of RNA pos. (days)", "Last day of RNA pos. (days)", 
                   "First day of infectious virus pos. (days)", "Peak day of infectious virus pos. (days)", "Last day of infectious virus pos. (days)")

Virus_cat_colors <- c(
  "Index like" = "#0D0887FF",
  "Alpha"      = "#6A00A8FF",
  "Beta"       = "#B12A90FF",
  "Delta"      = "#E16462FF",
  "Gamma"      = "#FCA636FF",
  "Omicron"    = "#F0F921FF"
)

route_colors <- c(
  "Inoculated"        = "#46337EFF",
  "Direct contact"    = "#31688EFF",
  "Indirect contact"  = "#35B779FF" 
)

Species_colors <- c(
  "Ferret"  = "#4777EFFF",
  "Mink"    = "#62FC6BFF"
)


for(i in 1:length(list_comparisons)){
  plot_for_loop <- function(dataset,x_var,y_var,y_max,y_label){
    
    x_levels_list <- list(
      "Virus_cat" = c("Index like", "Alpha", "Beta", "Delta", "Gamma", "Omicron"),
      "route" = c("Inoculated", "Direct contact", "Indirect contact"),
      "Species" = c("Ferret","Mink")
    )
    
   if (x_var %in% names(x_levels_list)) {
      dataset[[x_var]] <- factor(dataset[[x_var]], levels = x_levels_list[[x_var]])
   }
    
    ifelse(i==1,color_choice <- Virus_cat_colors,ifelse(i==2,color_choice <- route_colors,color_choice <- Species_colors))
    
    ## Filter only rows where y_var is not NA 
    valid_rows <- !is.na(dataset[[y_var]]) 
    ## Calculate counts for there rows 
    counts <- table(dataset[[x_var]][valid_rows]) 
    ## Make sure that counts follow the factor levels
    x_levels <- levels(dataset[[x_var]]) 
    new_labels <- paste0(x_levels, " (n=", counts[x_levels], ")") 
    names(new_labels) <- x_levels
    
   ggplot(dataset, aes(x = .data[[x_var]], y = .data[[y_var]], fill= .data[[x_var]])) + 
      geom_boxplot() +
      scale_fill_manual(values = color_choice, drop = FALSE) +
      geom_jitter(aes(color = Species, shape = Species), size = 0.8, alpha = 0.9) +
      
      scale_shape_manual(values = c("Ferret" = 16, "Mink" = 17)) +  # 16 = cirkel, 17 = trekant
      scale_color_manual(values = c("Ferret" ="#4777EFFF", "Mink" = "#62FC6BFF")) +
      
     theme_bw() +
     theme(
        legend.position="none",
        plot.title = element_text(size=11)
      ) +
      xlab("") +
      ylab(y_label) +
      scale_y_continuous(limits=c(0, y_max+1))+
      scale_x_discrete(labels = new_labels, guide = guide_axis(angle = 45))
  }
  
  plot_list <- colnames(sympt2[,c(5:12)]) %>%
    set_names( ~ paste0(colnames(sympt2)[list_comparisons[i]], "_", .x)) %>%
    map( ~ {
      y_index <- match(.x, colnames(sympt2[, 5:12]))
      p <- plot_for_loop(sympt2, colnames(sympt2)[list_comparisons[i]], .x, 
                    y_max_vec[y_index],
                    y_axis_titles[y_index])
      })
  
  all_plot_lists <- append(all_plot_lists, list(plot_list))
  
  for(m in 1:length(colnames(sympt2[,c(5:12)]))){
    groups <- unique(sympt2[[colnames(sympt2)[list_comparisons[i]]]])
    
    tableM <- as.data.frame(matrix(ncol = 9, nrow = 0))
    
    for(k in 1:length(groups)){
      select  <- cbind(sympt2[[colnames(sympt2[,c(5:12)])[m]]],sympt2[[colnames(sympt2[list_comparisons[i]])]])
      
      tableM <- rbind(tableM,
                      c(groups[k],
                        quantile(as.numeric(select[which(select[,2]==groups[k]),1]),c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=T),
                        sum(!is.na(select[which(select[,2]==groups[k]),1]))))
      colnames(tableM) <- c("Category","0","0.05","0.25","0.5","0.75","0.95","1","n")
    }
    write_xlsx(tableM,
               paste(colnames(sympt2[,c(5:12)])[m],
                     colnames(sympt2)[list_comparisons[i]],
                     '.xlsx',
                     sep=""))
    
  }
}

#####DESCR ANALYSES OF DAYS (DPI) OF EXCRETION, viral_load######################

## For comparison of inoculated viral load, we only use inoculated, and therefore reduce sympt2 to sympt3

sympt3 <- sympt2[which(sympt2$route=="Inoculated"),]

list_comparisons  <- c(3)

y_max_vec <- c(rep(max(sympt3[,c(5)],sympt3[,c(6)],na.rm=T),times=2),
               rep(max(sympt3[,c(7:12)],na.rm=T),times=6))

y_axis_titles <- c("First Symptoms (days)", "Duration of Symptoms (days)",
                   "First day of RNA pos. (days)", "Peak day of RNA pos. (days)", "Last day of RNA pos. (days)", 
                   "First day of infectious virus pos. (days)", "Peak day of infectious virus pos. (days)", "Last day of infectious virus pos. (days)")

ViralLoad_colors <- c(
  "<10^5"  = "#413D7BFF",
  "10^5"   = "#40B7ADFF",
  ">10^5"  = "#8AD9B1FF"
  )

for(i in 1:length(list_comparisons)){
  plot_for_loop <- function(dataset,x_var,y_var,y_max,y_label){
    
    x_levels_list <- list(
      "route" = c("Inoculated","Direct contact", "Indirect contact"),
      "Virus_cat" = c("Index like", "Alpha", "Beta", "Delta", "Gamma", "Omicron"),
      "ViralLoad" = c("<10^5", "10^5", ">10^5")
    )
    
    if (x_var %in% names(x_levels_list)) {
      dataset[[x_var]] <- factor(dataset[[x_var]], levels = x_levels_list[[x_var]])
    }
    
    ## Filter only rows where y_var is not NA 
    valid_rows <- !is.na(dataset[[y_var]]) 
    ## Calculate counts for there rows 
    counts <- table(dataset[[x_var]][valid_rows]) 
    ## Make sure that counts follow the factor levels
    x_levels <- levels(dataset[[x_var]]) 
    new_labels <- paste0(x_levels, " (n=", counts[x_levels], ")") 
    names(new_labels) <- x_levels
    
    ggplot(dataset, aes(x = .data[[x_var]], y = .data[[y_var]], fill= .data[[x_var]])) + 
      geom_boxplot() +
      
      scale_fill_manual(values = ViralLoad_colors, drop = FALSE) +
      geom_jitter(aes(color = Species, shape = Species), size = 0.8, alpha = 0.9) +
      
      scale_shape_manual(values = c("Ferret" = 16, "Mink" = 17)) +  # 16 = cirkel, 17 = trekant
      scale_color_manual(values = c("Ferret" ="#4777EFFF", "Mink" = "#62FC6BFF")) +
      
      theme_bw() +
      theme(
        legend.position="none",
        plot.title = element_text(size=11)
      ) +
      ggtitle(" ") +
      xlab("") +
      ylab(y_label) + 
      scale_y_continuous(limits=c(0, y_max+1))+
      scale_x_discrete(labels = new_labels, guide = guide_axis(angle = 45)) 
  }
  
  plot_list <- colnames(sympt3[,c(5:12)]) %>% 
    set_names( ~ paste0(colnames(sympt3)[list_comparisons[i]], "_", .x)) %>%
    map( ~ {
      y_index <- match(.x, colnames(sympt3[, 5:12]))
      plot_for_loop(sympt3, colnames(sympt3)[list_comparisons[i]], .x,
                    y_max_vec[y_index],
                    y_axis_titles[y_index])
    })
  
  all_plot_lists <- append(all_plot_lists, list(plot_list))
  
  for(m in 1:length(colnames(sympt3[,c(5:12)]))){
    groups <- unique(sympt3[[colnames(sympt3)[list_comparisons[i]]]])
    
    tableM <- as.data.frame(matrix(ncol = 9, nrow = 0))
    
    for(k in 1:length(groups)){
      select  <- cbind(sympt3[[colnames(sympt3[,c(5:12)])[m]]],sympt3[[colnames(sympt3[list_comparisons[i]])]])
      
      tableM <- rbind(tableM,
                      c(groups[k],
                        quantile(as.numeric(select[which(select[,2]==groups[k]),1]),c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=T),
                        sum(!is.na(select[which(select[,2]==groups[k]),1]))))
      colnames(tableM) <- c("Category","0","0.05","0.25","0.5","0.75","0.95","1","n")
    }
    write_xlsx(tableM,
               paste(colnames(sympt3[,c(5:12)])[m],
                     colnames(sympt3)[list_comparisons[i]],
                     '.xlsx',
                     sep=""))
    
  }
  
}

#####DESCR ANALYSES OF AMOUNTS EXCRETED#######################################
## List of datasets where we compare amount  of excreted virus

listDF         <- c('logPCR','logPCR','compRNA3',
                    'logVirus','logVirus','compVirus3')

y_names        <- c('logNasalPCR','logNasalPCR','logPCR',
                    'logNasalVirus','logNasalVirus','logVirus')

x_cat          <- c('Virus_cat','route','Species',
                    'ViralLoad','type')


Virus_cat_colors <- c(
  "Index like" = "#0D0887FF",
  "Alpha"      = "#6A00A8FF",
  "Beta"       = "#B12A90FF",
  "Delta"      = "#E16462FF",
  "Gamma"      = "#FCA636FF",
  "Omicron"    = "#F0F921FF"
)

route_colors <- c(
  "Inoculated"        = "#46337EFF",
  "Direct contact"    = "#31688EFF",
  "Indirect contact"  = "#35B779FF" 
)

Species_colors <- c(
  "Ferret"  = "#4777EFFF",
  "Mink"    = "#62FC6BFF"
)

type_colors <- c(
  "NasalSwabs"  = "#8B0AA5FF",
  "NasalWash"   = "#B93289FF",
  "OralSwabs"   = "#DB5C68FF"
)

for(i in 1:6){
  df1 <- get(listDF[i])
  
  df  <-  if (i %in% c(1, 4)) {df <- df1[which(!is.na(df1[colnames(df1[y_names[i]])])),]} else {
          if (i %in% c(2, 5)) {df <- df1[which(!is.na(df1[colnames(df1[y_names[i]])])&df1[colnames(df1[x_cat[2]])]=='Inoculated'),]} 
                               else {df <- df1}}
  
  list_comparisons  <-  if (i %in% c(1, 4)) {list_comparisons <- x_cat[1:3]} else {
                        if (i %in% c(2, 5)) {list_comparisons <- x_cat[4]} else {
                                             list_comparisons <- x_cat[5]}}
  
  list_measure      <-  y_names[[i]]
  
  y_axis_titles <- if (i %in% 1:3) {
      if (all(list_comparisons %in% c("Virus_cat", "route", "Species"))) {
        y_axis_titles <- rep("RNA measured in nasal swab/wash (10^x copies/mL or Log10TCID50eq/mL)", length(list_comparisons))
      } else if ("ViralLoad" %in% list_comparisons) {
        y_axis_titles <- rep("RNA measured in inoculated animals in nasal swab/wash (10^x copies/mL or TCID50/mL)", length(list_comparisons))
      } else if ("type" %in% list_comparisons) {
        y_axis_titles <- rep("RNA measured in nasal/oral swab/wash (10^x copies/mL or TCID50/mL)", length(list_comparisons))
      }
    } else if (i %in% 4:6) {
      if (all(list_comparisons %in% c("Virus_cat", "route", "Species"))) {
        y_axis_titles <- rep("Infectious virus measured in nasal swab/wash (Log10TCID50/mL)", length(list_comparisons))
      } else if ("ViralLoad" %in% list_comparisons) {
        y_axis_titles <- rep("Infectious virus measured in inoculated animals in nasal swab/wash (Log10TCID50/mL)", length(list_comparisons))
      } else if ("type" %in% list_comparisons) {
        y_axis_titles <- rep("Infectious virus measured in nasal/oral swab/wash (Log10TCID50/mL)", length(list_comparisons))
      }
    }
  
  
  for(j in 1:length(list_comparisons)){
    plot_for_loop <- function(df,x_var,y_var,y_label){
      
      x_levels_list <- list(
        "route" = c("Inoculated","Direct contact", "Indirect contact"),
        "Virus_cat" = c("Index like", "Alpha", "Beta", "Delta", "Gamma", "Omicron"),
        "Species" = c("Ferret","Mink"),
        "ViralLoad" = c("<10^5", "10^5", ">10^5"),
        "type" = c("OralSwabs", "NasalSwabs", "NasalWash")
      )
      
      if (x_var %in% names(x_levels_list)) {
        df[[x_var]] <- factor(df[[x_var]], levels = x_levels_list[[x_var]])
      }
      
      ifelse(colnames(df[list_comparisons[j]])=='Virus_cat',color_choice <- Virus_cat_colors,
        ifelse(colnames(df[list_comparisons[j]])=='route',color_choice <- route_colors,
          ifelse(colnames(df[list_comparisons[j]])=="Species",color_choice <- Species_colors,
            ifelse(colnames(df[list_comparisons[j]])=="ViralLoad",color_choice <- ViralLoad_colors,color_choice <- type_colors))))
      
      ## Filter only rows where y_var is not NA 
      valid_rows <- !is.na(df[[y_var]]) 
      ## Calculate counts for there rows 
      counts <- table(df[[x_var]][valid_rows]) 
      ## Make sure that counts follow the factor levels
      x_levels <- levels(df[[x_var]]) 
      new_labels <- paste0(x_levels, " (n=", counts[x_levels], ")") 
      names(new_labels) <- x_levels
      
    ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]], fill= .data[[x_var]])) + 
      geom_boxplot() +
      
      scale_fill_manual(values = color_choice, drop = FALSE) +
      geom_jitter(aes(color = Species, shape = Species), size = 0.8, alpha = 0.9) +
      
      scale_shape_manual(values = c("Ferret" = 16, "Mink" = 17)) +  # 16 = cirkel, 17 = trekant
      scale_color_manual(values = c("Ferret" ="#4777EFFF", "Mink" = "#62FC6BFF")) +
      theme_bw() +
      theme(
        legend.position="none",
        plot.title = element_text(size=11)
      ) +
      ggtitle(" ") +
      xlab("") +
      ylab(y_label) +
      scale_y_continuous()+
      scale_x_discrete(labels = new_labels, guide = guide_axis(angle = 45)) 
  }
  
    plot_list <- colnames(df[list_measure]) %>% 
      set_names( ~ paste0(colnames(df[list_comparisons[j]]), "_",colnames(df[list_measure]))) %>%
      map( ~ plot_for_loop(df, colnames(df[list_comparisons[j]]),
                           colnames(df[list_measure]),
                           y_axis_titles[j]))
    
    all_plot_lists <- append(all_plot_lists, list(plot_list))
  }
  
  for(m in 1:length(list_comparisons)){
    groups <- unique(df[[colnames(df[list_comparisons[m]])]])
    
    tableM <- as.data.frame(matrix(ncol = 9, nrow = 0))
    
    for(k in 1:length(groups)){
      select  <- data.frame(comparison =  df[[colnames(df[list_comparisons[m]])]],
                            measure    =  as.numeric(df[[colnames(df[list_measure])]]))
      
      tableM <- rbind(tableM,
                      c(groups[k],
                        quantile(select$measure[which(select[,1]==groups[k])],c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=T),
                        sum(!is.na(select$measure[which(select[,1]==groups[k])]))))
      colnames(tableM) <- c("Category","0","0.05","0.25","0.5","0.75","0.95","1","n")
    }
    write_xlsx(tableM,
               paste(colnames(df[list_comparisons[m]]),
                     colnames(df[list_measure]),
                     '.xlsx',
                     sep=""))
    
  }
}

#####DESCR ANALYSES OF PROB OF TRANSMISSION####################################

trans <- readings %>%
  filter(`Experimental/Field`=="Experimentally") %>%
  select(Title,Virus_cat,ViralLoad,route,
         `TransmissionRoute (direct, neighbour, airborne)`,N_symptoms,`proportion_RNApos (field studies)`,
         `proportion_virus_pos (field studies)`,`duration of contact (days)`,Species)

names(trans) <- c('Title','Virus_cat','ViralLoad','route',
                  'TransmRoute','N_symptoms','propPCR',
                  'propVirus','durContact','Species')

trans$PCRpos   <- NA
trans$PCR_N    <- NA
trans$Viruspos <- NA
trans$Virus_N  <- NA

for(i in 1:length(trans$Title)){
  # Split strengen ved "/"
  parts1 <- strsplit(trans$propPCR[i], "/")[[1]]
  parts2 <- strsplit(trans$propVirus[i], "/")[[1]]
  
  # Konverter til tal  
  trans$PCRpos[i]   <- as.numeric(parts1[1])
  trans$PCR_N[i]    <- as.numeric(parts1[2])
  trans$Viruspos[i] <- as.numeric(parts2[1])
  trans$Virus_N[i]  <- as.numeric(parts2[2])
}

trans$pPCR   <- trans$PCRpos/trans$PCR_N
trans$pVirus <- trans$Viruspos/trans$Virus_N

list_comparisons  <- c(2,2,4,10)

y_axis_titles     <- c(pPCR="Proportion of animals infected, measured by RNA",
                       pVirus="Proportion of animals infected, measured by infectious virus")

for(i in 1:length(list_comparisons)){
  trans1  <-  if (i==2) {trans1 <- trans[which(trans$route%in%c("Direct contact","Indirect contact")),]} else {
    {trans1 <- trans}}
  
  # Full or reduced dataset
  dataset_type <- if (i == 2) "reduced" else "full"
  
  plot_for_loop <- function(dataset,x_var,y_var,y_label){
    
    x_levels_list <- list(
      "route" = c("Inoculated", "Direct contact", "Indirect contact"),
      "Virus_cat" = c("Index like", "Alpha", "Beta", "Delta", "Gamma", "Omicron"),
      "ViralLoad" = c("<10^5", "10^5", ">10^5"),
      "Species" = c("Ferret","Mink")
    )
    
    if (x_var %in% names(x_levels_list)) {
      dataset[[x_var]] <- factor(dataset[[x_var]], levels = x_levels_list[[x_var]])
    }
    
    ifelse(colnames(dataset[list_comparisons[i]])=='Virus_cat',color_choice <- Virus_cat_colors,
           ifelse(colnames(dataset[list_comparisons[i]])=='route',color_choice <- route_colors,
                  ifelse(colnames(dataset[list_comparisons[i]])=="Species",color_choice <- Species_colors,
                         ifelse(colnames(dataset[list_comparisons[i]])=="ViralLoad",color_choice <- ViralLoad_colors,color_choice <- type_colors))))
    
    ## Filter only rows where y_var is not NA 
    valid_rows <- !is.na(dataset[[y_var]]) 
    ## Calculate counts for there rows 
    counts <- table(dataset[[x_var]][valid_rows]) 
    ## Make sure that counts follow the factor levels
    x_levels <- levels(dataset[[x_var]]) 
    new_labels <- paste0(x_levels, " (n=", counts[x_levels], ")") 
    names(new_labels) <- x_levels
    
    ggplot(dataset, aes(x = .data[[x_var]], y = .data[[y_var]], fill= .data[[x_var]])) + 
      geom_boxplot() +
      
      scale_fill_manual(values = color_choice, drop = FALSE) +
      geom_jitter(aes(color = Species, shape = Species), size = 0.8, alpha = 0.9) +
      
      scale_shape_manual(values = c("Ferret" = 16, "Mink" = 17)) +  # 16 = cirkel, 17 = trekant
      scale_color_manual(values = c("Ferret" ="#4777EFFF", "Mink" = "#62FC6BFF")) +
      theme_bw() +
      theme(
        legend.position="none",
        plot.title = element_text(size=11)
      ) +
      xlab("") +
      ylab(y_label) +
      scale_y_continuous(limits=c(0, 1))+
      scale_x_discrete(labels = new_labels, guide = guide_axis(angle = 45)) 
  }
  
  plot_list <- colnames(trans1[,c(15,16)]) %>% 
    set_names( ~ paste0(colnames(trans1)[list_comparisons[i]], "_",.x)) %>%
    map( ~ plot_for_loop(trans1, colnames(trans1)[list_comparisons[i]], .x, y_axis_titles[.x]))
  
  all_plot_lists <- append(all_plot_lists, list(plot_list))
  
  ## Create excell
  
  for(m in 1:length(colnames(trans1[,c(15,16)]))){
    groups <- unique(trans1[[colnames(trans1)[list_comparisons[i]]]])
    
    tableM <- as.data.frame(matrix(ncol = 9, nrow = 0))
    
    for(k in 1:length(groups)){
      select  <- cbind(trans1[[colnames(trans1[,c(15,16)])[m]]],trans1[[colnames(trans1[list_comparisons[i]])]])
      
      tableM <- rbind(tableM,
                      c(groups[k],
                        quantile(as.numeric(select[which(select[,2]==groups[k]),1]),c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=T),
                        sum(!is.na(select[which(select[,2]==groups[k]),1]))))
      colnames(tableM) <- c("Category","0","0.05","0.25","0.5","0.75","0.95","1","n")
    }
    write_xlsx(tableM,
               paste(dataset_type,
                     colnames(trans1[,c(15,16)])[m],
                     colnames(trans1)[list_comparisons[i]],
                     '.xlsx',
                     sep=""))
    
  }
  
}



#####POINT PLOTS PCR VERSUS INF VIRUS; ONLY FROM PUB.WITH BOTH MEASURES########

listDF2        <- c('PCR_nasal_oral')
x_names2       <- c('logOralPCR','logOralVirus')
y_names2       <- c('logNasalPCR', 'logNasalVirus')


## Point plot af oral/nasal i PCR i PCR_nasal_oral
## There are 0 observations of infectious virus measured in oral AND nasal

## Looking at observations with oral AND nasal meassures in PCR


for(i in 1:1){
  df <- get(listDF2[i])

  ggplot(df) +
  geom_point(aes(x = .data[[colnames(df[x_names2[i]])]], y = .data[[colnames(df[y_names2[i]])]],
                 color = Nasal_excr_measured_in))
   
  glm_n_o <- glm(df[[colnames(df[y_names2[i]])]]~df[[colnames(df[x_names2[i]])]]+df$Nasal_excr_measured_in)
  summary(glm_n_o)
  
  pointPlot <- ggplot(df, aes(x = .data[[colnames(df[x_names2[i]])]], y = .data[[colnames(df[y_names2[i]])]],
                              color = Nasal_excr_measured_in)) +
                      geom_point(size = 2) +
                      geom_smooth(method = "lm", se = FALSE) +  # En linje pr. gruppe
                      scale_color_manual(values = c("Wash" = "#FFB6C1", "Swabs" = "#B2DFEE")) +
                      geom_smooth(aes(color = NULL), method = "lm", se = FALSE, 
                                  color = "grey40") +  # Samlet linje

                          labs(title = "",
                           x = "Oral excretion measured by RNA (Log10 copies or TCID50/mL)",
                           y = "Nasal excretion measured by RNA (Log10 copies or TCID50/mL)")
  
 ggsave(paste(listDF2[i],
                 '.tiff',
                 sep="_"), 
        pointPlot,
        device = tiff,
        width = 10, height = 8, units = "in", dpi = 300)
 
 ggsave(paste(listDF2[i],
              '.pdf',
              sep="_"), 
        pointPlot,
        device = pdf,
        width = 10, height = 8, units = "in", dpi = 300)
 
 ggsave(paste(listDF2[i],
              '.eps',
              sep="_"), 
        pointPlot,
        device = cairo_ps,
        width = 10, height = 8, units = "in", dpi = 300)
}



#####ALL TABLE AS OVERALL (NO GROUPS)###########################################

tableM <- as.data.frame(matrix(ncol = 9, nrow = 0))
colnames(tableM) <- c("measure","0","0.05","0.25","0.5","0.75","0.95","1","n")  

liste <- c("dpi_firstSymptoms","durSympt",
           "firstRNA","peakRNA","lastRNA",
           "firstVirus","peakVirus","lastVirus")

for(i in 1:8){
  tableM <- rbind(tableM,
                    c(liste[i],
                      quantile(sympt2[,(i+4)],
                               c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=T),
                      sum(!is.na(sympt2[,(i+4)]))))
}
colnames(tableM) <- c("measure","0","0.05","0.25","0.5","0.75","0.95","1","n")  

tableM <- rbind(tableM,
                  c("logPCR",
                    quantile(c(logPCR$logNasalPCR,logPCR$logOralPCR),
                               c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=T),
                    sum(!is.na(c(logPCR$logNasalPCR,logPCR$logOralPCR)))))

colnames(tableM) <- c("measure","0","0.05","0.25","0.5","0.75","0.95","1","n")  

tableM <- rbind(tableM,
                c("logVirus",
                  quantile(c(logVirus$logNasalVirus,logVirus$logOralVirus),
                           c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=T),
                  sum(!is.na(c(logVirus$logNasalVirus,logVirus$logOralVirus)))))

colnames(tableM) <- c("measure","0","0.05","0.25","0.5","0.75","0.95","1","n")  


liste <- c("pPCR","pVirus")

for(i in 1:2){
  tableM <- rbind(tableM,
                  c(liste[i],
                    quantile(trans[,(i+14)],
                             c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=T),
                    sum(!is.na(trans[,(i+14)]))))
}
colnames(tableM) <- c("measure","0","0.05","0.25","0.5","0.75","0.95","1","n")  

write_xlsx(tableM,'GeneralMeasures.xlsx')



#####TRANSMISSION RATES ESTIMATION##############################################

library(tidyverse)
library(readxl)
library(readr)

# Notes:
# New infections (I_new) = dI/dt = (beta * S * I) / N * t
# Thus:
# beta = I_new / (S*I) / N * t

##### Vectorised beta estimation function: #####################################

beta_estimate <- function(S, I, I_new, t, mode = c("density", "frequency","sqrt","density.50","density.150")) {
  mode <- match.arg(mode)
  
  # Compatible lengths
  lens <- c(length(S), length(I), length(I_new), length(t))
  L <- max(lens)
  if (!all(lens == 1 | lens == L)) stop("All inputs must be same length or length 1.")
  
  S <- rep_len(S, L); I <- rep_len(I, L); I_new <- rep_len(I_new, L); t <- rep_len(t, L)
  
  # Force numeric
  to_num <- function(x) suppressWarnings(as.numeric(x))
  S <- to_num(S); I <- to_num(I); I_new <- to_num(I_new); t <- to_num(t)
  
  N <- S + I
  
  # Nonvalids -> NA
  valid_base <- !is.na(S) & !is.na(I) & !is.na(I_new) & !is.na(t) &
    S >= 0 & I >= 0 & I_new >= 0 & t > 0 & N > 0
  
  exposed <- rep(NA_real_, L)
  exposed[valid_base] <- if (mode == "density") {
    (S[valid_base] * I[valid_base]) * t[valid_base]
  } else  if (mode == "density.50") {
    (S[valid_base] * I[valid_base]) * (t[valid_base] * 0.5) # 50% exposure
  } else  if (mode == "density.150") {
    (S[valid_base] * I[valid_base]) * (t[valid_base] * 1.5) # 150% exposurre
  } else if (mode == "frequency") {
    (S[valid_base] * I[valid_base] / N[valid_base]) * t[valid_base]
  } else if (mode == "sqrt") {
    (S[valid_base] * I[valid_base] / sqrt(N[valid_base])) * t[valid_base]
  }
  
  out <- rep(NA_real_, L)
  ok  <- valid_base & is.finite(exposed) & exposed > 0
  out[ok] <- I_new[ok] / exposed[ok]
  
  
  out[!is.finite(out)] <- NA_real_
  
  out
}

##### Get data #################################################################

dat <- readings %>%
  filter(Type == "Transmission") %>%
  select(Title,Virus_cat,ViralLoad,route,Group,N_animals,
         `TransmissionRoute (direct, neighbour, airborne)`,N_symptoms,
         dpi_RNAdetection, dpi_infectious_virus,
         `proportion_RNApos (field studies)`,`proportion_virus_pos (field studies)`,
         `duration of contact (days)`,Species)

names(dat) <- c('Title','Virus_cat','ViralLoad','route','Group','N_animals',
                  'TransmRoute','N_symptoms',
                'dpi_RNAdetection','dpi_infectious_virus',
                'propPCR','propVirus',
                'durContact','Species')

dat$PCRpos   <- NA
dat$PCR_N    <- NA
dat$Viruspos <- NA
dat$Virus_N  <- NA

for(i in 1:length(dat$Title)){
  # Split strengen ved "/"
  parts1 <- strsplit(dat$propPCR[i], "/")[[1]]
  parts2 <- strsplit(dat$propVirus[i], "/")[[1]]
  
  # Konverter til tal  
  dat$PCRpos[i]   <- as.numeric(parts1[1])
  dat$PCR_N[i]    <- as.numeric(parts1[2])
  dat$Viruspos[i] <- as.numeric(parts2[1])
  dat$Virus_N[i]  <- as.numeric(parts2[2])
}

dat$pPCR   <- dat$PCRpos/dat$PCR_N
dat$pVirus <- dat$Viruspos/dat$Virus_N

dat2 <- dat 

dat3 <- dat2 |>
  # Gør felter numeriske
  mutate(across(
    c(durContact, dpi_RNAdetection, dpi_infectious_virus),
    ~ suppressWarnings(parse_number(as.character(.)))
  )) |>
  # Beregn eksponeringer
  mutate(
    expTime_rna = if_else(
      !is.na(dpi_RNAdetection) & dpi_RNAdetection < durContact,
      dpi_RNAdetection,
      durContact
    ),
    expTime_virus = if_else(
      !is.na(dpi_infectious_virus) & dpi_infectious_virus < durContact,
      dpi_infectious_virus,
      durContact
    )
  )

## Creating incoluated and exposed by PCR and inf Virus

dat3 <- dat3 |> 
  group_by(Group) %>%
  mutate(
    Infected_PCR = dplyr::first(
      PCRpos[`route` %in% "Inoculated"],
      default = NA
    ),
    Exposed_PCR = dplyr::first(
      PCR_N[`route` %in% c(DirContact,IndContact)],
      default = NA
    ),
    Infected_Virus = dplyr::first(
      Viruspos[`route` %in% "Inoculated"],
      default = NA
    ),
    Exposed_Virus = dplyr::first(
      Virus_N[`route` %in% c(DirContact,IndContact)],
      default = NA
    )
  ) |> 
  #dplyr::filter(Type == "Transmission") |> 
  dplyr::filter(!is.na(Group)) |>                       # Her tager jeg alle med Group = NA ud. 20251113
  ungroup()


for(i in unique(dat3$Group)){
  if(sum(dat3$Group==i)>2){
    
    tmp1 <- dat3$PCRpos[which(dat3$Group==i&
                                       dat3$TransmRoute%in%c("Direct contact","1st contact"))]
    
    tmp2 <- dat3$Viruspos[which(dat3$Group==i&
                                dat3$TransmRoute%in%c("Direct contact","1st contact"))]
    
    dat3$Infected_PCR[which(dat3$Group==i&
                              dat3$TransmRoute=="2nd contact")] <- tmp1
    
    dat3$Infected_Virus[which(dat3$Group==i&
                              dat3$TransmRoute=="2nd contact")] <- tmp2
  }
}

## Remove observations for inoculated
dat4 <- dat3 |> filter(route!='Inoculated')
table(dat4$route)
table(dat4$Group)

##### Daily transmission rates#####################################################
## Define combinations as tibble

combos <- tibble::tibble(
  S      = list(dat4$Exposed_PCR, dat4$Exposed_Virus,
                dat4$Exposed_PCR, dat4$Exposed_Virus,
                dat4$Exposed_PCR, dat4$Exposed_Virus,
                dat4$Exposed_PCR, dat4$Exposed_Virus),
  I      = list(dat4$Infected_PCR, dat4$Infected_Virus,
                dat4$Infected_PCR, dat4$Infected_Virus,
                dat4$Infected_PCR, dat4$Infected_Virus,
                dat4$Infected_PCR, dat4$Infected_Virus),
  I_new  = list(dat4$PCRpos, dat4$Viruspos,
                dat4$PCRpos, dat4$Viruspos,
                dat4$PCRpos, dat4$Viruspos,
                dat4$PCRpos, dat4$Viruspos),
  t      = list(dat4$expTime_rna, dat4$expTime_virus,
                0.5, 0.5,
                1, 1,
                2, 2),
  mode   = c("density", "density",
             "density", "density",
             "density", "density",
             "density", "density"),
  name   = c("beta_est_PCR_dens", "beta_est_vir_dens",
             "beta_est_PCR_dens_05", "beta_est_vir_dens_05",
             "beta_est_PCR_dens_1", "beta_est_vir_dens_1",
             "beta_est_PCR_dens_2", "beta_est_vir_dens_2")
)

# Brug purrr::pwalk til at iterere og tilføje kolonner
pwalk(combos, function(S, I, I_new, t, mode, name) {
  dat4[[name]] <<- beta_estimate(S = S, I = I, I_new = I_new,
                                 t = t, mode = mode)
})

tableBeta <- as.data.frame(matrix(ncol = 10, nrow = 0))

for(i in 1:8){
  tableBeta <- rbind(tableBeta,
                     c(combos$name[i],"Direct contact",
                       quantile(dat4[[combos$name[i]]][dat4$route == "Direct contact" ],c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=T),
                       sum(!is.na(dat4[[combos$name[i]]][dat4$route == "Direct contact" ]))),
                     c(combos$name[i],"Indirect contact",
                       quantile(dat4[[combos$name[i]]][dat4$route == "Indirect contact" ],c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=T),
                       sum(!is.na(dat4[[combos$name[i]]][dat4$route == "Indirect contact" ]))))
}
colnames(tableBeta) <- c("Beta","Route","0","0.05","0.25","0.5","0.75","0.95","1","n")

write_xlsx(tableBeta, path = "transmission_rates.xlsx")

#####BETA Boxplots by route#############################################
Virus_cat_colors <- c(
  "Index like" = "#0D0887FF",
  "Alpha"      = "#6A00A8FF",
  "Beta"       = "#B12A90FF",
  "Delta"      = "#E16462FF",
  "Gamma"      = "#FCA636FF",
  "Omicron"    = "#F0F921FF"
)

route_colors <- c(
  "Inoculated"        = "#46337EFF",
  "Direct contact"    = "#31688EFF",
  "Indirect contact"  = "#35B779FF" 
)

densities <- c("dens", "dens_05", "dens_1", "dens_2")

plots_PCR <- map(seq_along(densities), function(i) {
  d <- densities[i]
    
  ## Calculate counts 
  counts <- table(dat4$route[which(!is.na(dat4$beta_est_PCR_dens))]) 
  x_levels <- levels(as.factor(dat4$route)) 
  new_labels <- paste0(x_levels, " (n=", counts[x_levels], ")") 
  names(new_labels) <- x_levels
  
    ggplot(dat4, aes(x = route, 
                   y = .data[[paste0("beta_est_PCR_", d)]], 
                   fill = route)) +
    geom_boxplot() +
    geom_jitter(aes(color = Species, shape = Species), size = 0.8, alpha = 0.9) +
    scale_fill_manual(values = route_colors[2:3], drop = FALSE) +
    scale_shape_manual(values = c("Ferret" = 16, "Mink" = 17)) +
    scale_color_manual(values = c("Ferret" ="#4777EFFF", "Mink" = "#62FC6BFF")) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = 11)) +
    xlab("") +
    ylab("Transmission rate, infection detected by viral RNA") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_discrete(labels = new_labels, guide = guide_axis(angle = 45)) 
})

names(plots_PCR) <- paste0("route_PCR_beta_",densities)
  
all_plot_lists <- append(all_plot_lists, list(plots_PCR))

plots_vir <- map(seq_along(densities), function(i) {
  d <- densities[i]
  
  ## Calculate counts 
  counts <- table(dat4$route[which(!is.na(dat4$beta_est_vir_dens))]) 
  x_levels <- levels(as.factor(dat4$route[which(!is.na(dat4$beta_est_vir_dens))]))
  new_labels <- paste0(x_levels, " (n=", counts[x_levels], ")") 
  names(new_labels) <- x_levels
  
  ggplot(dat4, aes(x = route, 
                   y = .data[[paste0("beta_est_vir_", d)]], 
                   fill = route)) +
    geom_boxplot() +
    geom_jitter(aes(color = Species, shape = Species), size = 0.8, alpha = 0.9) +
    scale_fill_manual(values = route_colors[2:3], drop = FALSE) +
    scale_shape_manual(values = c("Ferret" = 16, "Mink" = 17)) +
    scale_color_manual(values = c("Ferret" ="#4777EFFF", "Mink" = "#62FC6BFF")) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = 11)) +
    xlab("") +
    ylab("Transmission rate, infection detected by infectious virus") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_discrete(labels = new_labels, guide = guide_axis(angle = 45))
})

names(plots_vir) <- paste0("route_vir_beta_",densities)

all_plot_lists <- append(all_plot_lists, list(plots_vir))

#####BETA Boxplots by Virus_cat########################################################
dat4 <- dat4 %>%
  mutate(Virus_cat = factor(Virus_cat,
                            levels = c("Index like", "Alpha", "Beta", "Delta", "Gamma", "Omicron")))

plots_PCR_vircat <- map(seq_along(densities), function(i) {
  d <- densities[i]
  
  ## Calculate counts 
  counts <- table(dat4$Virus_cat[which(!is.na(dat4$beta_est_PCR_dens))]) 
  x_levels <- levels(as.factor(dat4$Virus_cat[which(!is.na(dat4$beta_est_PCR_dens))]))
  new_labels <- paste0(x_levels, " (n=", counts[x_levels], ")") 
  names(new_labels) <- x_levels
  
  ggplot(dat4, aes(x = Virus_cat, 
                   y = .data[[paste0("beta_est_PCR_", d)]], 
                   fill = Virus_cat)) +
    geom_boxplot() +
    geom_jitter(aes(color = Species, shape = Species), size = 0.8, alpha = 0.9) +
    scale_fill_manual(values = Virus_cat_colors, drop = FALSE) +
    scale_shape_manual(values = c("Ferret" = 16, "Mink" = 17)) +
    scale_color_manual(values = c("Ferret" ="#4777EFFF", "Mink" = "#62FC6BFF")) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = 11)) +
    xlab("") +
    ylab("Transmission rate, infection detected by viral RNA") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_discrete(labels = new_labels, guide = guide_axis(angle = 45)) 
})

names(plots_PCR_vircat) <- paste0("VirusCat_PCR_beta_",densities)

all_plot_lists <- append(all_plot_lists, list(plots_PCR_vircat))

plots_vir_vircat <- map(seq_along(densities), function(i) {
  d <- densities[i]
  
  ## Calculate counts 
  counts <- table(dat4$Virus_cat[which(!is.na(dat4$beta_est_vir_dens))]) 
  x_levels <- levels(as.factor(dat4$Virus_cat[which(!is.na(dat4$beta_est_vir_dens))]))
  new_labels <- paste0(x_levels, " (n=", counts[x_levels], ")") 
  names(new_labels) <- x_levels
  
  ggplot(dat4, aes(x = Virus_cat, 
                   y = .data[[paste0("beta_est_vir_", d)]], 
                   fill = Virus_cat)) +
    geom_boxplot() +
    geom_jitter(aes(color = Species, shape = Species), size = 0.8, alpha = 0.9) +
    scale_fill_manual(values = Virus_cat_colors, drop = FALSE) +
    scale_shape_manual(values = c("Ferret" = 16, "Mink" = 17)) +
    scale_color_manual(values = c("Ferret" ="#4777EFFF", "Mink" = "#62FC6BFF")) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = 11)) +
    xlab("") +
    ylab("Transmission rate, infection detected by infectious virus") +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_discrete(labels = new_labels, guide = guide_axis(angle = 45))
})

names(plots_vir_vircat) <- paste0("VirusCat_vir_beta_",densities)

all_plot_lists <- append(all_plot_lists, list(plots_vir_vircat))

#####PLOTS######################################################################
for (i in seq_along(all_plot_lists)) { 
  plot_list <- all_plot_lists[[i]] 
  
  for (name in names(plot_list)) { 
    ggsave( filename = paste0(name, ".tiff"), 
            plot = plot_list[[name]], 
            device = tiff,
            width = 10, height = 8, dpi = 300) 
    ggsave( filename = paste0(name, ".pdf"), 
            plot = plot_list[[name]], 
            device = pdf,
            width = 10, height = 8, dpi = 300) 
    ggsave( filename = paste0(name, ".eps"), 
            plot = plot_list[[name]], 
            device = cairo_ps,
            width = 10, height = 8, dpi = 300) 
  } 
  }
    
library(cowplot)

#####FIGURE 2###################################################################
Fig2 <- plot_grid(plotlist = list(all_plot_lists[[2]][["route_firstRNA"]],
                                  all_plot_lists[[2]][["route_peakRNA"]],
                                  all_plot_lists[[2]][["route_lastRNA"]],
                                  all_plot_lists[[2]][["route_firstVirus"]],
                                  all_plot_lists[[2]][["route_peakVirus"]],
                                  all_plot_lists[[2]][["route_lastVirus"]]),
                  labels = "AUTO", # produces a, b, c, d, e, f 
                  label_size = 12, label_fontface = "bold",
                  ncol = 6)


p_values <- format(c(0.002,0.0034,0.0024,0.0018,0.0004,0.0001,0.0315), scientific=F)

Fig2_with_p <- ggdraw(Fig2) +
  draw_label(paste0("p = ", p_values[1]), x = 0.26, y = 0.9, size = 8, angle = 90) +
  draw_label(paste0("p = ", p_values[2]), x = 0.304, y = 0.9, size = 8, angle = 90) +
  draw_label(paste0("p = ", p_values[3]), x = 0.592, y = 0.9, size = 8, angle = 90) +
  draw_label(paste0("p = ", p_values[4]), x = 0.635, y = 0.9, size = 8, angle = 90) +
  draw_label(paste0("p = ", p_values[5]), x = 0.76, y = 0.9, size = 8, angle = 90) +
  draw_label(paste0("p = ", p_values[7]), x = 0.785, y = 0.8, size = 8) +
  draw_label(paste0("p < ", p_values[6]), x = 0.805, y = 0.9, size = 8, angle = 90)
  
ggsave(paste('route_','Fig2','.tiff',sep=""), 
       Fig2_with_p,
       device = tiff,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('route_','Fig2','.pdf',sep=""), 
       Fig2_with_p,
       device = pdf,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('route_','Fig2','.eps',sep=""), 
       Fig2_with_p,
       device = cairo_ps,
       width = 20, height = 8, units = "in", dpi = 300)

#####FIGURE 3###################################################################
Fig3 <- plot_grid(plotlist = list(all_plot_lists[[6]][["route_logNasalPCR"]],
                                  all_plot_lists[[11]][["route_logNasalVirus"]],
                                  all_plot_lists[[17]][["route_pPCR"]],
                                  all_plot_lists[[17]][["route_pVirus"]]),
                  labels = "AUTO", # produces a, b, c, d, e, f 
                  label_size = 12, label_fontface = "bold",
                  ncol = 2)

ggsave(paste('route_','Fig3','.tiff',sep=""), 
       Fig3,
       device = tiff,
       width = 10, height = 15, units = "in", dpi = 300)
ggsave(paste('route_','Fig3','.pdf',sep=""), 
       Fig3,
       device = pdf,
       width = 10, height = 15, units = "in", dpi = 300)
ggsave(paste('route_','Fig3','.eps',sep=""), 
       Fig3,
       device = cairo_ps,
       width = 10, height = 15, units = "in", dpi = 300)

#####FIGURE 4###################################################################
Fig4 <- plot_grid(plotlist = list(all_plot_lists[[19]][["route_PCR_beta_dens"]],
                                  all_plot_lists[[19]][["route_PCR_beta_dens_05"]],
                                  all_plot_lists[[19]][["route_PCR_beta_dens_1"]],
                                  all_plot_lists[[19]][["route_PCR_beta_dens_2"]],
                                  all_plot_lists[[20]][["route_vir_beta_dens"]],
                                  all_plot_lists[[20]][["route_vir_beta_dens_05"]],
                                  all_plot_lists[[20]][["route_vir_beta_dens_1"]],
                                  all_plot_lists[[20]][["route_vir_beta_dens_2"]]),
                  labels = "AUTO", # produces a, b, c, d, e, f 
                  label_size = 12, label_fontface = "bold",
                  ncol = 8)


ggsave(paste('route_','Fig4','.tiff',sep=""), 
       Fig4,
       device = tiff,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('route_','Fig4','.pdf',sep=""), 
       Fig4,
       device = pdf,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('route_','Fig4','.eps',sep=""), 
       Fig4,
       device = cairo_ps,
       width = 20, height = 8, units = "in", dpi = 300)

#####SUPPL FIGURE D############################################################
## D1
SupplFigD1 <- plot_grid(plotlist = list(all_plot_lists[[1]][["Virus_cat_firstRNA"]],
                                  all_plot_lists[[1]][["Virus_cat_peakRNA"]],
                                  all_plot_lists[[1]][["Virus_cat_lastRNA"]],
                                  all_plot_lists[[1]][["Virus_cat_firstVirus"]],
                                  all_plot_lists[[1]][["Virus_cat_peakVirus"]],
                                  all_plot_lists[[1]][["Virus_cat_lastVirus"]]),
                        labels = "AUTO", # produces a, b, c, d, e, f 
                        label_size = 12, label_fontface = "bold",
                        ncol = 6)


p_values <- format(c(0.032,0.027,0.0036), scientific=F)

SupplFigD1_with_p <- ggdraw(SupplFigD1) +
  draw_label(paste0("p = ", p_values[1]), x = 0.438, y = 0.45, size = 8, angle = 90) +
  draw_label(paste0("p = ", p_values[2]), x = 0.915, y = 0.9, size = 8, angle = 90) +
  draw_label(paste0("p = ", p_values[3]), x = 0.985, y = 0.9, size = 8, angle = 90)

ggsave(paste('Virus_cat_','SupplFigD1_with_p','.tiff',sep=""), 
       SupplFigD1_with_p,
       device = tiff,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Virus_cat_','SupplFigD1_with_p','.pdf',sep=""), 
       SupplFigD1_with_p,
       device = pdf,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Virus_cat_','SupplFigD1_with_p','.eps',sep=""), 
       SupplFigD1_with_p,
       device = cairo_ps,
       width = 20, height = 8, units = "in", dpi = 300)

## D2
SupplFigD2 <- plot_grid(plotlist = list(all_plot_lists[[5]][["Virus_cat_logNasalPCR"]],
                                        all_plot_lists[[10]][["Virus_cat_logNasalVirus"]]),
                        labels = "AUTO", # produces a, b, c, d, e, f 
                        label_size = 12, label_fontface = "bold",
                        ncol = 2)

ggsave(paste('Virus_cat_','SupplFigD2','.tiff',sep=""), 
       SupplFigD2,
       device = tiff,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Virus_cat_','SupplFigD2','.pdf',sep=""), 
       SupplFigD2,
       device = pdf,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Virus_cat_','SupplFigD2','.eps',sep=""), 
       SupplFigD2,
       device = cairo_ps,
       width = 20, height = 8, units = "in", dpi = 300)

## D3
SupplFigD3 <- plot_grid(plotlist = list(all_plot_lists[[15]][["Virus_cat_pPCR"]],
                                        all_plot_lists[[15]][["Virus_cat_pVirus"]]),
                        labels = "AUTO", # produces a, b, c, d, e, f 
                        label_size = 12, label_fontface = "bold",
                        ncol = 2)

ggsave(paste('Virus_cat_','SupplFigD3','.tiff',sep=""), 
       SupplFigD3,
       device = tiff,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Virus_cat_','SupplFigD3','.pdf',sep=""), 
       SupplFigD3,
       device = pdf,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Virus_cat_','SupplFigD3','.eps',sep=""), 
       SupplFigD3,
       device = cairo_ps,
       width = 20, height = 8, units = "in", dpi = 300)

## D4
SupplFigD4 <- plot_grid(plotlist = list(all_plot_lists[[16]][["Virus_cat_pPCR"]],
                                        all_plot_lists[[16]][["Virus_cat_pVirus"]]),
                        labels = "AUTO", # produces a, b, c, d, e, f 
                        label_size = 12, label_fontface = "bold",
                        ncol = 2)

ggsave(paste('Virus_cat_','SupplFigD4','.tiff',sep=""), 
       SupplFigD4,
       device = tiff,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Virus_cat_','SupplFigD4','.pdf',sep=""), 
       SupplFigD4,
       device = pdf,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Virus_cat_','SupplFigD4','.eps',sep=""), 
       SupplFigD4,
       device = cairo_ps,
       width = 20, height = 8, units = "in", dpi = 300)

## D5
SupplFigD5 <- plot_grid(plotlist = list(all_plot_lists[[21]][["VirusCat_PCR_beta_dens"]],
                                        all_plot_lists[[21]][["VirusCat_PCR_beta_dens_05"]],
                                        all_plot_lists[[21]][["VirusCat_PCR_beta_dens_1"]],
                                        all_plot_lists[[21]][["VirusCat_PCR_beta_dens_2"]],
                                        all_plot_lists[[22]][["VirusCat_vir_beta_dens"]],
                                        all_plot_lists[[22]][["VirusCat_vir_beta_dens_05"]],
                                        all_plot_lists[[22]][["VirusCat_vir_beta_dens_1"]],
                                        all_plot_lists[[22]][["VirusCat_vir_beta_dens_2"]]),
                        labels = "AUTO", # produces a, b, c, d, e, f 
                        label_size = 12, label_fontface = "bold",
                        ncol = 8)

ggsave(paste('Virus_cat_','SupplFigD5','.tiff',sep=""), 
       SupplFigD5,
       device = tiff,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Virus_cat_','SupplFigD5','.pdf',sep=""), 
       SupplFigD5,
       device = pdf,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Virus_cat_','SupplFigD5','.eps',sep=""), 
       SupplFigD5,
       device = cairo_ps,
       width = 20, height = 8, units = "in", dpi = 300)

#####SUPPL FIGURE E############################################################
## E1
SupplFigE1 <- plot_grid(plotlist = list(all_plot_lists[[3]][["Species_firstRNA"]],
                                        all_plot_lists[[3]][["Species_peakRNA"]],
                                        all_plot_lists[[3]][["Species_lastRNA"]],
                                        all_plot_lists[[3]][["Species_firstVirus"]],
                                        all_plot_lists[[3]][["Species_peakVirus"]],
                                        all_plot_lists[[3]][["Species_lastVirus"]]),
                        labels = "AUTO", # produces a, b, c, d, e, f 
                        label_size = 12, label_fontface = "bold",
                        ncol = 6)

ggsave(paste('Species_','SupplFigE1','.tiff',sep=""), 
       SupplFigE1,
       device = tiff,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Species_','SupplFigE1','.pdf',sep=""), 
       SupplFigE1,
       device = pdf,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Species_','SupplFigE1','.eps',sep=""), 
       SupplFigE1,
       device = cairo_ps,
       width = 20, height = 8, units = "in", dpi = 300)

## E2
SupplFigE2 <- plot_grid(plotlist = list(all_plot_lists[[7]][["Species_logNasalPCR"]],
                                        all_plot_lists[[12]][["Species_logNasalVirus"]]),
                        labels = "AUTO", # produces a, b, c, d, e, f 
                        label_size = 12, label_fontface = "bold",
                        ncol = 2)

ggsave(paste('Species_','SupplFigE2','.tiff',sep=""), 
       SupplFigE2,
       device = tiff,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Species_','SupplFigE2','.pdf',sep=""), 
       SupplFigE2,
       device = pdf,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Species_','SupplFigE2','.eps',sep=""), 
       SupplFigE2,
       device = cairo_ps,
       width = 20, height = 8, units = "in", dpi = 300)

## E3
SupplFigE3 <- plot_grid(plotlist = list(all_plot_lists[[18]][["Species_pPCR"]],
                                        all_plot_lists[[18]][["Species_pVirus"]]),
                        labels = "AUTO", # produces a, b, c, d, e, f 
                        label_size = 12, label_fontface = "bold",
                        ncol = 2)

ggsave(paste('Species_','SupplFigE3','.tiff',sep=""), 
       SupplFigE3,
       device = tiff,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Species_','SupplFigE3','.pdf',sep=""), 
       SupplFigE3,
       device = pdf,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Species_','SupplFigE3','.eps',sep=""), 
       SupplFigE3,
       device = cairo_ps,
       width = 20, height = 8, units = "in", dpi = 300)

## E4
SupplFigE4 <- plot_grid(plotlist = list(all_plot_lists[[3]][["Species_dpi_firstSymptoms"]],
                                        all_plot_lists[[3]][["Species_durSympt"]],
                                        all_plot_lists[[4]][["ViralLoad_dpi_firstSymptoms"]],
                                        all_plot_lists[[4]][["ViralLoad_durSympt"]]),
                        labels = "AUTO", # produces a, b, c, d, e, f 
                        label_size = 12, label_fontface = "bold",
                        ncol = 4)

ggsave(paste('Species_','SupplFigE4','.tiff',sep=""), 
       SupplFigE4,
       device = tiff,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Species_','SupplFigE4','.pdf',sep=""), 
       SupplFigE4,
       device = pdf,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Species_','SupplFigE4','.eps',sep=""), 
       SupplFigE4,
       device = cairo_ps,
       width = 20, height = 8, units = "in", dpi = 300)

#####SUPPL FIGURE F############################################################
## F1
SupplFigF1 <- plot_grid(plotlist = list(all_plot_lists[[4]][["ViralLoad_firstRNA"]],
                                        all_plot_lists[[4]][["ViralLoad_peakRNA"]],
                                        all_plot_lists[[4]][["ViralLoad_lastRNA"]],
                                        all_plot_lists[[4]][["ViralLoad_firstVirus"]],
                                        all_plot_lists[[4]][["ViralLoad_peakVirus"]],
                                        all_plot_lists[[4]][["ViralLoad_lastVirus"]]),
                        labels = "AUTO", # produces a, b, c, d, e, f 
                        label_size = 12, label_fontface = "bold",
                        ncol = 6)

ggsave(paste('ViralLoad_','SupplFigF1','.tiff',sep=""), 
       SupplFigF1,
       device = tiff,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('ViralLoad_','SupplFigF1','.pdf',sep=""), 
       SupplFigF1,
       device = pdf,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('ViralLoad_','SupplFigF1','.eps',sep=""), 
       SupplFigF1,
       device = cairo_ps,
       width = 20, height = 8, units = "in", dpi = 300)

## F2
SupplFigF2 <- plot_grid(plotlist = list(all_plot_lists[[8]][["ViralLoad_logNasalPCR"]],
                                        all_plot_lists[[13]][["ViralLoad_logNasalVirus"]]),
                        labels = "AUTO", # produces a, b, c, d, e, f 
                        label_size = 12, label_fontface = "bold",
                        ncol = 2)

ggsave(paste('ViralLoad_','SupplFigF2','.tiff',sep=""), 
       SupplFigF2,
       device = tiff,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('ViralLoad_','SupplFigF2','.pdf',sep=""), 
       SupplFigF2,
       device = pdf,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('ViralLoad_','SupplFigF2','.eps',sep=""), 
       SupplFigF2,
       device = cairo_ps,
       width = 20, height = 8, units = "in", dpi = 300)

#####SUPPL FIGURE G############################################################
SupplFigG1 <- plot_grid(plotlist = list(all_plot_lists[[9]][["type_logPCR"]],
                                        all_plot_lists[[14]][["type_logVirus"]]),
                        labels = "AUTO", # produces a, b, c, d, e, f 
                        label_size = 12, label_fontface = "bold",
                        ncol = 2)

ggsave(paste('Type_','SupplFigG1','.tiff',sep=""), 
       SupplFigG1,
       device = tiff,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Type_','SupplFigG1','.pdf',sep=""), 
       SupplFigG1,
       device = pdf,
       width = 20, height = 8, units = "in", dpi = 300)
ggsave(paste('Type_','SupplFigG1','.eps',sep=""), 
       SupplFigG1,
       device = cairo_ps,
       width = 20, height = 8, units = "in", dpi = 300)
#####COUNTING OBSERVATIONS######################################################

## Counting numbers of observation measuring viral RNA
sum(!is.na(readings$`Oral/Throat,MaxVirusExcreted_PCR_original`)|!is.na(readings$`Nasal,MaxVirusExcreted_PCR_original`),na.rm=T)
sum(readings$`Oral/Throat,MaxVirusExcreted_PCR_original`=="ND"|readings$`Nasal,MaxVirusExcreted_PCR_original`=="ND",na.rm=T)

sum(!is.na(logPCR$logNasalPCR)|!is.na(logPCR$logOralPCR))
sum(!is.na(logPCR$logNasalPCR))
sum(!is.na(logPCR$logOralPCR))
sum(!is.na(logPCR$logNasalPCR)&!is.na(logPCR$logOralPCR))

## Counting numbers of observation measuring infectious virus
sum(!is.na(readings$`Oral/Throat, MaxVirusExcreted_INFvirus_original`)|!is.na(readings$`Nasal, MaxVirusExcreted_INFvirus_original`),na.rm=T)
sum(readings$`Oral/Throat, MaxVirusExcreted_INFvirus_original`=="ND"|readings$`Nasal, MaxVirusExcreted_INFvirus_original`=="ND",na.rm=T)

sum(!is.na(logVirus$logNasalVirus)|!is.na(logVirus$logOralVirus))
sum(!is.na(logVirus$logNasalVirus))
sum(!is.na(logVirus$logOralVirus))
sum(!is.na(logVirus$logNasalVirus)&!is.na(logVirus$logOralVirus))

## Counting prob of inf
sum(!is.na(trans$pPCR))
sum(!is.na(trans$pVirus))

## Counting symptoms
sum(!is.na(readings$dpi_firstSymptoms)) ## 24
sum(is.na(readings$dpi_firstSymptoms)) ## 
sum(readings$`Duration of symptoms`=="NA") ## 68
sum(readings$`Duration of symptoms`=="ND") ## 51 these are observed for clin signs, but none detected 
138-68-51 ## 19 obs of duration

## Counting RNA detection/peak/duration
sum(!is.na(readings$dpi_RNAdetection))
sum(!is.na(readings$dpi_PeakRNA))
sum(!is.na(readings$`end_of_Shedding_RNA (last day pos)`))

## Counting inf virus detection/peak/duration
sum(!is.na(readings$dpi_infectious_virus))
sum(!is.na(readings$dpi_PeakVirus))
sum(!is.na(readings$end_of_Shedding_virus))



