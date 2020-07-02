# get isolate locations (for Figure S11)

library(lubridate)

# Change names of locations
change_loc_names = function(metadata_location_info) {
  
  sapply(metadata_location_info, function(x){
    if(grepl('Canada',x)){
      'CAN'
    } else if(grepl('Australia',x)){
      'AUS'
    } else if(grepl('Israel',x) | grepl('Lebanon',x) | grepl('United Arab Emirates',x)
              | grepl('Saudi Arabia',x)){
      'MIDE'
    } else if(grepl('Brazil',x) | grepl('Colombia',x) | grepl('Guatemala',x) | grepl('Honduras',x) 
              | grepl('Argentina',x)) {
      'SA'
    } else if(grepl('Denmark',x) | grepl('Finland',x) | grepl('Greece',x) | grepl('Italy',x) | grepl('Poland',x)
              | grepl('Netherlands',x) | grepl('Spain',x) | grepl('Norway',x) | grepl('Germany',x)
              | grepl('France',x) | grepl('United Kingdom',x) | grepl('Portugal',x) | grepl('BSAC',x)) {
      'EUR'
    } else if(grepl('Hong Kong',x) | grepl('India',x) | grepl('Indonesia',x) | grepl('Korea',x) | grepl('Malaysia',x) 
              | grepl('Philippines',x) | grepl('Singapore',x) | grepl('Thailand',x) | grepl('China',x) 
              | grepl('Taiwan',x) | grepl('Nepal',x)) {
      'ASIA'
    } else if(grepl('South Africa',x)) {
      'AFR'
    }else if(grepl('California',x) | grepl('CA',x)){
      'US-CA'
    } else if(grepl('Arizona',x) | grepl('Colorado',x) | grepl('New Mexico',x) 
              | grepl('Oregon',x) | grepl('Washington',x) | grepl('WA',x)) { # MIGHT BE AN ISSUE WITH WASHINGTON DC
      'US-W'
    } else if(grepl('New Jersey',x) | grepl('New York',x) | grepl('NY',x) | grepl('Pennsylvania',x) | grepl('New Hampshire',x) | grepl('Boston',x) | grepl('Massachusetts',x) | grepl('MA',x)
              | grepl('Rhode Island',x) | grepl('Connecticut',x)) {
      'US-NE'
    } else if(grepl('North Carolina',x) | grepl('Florida',x) | grepl('Georgia',x) | grepl('West Virginia',x) | grepl('Virginia',x) 
              | grepl('Texas',x) | grepl('North Carolina',x) | grepl('Maryland',x) | grepl('TX',x)
              | grepl('San Antonio',x) | grepl('Delaware',x) | grepl('MD',x) | grepl('Orlando',x) 
              | grepl('Walter Reed',x) | grepl('Houston',x) | grepl('DC',x) | grepl('Louisville',x)
              | grepl('Tampa',x) | grepl('FL',x) | grepl('NIH',x)) {
      'US-S'
    } else if(grepl('Illinois',x) | grepl('Michigan',x) | grepl('Missouri',x) | grepl('North Dakota',x) 
              | grepl('South Dakota',x) | grepl('Massachusetts',x) | grepl('IA',x) | grepl('OH',x) 
              | grepl('Mid Atlantic',x) | grepl('Chicago',x) | grepl('HHS Region 5',x)) {
      'US-MW'
    } else if(x == '' | grepl('unknown',x) | grepl('not applicable',x) | grepl('not collected',x)) {
      'unknown'
    } else {
      paste(x)
    }
  })
}


#Change date at end of string after underscore to number of years
date_to_decimal_year = function(names_with_date,ymd='ymd') {
  sapply(1:length(names_with_date), function(x){
    if(ymd == 'ymd'){
      date = sub("[^_]*$", decimal_date(ymd(regmatches(names_with_date[x],
                                                       gregexpr("[^_]*$", names_with_date[x])),quiet=T)),
                 names_with_date[x])
      if(is.na(date)){names_with_date[x]} else {date}
    } else if(ymd =='mdy') {
      date = sub("[^_]*$", decimal_date(mdy(regmatches(names_with_date[x],
                                                       gregexpr("[^_]*$", names_with_date[x])),quiet=T)),
                 names_with_date[x]) 
      if(is.na(date)){names_with_date[x]} else {date}
    } else if(ymd =='dmy') {
      date = sub("[^_]*$", decimal_date(dmy(regmatches(names_with_date[x],
                                                       gregexpr("[^_]*$", names_with_date[x])),quiet=T)),
                 names_with_date[x]) 
      if(is.na(date)){names_with_date[x]} else {date}
    }
  })}


get_locations = function(genome_ids){
  
  sapply(genome_ids, function(x){
    if(sum(x %in% sample_key$isolate_no) > 0) {
      
      loc = sample_key$LTACH_col[sample_key$isolate_no == x]
      loc = change_loc_names(loc)
      
    } else if(sum(x %in% limbago_data$Run_s) > 0) {
      
      loc = paste(limbago_data$geo_loc_name_s[limbago_data$Run_s == x])
      loc = change_loc_names(loc)
      
    } else if(sum(x %in% srr_biosample$SRR) > 0) {
      loc = paste(patric_data$Geographic.Location[patric_data$BioSample.Accession == as.character(srr_biosample$BioSamp[srr_biosample$SRR == x])])
      loc = change_loc_names(loc)
      
    } else if(x == 'gi_661922017_gb_CP008827'){
      
      loc = 'NIH'
      loc = change_loc_names(loc)
      
    } else if(sum(x %in% sample_key$isolate_no) == 0 & sum(x %in% patric_data$Genome.ID) == 0){
      
      paste('')
    }
  })
  
} 
