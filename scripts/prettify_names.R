# make prettier names for features

library(tidyverse)
library(Hmisc)

source('scripts/get_model_matrix.R')

mat = data.frame(get_model_matrix())

prettify_names = function(x){
    x = gsub('_',' ',x) %>% 
      gsub('30yes',' (past 30 days',.) %>%
      gsub('yes','',.) %>% 
      gsub('prior admit','Admission to an LTACH in the year prior to culture',.) %>% 
      gsub('central line','Presence of a central venous catheter',.) %>%
      gsub('foley','Presence of a urinary catheter',.) %>%
      gsub('LOSbeforeCx','Length of stay before culture',.) %>%
      gsub('sexmale','Sex (male',.) %>%
      gsub('trach','Tracheostomy',.) %>%
      gsub('CHF','Congestive heart failure',.) %>%
      gsub('resp failure','Acute or chronic respiratory failure',.) %>%
      gsub('AKI','Acute kidney injury',.) %>%
      gsub('gastro','Presence of a gastrostomy tube',.) %>%
      gsub('underweight','Malnourished/underweight',.) %>%
      gsub('txp','Transplant',.) %>%
      gsub('CKD','Severe chronic kidney disease',.) %>%
      gsub('lung dz','Pulmonary disease',.) %>%
      gsub('VDRF','Ventilator-dependent respiratory failure',.) %>% 
      gsub('decub','Stage IV/V decubitus ulcer',.) %>%
      gsub('malignancy','Malignancy (solid or liquid',.) %>%
      gsub('obese','Obesity',.) %>%
      gsub('bactrim','Trimethoprim/sulfamethoxazole',.) %>%
      gsub('ctx','Ceftriaxone',.) %>%
      gsub('cipro','Ciprofloxacin',.) %>%
      gsub('dapto','Daptomycin',.) %>% 
      gsub('erta','Ertapenem',.) %>%
      gsub('flagyl','Metronidazole',.) %>%
      gsub('gent','Gentamicin',.) %>%
      gsub('levo','Levofloxacin',.) %>%
      gsub('tobra','Tobramycin',.) %>%
      gsub('zosyn','Piperacillin/tazobactam',.) %>%
      gsub('vanco','Intravenous vancomycin',.) %>%
      gsub('poly colistin','Polymyxin or colistin',.) %>% 
      gsub('AG ','Aminoglycoside ',.) %>%
      gsub('third gen ceph','Third-generation cephalosporin',.) %>%
      gsub('anti pseudo','Anti-pseudomonal antibiotic',.) %>%
      gsub('FQ','Fluoroquinolone',.) %>%
      gsub('pseudoPCN','Piperacillin/tazobactam or ceftazidime',.) %>%
      gsub('.ybt 17',' (ybt17',.) %>%
      gsub('.ybt 0',' (ybt0',.) %>%
      gsub('.ybt unknown',' (ybt unknown',.) %>%
      gsub('Yersiniabactin. ICEKp10','ICEKp10',.) %>%
      gsub('Colibactinclb.3','Colibactin (clb 3',.) %>%
      gsub('Colibactinclb.unknown','Colibactin (clb unknown',.) %>%
      gsub('Aerobactiniuc.5','Aerobactin (iuc 5',.) %>%
      gsub('^wzi','',.) %>%
      gsub('K locusKL','K locus ',.) %>%
      gsub('K locus missing genes.KL','Missing K locus ',.) %>%
      gsub('O locusO','O locus ',.) %>%
      gsub('O locus missing genes.O1.','Missing O locus O1/',.) %>%
      gsub('AGly.StrB\\|AGly.StrA|AGly.StrA\\|AGly.StrB','Aminoglycoside res (StrA & StrB',.) %>%
      gsub('Tmt.DfrA17\\|AGly.AadA5|AGly.AadA5\\|Tmt.DfrA17', '    Trimethoprim res (DfrA17) & Aminoglycoside res (AadA5',.) %>% 
      gsub('AGly\\.','Aminoglycoside res (',.) %>%
      gsub('Flq\\.','Fluoroquinolone res (',.) %>%
      gsub('ColMgrB.trunc','Truncated MgrB',.) %>%
      gsub('MLS\\.','Macrolide res (',.) %>%
      gsub('Phe\\.','Phenicol res (',.) %>%
      gsub('Rif','Rifampin res (',.) %>% 
      gsub('Sul\\.','Sulfonamide res (',.) %>% 
      gsub('^Tet','Tetracycline res (',.) %>% 
      gsub('Tmt\\.','Trimethoprim res (',.) %>% 
      gsub('Bla\\.', 'Bla (',.) %>% 
      gsub('Bla Carb', 'Carbapenemase (',.) %>% 
      gsub('Bla ESBL', 'ESBL (',.) %>% 
      gsub('Bla broad', 'Bla broad (',.) %>% 
      gsub('Bla broad \\( inhR\\.', 'Bla broad inhR (',.) %>%
      gsub('KPC.','KPC-',.) %>%
      gsub('SHV.','SHV-',.)
    x = sapply(x, function(y) ifelse(grepl('\\(',y),paste0(y,')'),y))
    x = capitalize(x)
    return(unname(x))
}


