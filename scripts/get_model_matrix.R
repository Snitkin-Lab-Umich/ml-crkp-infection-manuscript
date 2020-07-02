# Get binary model matrix to more easily plot results

get_model_matrix = function(matpath = 'data/combined/infection_patient-kleborate.tsv',incl_out=F){
  mat = read.delim(matpath)
  mat = data.frame(t(mat))
  if(grepl('infection',matpath)){
    outcome_name='infection'
    outcome=mat$infection
  }else if(grepl('blood',matpath)){
    outcome_name='blood'
    outcome=mat$blood
  }
  for(i in c('LOSbeforeCx','age','num_resistance_classes','num_resistance_genes','virulence_score','resistance_score')){
    if(i %in% colnames(mat)){
      mat[,i] = as.numeric(as.character(mat[,i]))
    }
  }
  # remove rows with no variation
  mat = mat[,sapply(mat,function(x) length(unique(x)) != 1)]
  if(grepl('infection',matpath)) mat = model.matrix(infection~.,dat=mat)
  if(grepl('blood',matpath)) mat = model.matrix(blood~.,dat=mat)
  colnames(mat) = gsub('1$','',colnames(mat))
  colnames(mat)[colnames(mat) %in% c("Bla_broadSHV-1")] = paste0(colnames(mat)[colnames(mat) %in% c("Bla_broadSHV-1")],'1')
  mat = mat[,2:ncol(mat)]
  mat = data.frame(cbind(ifelse(outcome=='yes',1,0),mat))
  colnames(mat)[1] = outcome_name
  if(incl_out) return(mat)
  return(mat[,2:ncol(mat)])
}
