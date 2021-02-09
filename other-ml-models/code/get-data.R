dat = read.delim('../../../../../../esnitkin/Project_Penn_KPC/Analysis/Patient_outcome/2019-10-08_ml-inf-lr/data/combined/infection_patient-kleborate.tsv')
dat = t(dat)

ltachs = read.delim('../../../../../../esnitkin/Project_Penn_KPC/Analysis/Patient_outcome/2019-10-08_ml-inf-lr/data/ltachs.tsv')
ltachs_names = rownames(ltachs)
ltachs = ltachs[,1]
names(ltachs) = ltachs_names
ltachs = ltachs[rownames(dat)]

write.csv(dat, 'data/infection_patient-kleborate.csv', row.names = F, quote = F)
write.csv(data.frame(ltach=ltachs), 'data/ltachs.csv', row.names = F, quote = F)
