data1<-rnorm(100)

bi<-mean(data1)

sdd<-sd(data1)

out<-data.frame(bi,sdd)

Sys.sleep(90)

save(out,file='save/save.rdata')
write.table(out,file="save/Ni.csv",sep = ",")
