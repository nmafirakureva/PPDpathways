library(here) #NOTE this uses petedodd/ecrins: keep versions in sync!
source(here('transmission/utilities.R'))

## example output
tt <- seq(from=0, to=10, by=0.1)  #time frame to run over
y <- runmodel(tt,parms)           #run model
head(y)

## example use of using 1st row of churn samples
parms <- revise.flow.parms(parms,smpsd,1) #change flow parameter to 1st row of samples
y <- runmodel(tt,parms)           #run model
tail(y)

## trying with tree parameters also
parms <- revise.HE.parms(parms,DR,1,arm='soc') #change flow parameter to 1st row of samples
y <- runmodel(tt,parms)           #run model
tail(y)

parms <- revise.HE.parms(parms,DR,1,arm='int') #change flow parameter to 1st row of samples
y <- runmodel(tt,parms)           #run model
tail(y)


## run test NOTE breaks
parms$staticfoi <- -1
tt <- seq(from=0, to=50, by=0.1)  #time frame to run over
y <- runmodel(tt,parms)           #run model
ydm <- output2dt(y)
## ydm[!is.na(population),unique(variable)] #all TB states

## PPD look:
pop <- ydm[!is.na(population),.(value=sum(value)),by=.(t,population)]
popt <- pop[population!='previously detained',.(value=sum(value)),by=t] #total
ggplot(popt,aes(t,value))+geom_line()+ylim(c(0,NA))

ggplot(pop[population!='previously detained'],aes(t,value,col=population))+
  geom_line()+ylim(c(0,NA))

## TB look:
tbpop <- ydm[!is.na(population) & population!='previously detained',
             .(value=sum(value)),
             by=.(t,variable)]
ggplot(tbpop,aes(t,value,col=variable))+
  geom_line()+facet_wrap(~variable,scales='free')

tbpopf <- ydm[!is.na(population) & population!='previously detained',
             .(value=sum(value)),
             by=.(t,variable,population)]
ggplot(tbpopf,aes(t,value,col=population))+
  geom_line()+facet_wrap(~variable,scales='free')

## TPT look:
tptpop <- ydm[!is.na(population) & population!='previously detained',
             .(value=sum(value)),
             by=.(t,tpt,population)]
ggplot(tptpop,aes(t,value,col=tpt))+
  geom_line()+facet_wrap(~population,scales='free')


## Other outputs:
ggplot(ydm[variable %in% others],aes(t,value,col=variable))+
  geom_line()+facet_wrap(~variable,scales='free')

## ============= diff plot  =============

parms$staticfoi <- -1            #dynamic
parms$int_time <- 50             #fix
ydb <- diffdata(parms)

ggplot(ydb[variable %in% others],aes(t,value,col=variable,lty=arm))+
  geom_line()+facet_wrap(~variable,scales='free')

ggsave(here('transmission/plots/p_illustrate_dyn.png'),w=7,h=5)

## ============= HE workflow =============

## NOTE LTBI & foi parameters
parms$parm_frac_L <- 0.2
summary(qfun(runif(1e3),hyperparms$foi))

## RUN PSA:
RES <- PSAloop(Niter=2e3,parms,smpsd,DR) #TODO update this as we go



ggplot(RES,aes(Q.soc - Q.int,int.CC-soc.CC,col=mid.notes<100 & mid.notes>30))+
  geom_point(shape=1) +
  scale_y_continuous(label=comma,'Incremental cost',limits=c(0,NA))+
  scale_x_continuous(label=comma,'Health gain',limits=c(0,NA))+
  geom_abline(intercept=0,slope=30e3,col=2)+
  theme_classic()+ggpubr::grids()+
  theme(legend.position='top')

ggsave(file=here('transmission/plots/p_CE1.png'),w=7,h=7)

RES[,mean(int.CC-soc.CC)]
RES[,mean(Q.int - Q.soc)]
RES[,mean(int.CC-soc.CC)/mean(Q.int - Q.soc)]

RES[mid.notes<100 & mid.notes>30,mean(int.CC-soc.CC)]
RES[mid.notes<100 & mid.notes>30,mean(Q.int - Q.soc)] #BUG?
RES[mid.notes<100 & mid.notes>30,mean(int.CC-soc.CC)/mean(Q.int - Q.soc)]

## make CEAC
ll <- seq(from=0,to=50e3,by=100)
cc <- make.ceac(RES[,.(Q=dQ,P=int.CC-soc.CC)],ll)
ccr <- make.ceac(RES[,## mid.notes<100 & mid.notes>30,
                     .(Q=dQ,P=int.CC-soc.CC)],
                 ll)
D <- data.table(x=ll,all=cc,restricted=ccr)
D <- melt(D,id='x')

## plot
xpad <- 2.5e3
GP <- ggplot(D,aes(x,value,col=variable)) +
  geom_line() +
  theme_classic() +
  theme(legend.position = 'top')+
  scale_y_continuous(label=percent,
                     limits=c(0,1),
                     expand = c(0,0))+
  scale_x_continuous(label=comma,
                     limits=c(0,max(D$x)),
                     expand = c(0,xpad))+
  ## scale_color_colorblind()+
  scale_colour_manual(values=cbPalette)+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD/DALY)')+
  ggpubr::grids()
GP


ggsave(file=here('transmission/plots/p_CEAC1.png'),w=7,h=7)
