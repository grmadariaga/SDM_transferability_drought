library(dplyr)
library(ggplot2)
library(lme4)
library(MuMIn)
library(randomForest)
library(tidyverse)
library(lmerTest)
library(AICcmodavg)
library(report)
library(optimx)
library(performance)
library(AICcmodavg)
library(jcolors)
library(see)
library(emmeans)
library(interactions)


###Create a variable with the name of the folderwhere the files were downloaded

Myfiles<- "Your/folder/"

setwd(Myfiles)

model_df <- readRDS("df_models_all_or.rds")

Knames <- readRDS("Knames.rds")

#disperse <- read.csv("R_disperse.csv")

#disperse_rec <- disperse%>%
#  mutate(across(5:44, ~ ifelse(. %in% c(0), 0, 1)))
#pco_12 <- readRDS("pco.12.rds")
#pco_12$Genus <- rownames(pco_12)

Freshweco<-read.csv("Freshwatereco.csv")

#density <- read.csv("density_train.csv", header = F)
#density<- density[, c(1, 136)]
#names(density)<- c("Species", "density")


Freshweco<- Freshweco%>%
  mutate(tolerance= case_when(
    tol==1 ~"tol",
    mtol==1 ~ "mtol",
    msen==1 ~"msen", 
    sen==1~ "sen", 
    hsen==1 ~ "hsen"
  ))%>%
  mutate(Genus= str_extract(Taxon, "^(\\w+)"),
         spname=str_trim(str_extract(Taxon, "\\s(\\w+)")))%>%
  mutate(Family= str_to_title(Family))%>%
  #distinct(Family, .keep_all = TRUE)%>%
  as.data.frame()

model_df <- model_df %>%
  mutate(Species2=Species)%>%
  separate(Species2, into = c('Genus', 'spname'), sep = "\\.") %>%
  mutate(AUC_diff = AUC_tr-AUC_ev)%>%
  mutate(TSS_diff=TSS_tr-TSS_ev)%>%
  group_by(Species, model)%>%
  mutate(Var_AUC_tr= var(AUC_tr),
         Var_TSS_tr= var(TSS_tr),
         Var_AUC_ev= var(AUC_ev),
         Var_TSS_ev= var(TSS_ev),
         Var_AUC_diff= var(abs(AUC_diff)),
         Var_TSS_diff= var(abs(TSS_diff)))%>%
  mutate(Fold= as.factor(Fold))%>%
  as.data.frame()

#model_df_o <- left_join(model_df, disperse, by= "Genus")
#model_df_pco <- left_join(model_df, pco_12, by= "Genus")
#model_df_bin <- left_join(model_df_pco, disperse_rec, by= "Genus")
model_df_fwei <- left_join(model_df, Freshweco%>%
                             filter(!Taxon=="Asellus aquaticus (karstic type)")%>%
                             filter(!Taxon=="Oulimnius tuberculatus Lv.")%>%
                             filter(!Taxon=="Limnius volckmari Lv.")%>%
                             filter(!Taxon=="Orectochilus villosus Lv.")%>%
                             distinct(Genus, spname, .keep_all = TRUE),
                           by= c("Genus", "spname"))%>%
  mutate(
    AUC_diff_sc= scale(AUC_diff),
    TSS_diff_sc= scale(TSS_diff)
  )


filt <-model_df_fwei%>%filter(!is.na(tolerance))

#filt<- filt%>%filter(!Species=="Aphelocheirus.aestivalis")

#filt <- filt%>%filter(Species!="Calopteryx.virgo")


filt$gap_perc <- (filt$AUC_diff*100)/filt$AUC_tr


filt%>%
  pull(AUC_tr)%>%
  mean()

filt%>%
  pull(AUC_tr)%>%
  sd()



cor.test(filt$AUC_tr, filt$TSS_tr,  method = "pearson", use = "complete.obs")



wilcox_data <- pivot_longer(filt, cols= c("AUC_tr", "AUC_ev"), names_to = "Type", values_to = "AUC_tr_ev")

wilcox.test(AUC_tr_ev ~ Type, data= wilcox_data)

table(filt$model)

qqnorm(filt$AUC_diff)
qqline(filt$AUC_diff)
shapiro.test(filt$AUC_diff)



kruskal.test(AUC_tr~model, data = filt) ###sig difference in models BEE ac to model
kruskal.test(AUC_tr~tolerance, data = filt) ###no difference in models BEE ac to tolerance

dunn.test::dunn.test(filt$AUC_tr, filt$model, method = "BH")
dunn.test::dunn.test(filt$AUC_tr, filt$tolerance, method = "BH")

kruskal.test(AUC_ev~model, data = filt) ###sig difference in models BEE ac to model
kruskal.test(AUC_ev~tolerance, data = filt)

filt%>%
  filter(AUC_ev>0.5)%>%
  nrow(.)

###426 models better than random ~81%

filt%>%
  filter(AUC_ev>=0.7)%>%
  nrow(.)
#### 205 good models ~ 39%

filt%>%
  pull(AUC_ev)%>%
  mean()


filt%>%
  pull(AUC_ev)%>%
  sd()


cor.test(filt$AUC_ev, filt$TSS_ev,  method = "pearson", use = "complete.obs")


filt%>%
  pull(AUC_diff)%>%
  mean()


filt%>%
  pull(AUC_diff)%>%
  sd()

filt%>%
  pull(gap_perc)%>%
  mean()


lm(AUC_diff~AUC_tr, data=filt)%>%summary()
lm(AUC_diff~AUC_tr, data=filt)%>%summary()


Preval_df <-read.csv("prevalence.csv")[,-1]

filt_prev<-left_join(filt, Preval_df, by="Species")

filt_prev$Percprev <- filt_prev$Presences/146

lm(AUC_diff~Percprev, data=filt_prev)%>%summary()### high significance, low rsquared
lm(AUC_diff~Presences, data=filt_prev)%>%summary()

kruskal.test(AUC_diff~model, data = filt) ###sig difference in AUCgap ac to model
dunn.test::dunn.test(filt$AUC_diff, filt$model, method = "BH")



filt%>%
  filter(model=="ssn")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(model=="ssn")%>%
  pull(AUC_diff)%>%
  sd()


filt%>%
  filter(model=="RF")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(model=="RF")%>%
  pull(AUC_diff)%>%
  sd()

filt%>%
  filter(model=="glm")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(model=="glm")%>%
  pull(AUC_diff)%>%
  sd()

filt%>%
  filter(model=="Maxent")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(model=="Maxent")%>%
  pull(AUC_diff)%>%
  sd()


kruskal.test(AUC_diff~tolerance, data = filt) ###sig difference in AUCgap ac to tolerance
dunn.test::dunn.test(filt$AUC_diff, filt$tolerance, method = "BH")




filt%>%
  filter(tolerance=="tol")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(tolerance=="tol")%>%
  pull(AUC_diff)%>%
  sd()


filt%>%
  filter(tolerance=="mtol")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(tolerance=="mtol")%>%
  pull(AUC_diff)%>%
  sd()

filt%>%
  filter(tolerance=="msen")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(tolerance=="msen")%>%
  pull(AUC_diff)%>%
  sd()


filt%>%
  filter(tolerance=="sen")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(tolerance=="sen")%>%
  pull(AUC_diff)%>%
  sd()

lm(AUC_diff~Species, data=filt_prev)%>%summary()###higher significance, very high r squared





#lm(AUC_diff~Species, data=filt)%>%summary()
#lm(AUC_diff~AUC_tr, data=filt)%>%summary()


###fit all relevant models

Alg_off <- lmer(AUC_diff~ model+ (1|Species), 
                control= lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                data=filt%>% 
                  filter(!is.na(tolerance)),
                REML= FALSE)

Tolerance_off <-  lmer(AUC_diff~ tolerance+ (1|Species), 
                       control= lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                       data=filt%>% 
                         filter(!is.na(tolerance)),
                       REML= FALSE)

Tol.Alg_off <- lmer(AUC_diff~ model+ tolerance + (1|Species), 
                    control= lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                    data=filt%>% 
                      filter(!is.na(tolerance)),
                    REML= FALSE)

TolxAlg_off <- lmer(AUC_diff~ model* tolerance+(1|Species), 
                    control= lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                    data=filt%>% 
                      filter(!is.na(tolerance)),
                    REML= FALSE)

Rand.mod_off <-lmer(AUC_diff~ (1|Species), 
                    control= lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                    data=filt%>% 
                      filter(!is.na(tolerance)),
                    REML= FALSE)






mods.list <- list(Alg_off, Tolerance_off, Tol.Alg_off,  TolxAlg_off, Rand.mod_off)

mods.names <- c("Model", "Tolerance", "Tol+Mod", "Tol*Mod", "Random")

aictab(cand.set = mods.list, modnames = mods.names)


Alg <- lmer(AUC_diff~ model+ (1|Species), 
            control= lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
            data=filt%>% 
              filter(!is.na(tolerance)),
            REML= TRUE)

Tolerance <-  lmer(AUC_diff~ tolerance+ (1|Species), 
                   control= lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                   data=filt%>% 
                     filter(!is.na(tolerance)),
                   REML= TRUE)


Tol.Alg <- lmer(AUC_diff~ model+ tolerance + (1|Species), 
                control= lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                data=filt%>% 
                  filter(!is.na(tolerance)),
                REML= TRUE)


TolxAlg <- lmer(AUC_diff~ model* tolerance+(1|Species), 
                control= lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                data=filt%>% 
                  filter(!is.na(tolerance)),
                REML= TRUE)


Rand.mod <-lmer(AUC_diff~ 1+ (1|Species), 
                control= lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                data=filt%>% 
                  filter(!is.na(tolerance)),
                REML= TRUE)


###Check residuals


plot(fitted(TolxAlg),resid(TolxAlg))
acf(resid(TolxAlg))

qqnorm(ranef(TolxAlg)$Species[,"(Intercept)"])##Looks good, small outliers
qqline(ranef(TolxAlg)$Species[,"(Intercept)"])##Looks good, small outliers


compare_performance(Alg, Tolerance, Tol.Alg, TolxAlg, Rand.mod, rank = TRUE, metrics = c("R2", "AIC"))


TolxAlgAnova <- lmer(AUC_diff~ model* tolerance+(1|Species), 
                     control= lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                     contrasts=list(tolerance=contr.sum, model=contr.sum),
                     data=filt%>% 
                       filter(!is.na(tolerance)),
                     REML= TRUE)

car::Anova(TolxAlg, type=3, test.statistic="Chisq")

plot(TolxAlgAnova)

emmip(TolxAlg, model ~ tolerance, CIs=TRUE,
      dodge=0.5,
      dotarg= list(size=6, shape=19),
      linearg=list(size=3, linetype=1),
      CIarg= list(lwd=6, alpha=0.5, position= "dodge"))+
  labs(x="Tolerance", y= "Predicted AUCgap", col="Method")+
  theme_light()+
  scale_fill_manual(labels = c("GLM", "MAXENT", "RF", "SSN"), values = c("#66CD00","#009ACD", "#FF4500","#9A32CD"  ))+
  scale_color_manual(labels = c("GLM", "MAXENT", "RF", "SSN"), values = c("#66CD00","#009ACD", "#FF4500","#9A32CD"  ))+
  theme(axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size=18, face="bold"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        axis.line = element_line(size = 1.2)
  )

model_emm <- emmeans(TolxAlgAnova, pairwise ~ model | tolerance)

pairs(model_emm, adjust="BH")

color.comp <- rep(c("#66CD00","#009ACD","#66CD00", "#FF4500","#66CD00","#9A32CD", "#009ACD", "#FF4500", "#009ACD", "#9A32CD", "#FF4500", "#9A32CD"), 4)
lab.cols <- c("#66CD00", "#009ACD", "#FF4500", "#9A32CD") 

my.aes <- list(label = list(size = 5), color= lab.cols,
               point = list(shape = "circle", color=color.comp), 
               segment = list(linewidth = 2), color=color.comp)

mod.colors <- c(glm="#66CD00", Maxent="#009ACD", RF="#FF4500",ssn="#9A32CD"  )

pwpp(model_emm, sort=FALSE, 
     aes = my.aes, 
     ylab = c(""))+
  theme_light()+
  geom_point(aes(size = 3)) + 
  theme(text = element_text(size = 25),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=19)) +
  scale_color_manual(values=mod.colors)+
  scale_size(guide = F)+
  facet_grid(factor(tolerance, levels= c("tol", "mtol", "msen", "sen"))~., labeller = labeller(c("msen", "mtol", "sen", "tol")))+
  scale_y_discrete(labels=c("GLM", "MAXENT", "RF", "SSN"))




####GRAPHS
gd <- filt %>% 
  group_by(model) %>% 
  summarise(
    AUC_tr = mean(AUC_tr),
    AUC_ev = mean(AUC_ev),
    TSS_tr = mean(TSS_tr),
    TSS_ev= mean(TSS_ev),
    AUC_diff= mean(AUC_diff),
    TSS_diff=mean(TSS_diff)
  )

tol_order <- c("tol", "mtol", "msen", "sen")




library(tidyr)
library(ggplot2)


ggplot(filt, aes(x = AUC_tr, y = AUC_diff)) + 
  geom_point(col="gray2") +
  stat_smooth(method = "lm", col = "blue")+
  xlab("Accuracy Period1")+
  ylab("Accuracy Gap")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45))+
  #scale_x_continuous(breaks=c(2007:2021), limits=c(2006, 2022))+
  theme_light()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x=element_text(size=16, hjust = 1), 
        axis.text.y =element_text(size = 14), 
        axis.title = element_text(size=18),
        axis.line = element_line(size = 1.2))



ggplot(filt_prev, aes(x = Presences/146, y = AUC_diff)) + 
  geom_point(col="gray2", position = "jitter") +
  stat_smooth(method = "lm", col = "blue")+
  xlab("Prevalence")+
  ylab("Accuracy Gap")+
  scale_x_continuous(labels = scales::percent)+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45))+
  #scale_x_continuous(breaks=c(2007:2021), limits=c(2006, 2022))+
  theme_light()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x=element_text(size=16, hjust = 1), 
        axis.text.y =element_text(size = 14), 
        axis.title = element_text(size=18),
        axis.line = element_line(size = 1.2))

okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

Before <- ggplot(filt%>%filter(
  !is.na(tolerance)),
  aes(x = AUC_tr,
      y = TSS_tr, 
      color= model,
      fill=model)) +
  geom_point(alpha = .5, size=2.5, shape= 8) +
  geom_point(data = gd, size = 8, shape= 19, alpha= 0.85)+
  theme_minimal() +
  labs(x = "AUC", y = "TSS",
       #title = "Accuracy metrics (AUC and TSS) of models before extreme events",
       col= "Method")+
  ylim(0.25,1)+
  xlim(0.25,1)+
  scale_fill_manual(labels = c("GLM", "MAXENT", "RF", "SSN"), values = c("#66CD00","#009ACD", "#FF4500","#9A32CD"  ))+
  scale_color_manual(labels = c("GLM", "MAXENT", "RF", "SSN"), values = c("#66CD00","#009ACD", "#FF4500","#9A32CD"  ))+
  theme(axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size=18, face="bold"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        axis.line = element_line(size = 1.2)
  )+
  guides(fill=guide_legend(title="Method"))+
  theme(
    #panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #legend.background = element_rect(fill='transparent'),
    #legend.box.background = element_rect(fill='transparent')
  )

Before
#ggsave("Images/Manuscript/Before_EE.png", Before, bg="transparent")

After <- ggplot(filt%>%filter(
  !is.na(tolerance)),
  aes(x = AUC_ev,
      y = TSS_ev, 
      color= model,
      fill=model)) +
  geom_point(alpha = .5, size=2.5, shape= 8) +
  geom_point(data = gd, size = 8, shape= 19, alpha= 0.85)+
  theme_minimal() +
  labs(x = "AUC", y = "TSS",
       #title = "Accuracy metrics (AUC and TSS) of models after extreme events",
       col= "Method")+
  ylim(0.25,1)+
  xlim(0.25,1)+
  scale_fill_manual(labels = c("GLM", "MAXENT", "RF", "SSN"), values = c("#66CD00","#009ACD", "#FF4500","#9A32CD"  ))+
  scale_color_manual(labels = c("GLM", "MAXENT", "RF", "SSN"), values = c("#66CD00","#009ACD", "#FF4500","#9A32CD"  ))+
  theme(axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size=18, face="bold"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        axis.line = element_line(size = 1.2)
  )+
  guides(fill=guide_legend(title="Method"))+
  theme(
    #panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #legend.background = element_rect(fill='transparent'),
    #legend.box.background = element_rect(fill='transparent')
  )

After

#ggsave("Images/Manuscript/After_EE.png", After, bg="transparent")

filt_cond <- filt%>%rename(Drought.free=AUC_tr,
                           Drought.influenced=AUC_ev)%>%
  pivot_longer(cols = c("Drought.free", "Drought.influenced"), 
               names_to="Conditions")

filt_cond$alpha<- ifelse(filt_cond$Conditions=="Drought.free", 0.8, 0.6)

gd_cond  <- gd%>%rename(Drought.free=AUC_tr,
                        Drought.influenced=AUC_ev)%>%
  pivot_longer(cols = c("Drought.free", "Drought.influenced"), 
               names_to="Conditions")

model_names <- c(
  'glm'="GLM",
  'ssn'="SSN",
  'Maxent'="MAXENT",
  'RF'="RF"
)

panelmodel <-ggplot(filt_cond%>%filter(
  !is.na(tolerance)),
  aes(y = value,
      #group=model,
      x = Conditions, 
      color= model,
      fill=model, 
      alpha=Conditions)) +
  geom_boxplot(linewidth=1.5)+
  geom_point(data = gd_cond, size = 8, shape= 19, alpha= 0.85)+
  theme_light() +
  facet_wrap(~model, labeller = as_labeller(model_names), nrow=1 )+
  labs(y = "AUC",
       #title = "Accuracy metrics (AUC and TSS) of models before extreme events",
       col= "Method")+
  ylim(0.25,1)+
  scale_x_discrete(labels=c("Drought \nfree", "Drought \ninfluenced"))+
  scale_fill_manual(labels = c("GLM", "MAXENT", "RF", "SSN"), values = c("#66CD00","#009ACD", "#FF4500","#9A32CD"  ))+
  scale_color_manual(labels = c("GLM", "MAXENT", "RF", "SSN"), values = c("#66CD00","#009ACD", "#FF4500","#9A32CD"  ))+
  scale_alpha_manual(values = c(0.4, 0.15), guide = FALSE)+
  theme(axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16,face="bold"),
        axis.title.x=element_blank(),
        strip.text = element_text(size=18, face="bold"),
        plot.title = element_text(size=18, face="bold"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        # axis.line = element_line(size = 1.2)
  )+
  guides(fill=guide_legend(title="Method"))+
  theme(
    #panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #legend.background = element_rect(fill='transparent'),
    #legend.box.background = element_rect(fill='transparent')
  )

ggsave("Panel_methods.png",
       panelmodel, width = 12, height = 8,
       units = "in", dpi = 300)





filt_long <- filt%>%
  pivot_longer(cols = c("AUC_tr", "AUC_ev"),
               names_to = "Dataset",
               values_to = "Acc")



AUC_gap <-filt_prev %>% 
  ggplot(aes(x = reorder(Species, -Presences) , y = abs(AUC_diff), col= model)) +
  stat_summary(fun = mean,
               geom = "pointrange",
               position = position_dodge(0.8),
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)))+
  theme_light() +
  labs(x = "Species", 
       y = "Accuracy gap (AUC)", 
       title = "Difference in accuracy (AUC) between models before and after extreme events",
       col= "Method:")+
  scale_fill_manual(labels = c("GLM", "MAXENT", "RF", "SSN"), values = c("#66CD00","#009ACD", "#FF4500","#9A32CD"  ))+
  scale_color_manual(labels = c("GLM", "MAXENT", "RF", "SSN"), values = c("#66CD00","#009ACD", "#FF4500","#9A32CD"  ))+
  facet_grid(~factor(tolerance, levels= c("tol", "mtol", "msen", "sen")), scales="free_x", space="free_x")+
  theme(axis.text.x=element_text(angle=70,size=11, hjust = 1),
        axis.text.y=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size=18, face="bold"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16))+
  guides(fill=guide_legend(title="Method:"))+
  theme(axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size=18, face="bold"),
        strip.text.x=element_text(size=14, face = "bold"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        axis.line = element_line(size = 1.2)
  )+
  theme(
    #panel.background = element_rect(fill='transparent'),
    #plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    #legend.background = element_rect(fill='transparent'),
    #legend.box.background = element_rect(fill='transparent')
  )

AUC_gap








group.colors <- c(tol = "#0072B2", mtol = "#D55E00", msen = "#CC79A7", sen = "#009E73")



ggplot(filt%>%filter(
  !is.na(tolerance)),
  aes(x = factor(model),
      y = AUC_diff, 
      fill= model,
      col=model)) +
  geom_violinhalf(alpha=0.2, linewidth=1.2, width= 1.3, trim=FALSE) +
  geom_jitter(alpha=0.4, size=0.4)+
  stat_summary(fun = mean,
               size=1,
               geom = "pointrange",
               color="black",
               linewidth=1.2,
               position = position_nudge(-0.08,0),
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x))+
  theme_minimal() +
  labs(x = "Method", y = "AUC accuracy gap", 
       #title = "Accuracy gap (AUC) per model across species tolerances", 
       fill="Model", col="Model")+
  scale_fill_manual(values=mod.colors, labels= c("GLM", "Maxent", "RF", "SSN"))+
  scale_color_manual(values=mod.colors, labels= c("GLM", "Maxent", "RF", "SSN"))+
  labs(color = "Method", fill="Method")+
  scale_y_continuous(breaks= c(-0.2, 0, 0.2, 0.4, 0.6, 0.8))+
  scale_x_discrete(labels=c("GLM", "Maxent", "RF", "SSN"))+
  theme(axis.text.x=element_text(size=18, angle = 0, hjust= -0.1),
        axis.text.y=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(size = 1.2))


Tolerance_accuracy <- ggplot(filt%>%filter(
  !is.na(tolerance)),
  aes(x = factor(tolerance, levels= tol_order),
      y = AUC_diff, 
      fill= tolerance,
      col=tolerance)) +
  geom_violinhalf(alpha=0.2, linewidth=1.2, width= 1.3, trim=FALSE) +
  geom_jitter(alpha=0.4, size=0.4)+
  stat_summary(fun = mean,
               size=1,
               geom = "pointrange",
               color="black",
               linewidth=1.2,
               position = position_nudge(-0.08,0),
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x))+
  theme_minimal() +
  labs(x = "Tolerance score", y = "AUC accuracy gap", 
       #title = "Accuracy gap (AUC) per model across species tolerances", 
       fill="Tolerance score", col="Tolerance score")+
  scale_fill_manual(breaks= c("tol", "mtol", "msen", "sen"), values=group.colors)+
  scale_color_manual(breaks= c("tol", "mtol", "msen", "sen"),values=group.colors)+
  scale_y_continuous(breaks= c(-0.2, 0, 0.2, 0.4, 0.6, 0.8))+
  theme(axis.text.x=element_text(size=18, angle = 0, vjust= 1.4),
        axis.text.y=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(size = 1.2))

Tolerance_accuracy


#ggsave("Images/Manuscript/Tolerance.png", Tolerance_accuracy, bg="transparent")


AUC_model_tol<- ggplot(filt%>%filter(
  !is.na(tolerance)),
  aes(x = factor(tolerance, levels= tol_order),
      y = AUC_ev, 
      color = model, fill = model)) +
  geom_boxplot(alpha=0.2, size=1.2) +
  stat_summary(fun=mean, geom="point", 
               shape=18, size=5, position = position_dodge2(width = 0.8))+
  theme_minimal() +
  labs(x = "Tolerance", y = "AUC accuracy gap", title = "Accuracy gap (AUC) per model across species tolerances",col= "Method", fill="Method")+
  scale_fill_manual(labels = c("GLM", "MAXENT", "RF", "SSN"), values = c("#66CD00","#009ACD", "#FF4500","#9A32CD"  ))+
  scale_color_manual(labels = c("GLM", "MAXENT", "RF", "SSN"), values = c("#66CD00","#009ACD", "#FF4500","#9A32CD"  ))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
  theme(axis.text.x=element_text(size=16),
        axis.text.y=element_text(size=16),
        axis.title=element_text(size=16,face="bold"),
        plot.title = element_text(size=18, face="bold"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.key.size = unit(1.5, "cm"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(size = 1.2))


AUC_model_tol

#ggsave("Images/Manuscript/Interaction.png", AUC_model_tol, bg="transparent")


#ggsave("Images/Manuscript/Method.png", AUCgap_model, bg="transparent")



B_A_auc <-ggplot(data=filt_long,
                 aes(x = Dataset,
                     y = Acc, 
                     color = model, fill = model)) +
  geom_boxplot(alpha=0.2, size=1.2) +
  stat_summary(fun=mean, geom="point", 
               shape=18, size=5, position = position_dodge2(width = 0.8))+
  theme_minimal()+
  labs(x = "", y = "Accuracy (AUC)", col= "Method", fill="Method")+
  scale_x_discrete(limits= c("AUC_tr", "AUC_ev"),labels = c("", ""))+
  scale_fill_manual(labels = c("GLM", "MAXENT", "RF", "SSN"), values = c("#66CD00","#009ACD", "#FF4500","#9A32CD"  ))+
  scale_color_manual(labels = c("GLM", "MAXENT", "RF", "SSN"), values = c("#66CD00","#009ACD", "#FF4500","#9A32CD"  ))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, color="gray56"))+
  theme(axis.text.x=element_text(size=16, color="gray56"),
        axis.text.y=element_text(size=16, color="gray56"),
        axis.title=element_text(size=16,face="bold", color="gray56"),
        plot.title = element_text(size=18, face="bold", color="gray56"),
        legend.title = element_text(face= "bold", size=16, color="gray56"),
        legend.text = element_text(face="bold", size=16, color="gray56"),
        legend.key.size = unit(1.5, "cm"),
        rect = element_rect(fill = "transparent"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth = 1.8, color="gray56"), axis.line = element_line(size = 1.2, color="gray56"))

ggsave("G:/RESIST/Graciela/Presentations/CRC_poster/BA_AUC.png", B_A_auc, bg="transparent")


ggplot(data=wilcox_data,
       aes(x = Type,
           y = AUC_tr_ev)) +
  geom_boxplot(alpha=0.2, size=1.2) +
  geom_jitter(alpha=0.2)+
  theme_minimal()+
  scale_x_discrete(breaks=c("AUC_tr","AUC_ev"), limits=c("AUC_tr","AUC_ev"), labels= c("Period1", "Period2"))+
  labs(title = "",
       x = "Dataset",
       y = "AUC")+
  theme_light()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x=element_text(size=16, hjust = 1), 
        axis.text.y =element_text(size = 14), 
        axis.title = element_text(size=18),
        axis.line = element_line(size = 1.2))





filt %>%
  distinct(Species, tolerance) %>%
  group_by(tolerance) %>%
  summarise(species_count = n_distinct(Species))

table(filt$tolerance)



ggplot(filt%>%filter(
  !is.na(tolerance)),
  aes(x = factor(model),
      y = AUC_diff, 
      fill= model,
      col=model)) +
  geom_violinhalf(alpha=0.2, linewidth=1.2, width= 1.3, trim=FALSE) +
  geom_jitter(alpha=0.4, size=0.4)+
  stat_summary(fun = mean,
               size=1,
               geom = "pointrange",
               color="black",
               linewidth=1.2,
               position = position_nudge(-0.08,0),
               fun.max = function(x) mean(x) + sd(x),
               fun.min = function(x) mean(x) - sd(x))+
  theme_minimal() +
  labs(x = "Model", y = "AUC accuracy gap", 
       #title = "Accuracy gap (AUC) per model across species tolerances", 
       fill="Model", col="Model")+
  scale_fill_manual(values=mod.colors, labels= c("GLM", "Maxent", "RF", "SSN"))+
  scale_color_manual(values=mod.colors, labels= c("GLM", "Maxent", "RF", "SSN"))+
  labs(color = "Method", fill="Method")+
  scale_y_continuous(breaks= c(-0.2, 0, 0.2, 0.4, 0.6, 0.8))+
  scale_x_discrete(labels=c("GLM", "Maxent", "RF", "SSN"))+
  theme(axis.text.x=element_text(size=18, angle = 0, hjust= -0.1),
        axis.text.y=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(size = 1.2))



##Extra computations

(filt_prev%>%
    filter(tolerance=="msen")%>%
    pull(Presences)%>%
    mean())/208


filt%>%
  filter(model=="glm")%>%
  filter(tolerance=="tol")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(model=="glm")%>%
  filter(tolerance=="mtol")%>%
  pull(AUC_diff)%>%
  mean()


filt%>%
  filter(model=="glm")%>%
  filter(tolerance=="msen")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(model=="glm")%>%
  filter(tolerance=="sen")%>%
  pull(AUC_diff)%>%
  mean()


filt%>%
  filter(model=="RF")%>%
  filter(tolerance=="tol")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(model=="RF")%>%
  filter(tolerance=="mtol")%>%
  pull(AUC_diff)%>%
  mean()


filt%>%
  filter(model=="RF")%>%
  filter(tolerance=="msen")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(model=="RF")%>%
  filter(tolerance=="sen")%>%
  pull(AUC_diff)%>%
  mean()



filt%>%
  filter(model=="Maxent")%>%
  filter(tolerance=="tol")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(model=="Maxent")%>%
  filter(tolerance=="mtol")%>%
  pull(AUC_diff)%>%
  mean()


filt%>%
  filter(model=="Maxent")%>%
  filter(tolerance=="msen")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(model=="Maxent")%>%
  filter(tolerance=="sen")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(model=="ssn")%>%
  filter(tolerance=="tol")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(model=="ssn")%>%
  filter(tolerance=="mtol")%>%
  pull(AUC_diff)%>%
  mean()


filt%>%
  filter(model=="ssn")%>%
  filter(tolerance=="msen")%>%
  pull(AUC_diff)%>%
  mean()

filt%>%
  filter(model=="ssn")%>%
  filter(tolerance=="sen")%>%
  pull(AUC_diff)%>%
  mean()


model_df_fwei%>%filter(tolerance=="tol")%>%distinct(Species)%>%nrow()

model_df_fwei%>%filter(tolerance=="mtol")%>%distinct(Species)%>%nrow()

model_df_fwei%>%filter(tolerance=="msen")%>%distinct(Species)%>%nrow()

model_df_fwei%>%filter(tolerance=="sen")%>%distinct(Species)%>%nrow()



emmip_tolalg <-emmip(TolxAlg, model ~ tolerance, CIs=TRUE,
                     dodge=0.5,
                     dotarg= list(size=6, shape=19),
                     linearg=list(size=3, linetype=1),
                     CIarg= list(lwd=5, alpha=0.5, position= "dodge"))+
  labs(x="Tolerance", y= "Transferability", col="Method")+
  theme_light()+
  scale_fill_manual(labels = c("GLM", "MAXENT", "RF", "SSN"), values = c("#66CD00","#009ACD", "#FF4500","#9A32CD"  ))+
  scale_color_manual(labels = c("GLM", "MAXENT", "RF", "SSN"), values = c("#66CD00","#009ACD", "#FF4500","#9A32CD"  ))+
  theme(axis.text.x=element_text(size=16, colour = "gray56", face= "bold"),
        axis.text.y=element_text(size=16, colour = "gray56", face = "bold"),
        axis.title=element_text(size=16,face="bold", colour = "gray56"),
        plot.title = element_text(size=18, face="bold", colour = "gray56"),
        legend.title = element_text(size=16, colour = "gray56", face= "bold"),
        legend.key = element_rect(fill = "transparent"),
        legend.text = element_text(size=16, colour = "gray56", face= "bold"),
        legend.background = element_rect(fill = "transparent", colour = NA),
        axis.line = element_line(size = 1.2, colour = "gray56"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA, linewidth = 1.8, color= "gray56"), 
        plot.background = element_rect(fill = "transparent", colour = NA)
  )

ggsave("G:/RESIST/Graciela/Presentations/CRC_poster/Tolalg.png", emmip_tolalg, bg="transparent")
