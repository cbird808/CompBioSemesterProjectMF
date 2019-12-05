rm(list=ls())

AlbInFile <- "./AlbatrossBarcodes.csv"
ContInFile <- "./ContempBarcodes.csv"


library(tidyverse)
library(magrittr)
library(splines)
library(emmeans)


# Albatross %>% str
# 
# Contemp %>% str

Albatross <- read_csv(AlbInFile)
Contemp <- read_csv(ContInFile)

a<-1;b<-1

#take two input tibbles and create a tidy tibble that identifies errors
sequencing_error_data <- bind_rows(Albatross, Contemp, .id = 'Sampling') %>%
  mutate(Sampling = case_when(Sampling == '1' ~ 'Albatross',
                              Sampling == '2' ~ 'Contemp',
                              TRUE ~ 'error')) %>%
  mutate(Barcode = str_c('GG', Barcode, "GG", sep = '')) %>% #second pair is for the overhang of the cutsite
  gather(position, base, -Sampling:-Indels) %>%
  mutate(base_position = str_extract(position, "[0-9]+") %>% as.integer,
         true_base = str_sub(Barcode, base_position, base_position),
         true_base = if_else(true_base == "", NA_character_, true_base)) %>%
  group_by(Sampling, base_position, Barcode, Read, Indels) %>%
  summarise(correct = sum(base == true_base),
            error = sum(base != true_base),
            total = n()) %>%
  ungroup %>%
  mutate_at(vars(Read, Barcode, Sampling, base_position, Indels), as.factor)

error_model <- glm(cbind(correct, error) ~ Read*Barcode*Sampling*base_position, 
    data = sequencing_error_data, family = binomial)
summary(error_model)
anova(error_model, test = 'LRT')

error_emm <- emmeans(error_model, ~ Read*Barcode*Sampling*base_position)
pairwise_postHoc <- contrast(error_emm, method = 'consec', simple = 'each', combine = FALSE, adjust = 'holm')
pairwise_postHoc$`simple contrasts for base_position`
pairwise_postHoc$`simple contrasts for Sampling`

predict(error_model, se.fit = TRUE, type = 'response') %>%
  as_tibble() %>%
  select(-residual.scale) %>%
  bind_cols(sequencing_error_data, .) %>%
  
  #filter(base_position < 15) %>% #Temporary
  mutate(prop_correct = (correct + a)/(total + a + b),
            error_prop_lwr = qbeta(0.025, correct + a, total - correct + b),
            error_prop_upr = qbeta(0.975, correct + a, total - correct + b)) %>%
  ggplot(aes(x = base_position, ymin = fit - se.fit, ymax = fit + se.fit,
             colour = Sampling)) +
  #geom_point(aes(y = prop_correct), size = 5, shape = 'square') +
  geom_linerange(size = 1) +
  geom_point(aes(y = fit), size = 3, shape = 'circle') +
  facet_grid(Barcode~Read)



################################################################
sequencing_error_data %>% mutate(error_rate = error/total) -> sequencing_error_data
sequencing_error_data %>% mutate(success_rate = correct/total) -> sequencing_error_data
pdf(file="./SequencingErrorData1.pdf")
ggplot(data=sequencing_error_data, aes(x=base_position, y=success_rate, fill=Sampling)) + 
  geom_boxplot() +
  facet_grid(Barcode ~ Read)
dev.off()
pdf(file="./SequencingErrorData2.pdf")
ggplot(data=sequencing_error_data, aes(x=base_position, y=correct, fill=Sampling)) + 
  geom_boxplot() +
  facet_grid(Barcode ~ Read)
dev.off()
pdf(file="./SequencingErrorDataAll.pdf")
ggplot(data=sequencing_error_data, aes(x=base_position, y=error_rate, fill=Sampling)) + 
  geom_boxplot() +
  facet_grid(Barcode ~ Indels)
dev.off()

# Subset Two deletions sequences 
TwoDel <- sequencing_error_data[sequencing_error_data$Indels == "2del",]
pdf(file="./SequencingErrorData2del.pdf")
ggplot(data=TwoDel, aes(x=base_position, y=error_rate, fill=Sampling)) +
  ggtitle("2del") + 
  geom_boxplot() +
  facet_grid(Barcode ~ Read)
dev.off()

# Subset one deletion sequences 
OneDel <- sequencing_error_data[sequencing_error_data$Indels == "1del",]
pdf(file="./SequencingErrorData1del.pdf")
ggplot(data=OneDel, aes(x=base_position, y=error_rate, fill=Sampling)) +
  ggtitle("1del") + 
  geom_boxplot() +
  facet_grid(Barcode ~ Read)
dev.off()

# Subset No indel sequences
NoInd <- sequencing_error_data[sequencing_error_data$Indels == "0Ind",]
pdf(file="./SequencingErrorData0Ind.pdf")
ggplot(data=NoInd, aes(x=base_position, y=error_rate, fill=Sampling)) + 
  ggtitle("0Ind") +
  geom_boxplot() +
  facet_grid(Barcode ~ Read)
dev.off()

# Subset one insertion sequences
OneIns <- sequencing_error_data[sequencing_error_data$Indels == "1Ins",]
pdf(file="./SequencingErrorData1Ins.pdf")
ggplot(data=OneIns, aes(x=base_position, y=error_rate, fill=Sampling)) + 
  ggtitle("1Ins") +
  geom_boxplot() +
  facet_grid(Barcode ~ Read)
dev.off()

# Subset two insertion sequences
TwoIns <- sequencing_error_data[sequencing_error_data$Indels == "2Ins",]
pdf(file="./SequencingErrorData2Ins.pdf")
ggplot(data=TwoIns, aes(x=base_position, y=error_rate, fill=Sampling)) + 
  ggtitle("2Ins") +
  geom_boxplot() +
  facet_grid(Barcode ~ Read)
dev.off()

# Subset three insertion sequences
ThreeIns <- sequencing_error_data[sequencing_error_data$Indels == "3Ins",]
pdf(file="./SequencingErrorData3Ins.pdf")
ggplot(data=ThreeIns, aes(x=base_position, y=error_rate, fill=Sampling)) + 
  ggtitle("3Ins") +
  geom_boxplot() +
  facet_grid(Barcode ~ Read)
dev.off()

# Subset four insertion sequences
FourIns <- sequencing_error_data[sequencing_error_data$Indels == "4Ins",]
pdf(file="./SequencingErrorData4Ins.pdf")
ggplot(data=FourIns, aes(x=base_position, y=error_rate, fill=Sampling)) + 
  ggtitle("4Ins") +
  geom_boxplot() +
  facet_grid(Barcode ~ Read)
dev.off()

# Subset five insertion sequences
FiveIns <- sequencing_error_data[sequencing_error_data$Indels == "5Ins",]
pdf(file="./SequencingErrorData5Ins.pdf")
ggplot(data=FiveIns, aes(x=base_position, y=error_rate, fill=Sampling)) + 
  ggtitle("5Ins") +
  geom_boxplot() +
  facet_grid(Barcode ~ Read)
dev.off()

# Subset six insertion sequences
SixIns <- sequencing_error_data[sequencing_error_data$Indels == "6Ins",]
pdf(file="./SequencingErrorData6Ins.pdf")
ggplot(data=SixIns, aes(x=base_position, y=error_rate, fill=Sampling)) + 
  ggtitle("6Ins") +
  geom_boxplot() +
  facet_grid(Barcode ~ Read)
dev.off()

# Subset seven insertion sequences
SevenIns <- sequencing_error_data[sequencing_error_data$Indels == "7Ins",]
pdf(file="./SequencingErrorData7Ins.pdf")
ggplot(data=SevenIns, aes(x=base_position, y=error_rate, fill=Sampling)) + 
  ggtitle("7Ins") +
  geom_boxplot() +
  facet_grid(Barcode ~ Read)
dev.off()

# Subset eight insertion sequences
EightIns <- sequencing_error_data[sequencing_error_data$Indels == "8Ins",]
pdf(file="./SequencingErrorData8Ins.pdf")
ggplot(data=EightIns, aes(x=base_position, y=error_rate, fill=Sampling)) + 
  ggtitle("8Ins") +
  geom_boxplot() +
  facet_grid(Barcode ~ Read)
dev.off()

#######################################################################################################################


