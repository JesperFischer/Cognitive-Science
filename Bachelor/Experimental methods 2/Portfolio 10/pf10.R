a = read.csv("emp_all_all.csv")

library(ggpubr)


colnames(a) <- c(1:73)
rownames(a) <- c(1:104)

parallel = fa.parallel(gammel, fm = "minres", fa = "fa", sim = FALSE)
parallel = fa.parallel(a, fm = "minres", fa = "fa", sim = FALSE)


#when our boss asks us what we think of his questionare, we might first check the data 
#with and without his test wether we get a different amount of factors. 

gammel = a[,1:58]

parallel = fa.parallel(gammel, fm = "minres", fa = "fa", sim = FALSE)
parallel = fa.parallel(a, fm = "minres", fa = "fa", sim = FALSE)

#as can be seen it seems like there are 4 factors if we run the analysis with his data
#or without indicating that it does not measure anything new. However what if we compared
#a model with 4 factors to one of 5. Because the scree plots seem to be around 4-5 factors.

#boring code:
fagamel = fa(r = gammel, nfactors = 4, rotate = "oblimin", fm = "minres")
fagamel5 = fa(r = gammel, nfactors = 5, rotate = "oblimin", fm = "minres")
F1 = fa(r = a, nfactors = 4, rotate = "oblimin", fm = "minres")
F2 = fa(r = a, nfactors = 5, rotate = "oblimin", fm = "minres")

print(F3, cut = 0.3)


load = as.data.frame(F1$loadings[,1:4])
load1 = as.data.frame(fagamel$loadings[,1:4])
load11 = as.data.frame(F2$loadings[,1:5])
load111 = as.data.frame(fagamel5$loadings[,1:5])

load$question_number = seq(73)
load1$question_number = seq(58)
load11$question_number = seq(73)
load111$question_number = seq(58)

d_MR1111 = load111 %>% select(c(1,6))
d_MR2111 = load111 %>% select(c(4,6))
d_MR5111 = load111 %>% select(c(3,6))
d_MR4111 = load111 %>% select(c(2,6))
d_MR3111 = load111 %>% select(c(5,6))


d_MR111 = load11 %>% select(c(1,6))
d_MR211 = load11 %>% select(c(3,6))
d_MR511 = load11 %>% select(c(4,6))
d_MR411 = load11 %>% select(c(2,6))
d_MR311 = load11 %>% select(c(5,6))

d_MR11 = load1 %>% select(c(1,5))
d_MR21 = load1 %>% select(c(3,5))
d_MR31 = load1 %>% select(c(4,5))
d_MR41 = load1 %>% select(c(2,5))


d_MR1 = load %>% select(c(1,5))
d_MR2 = load %>% select(c(3,5))
d_MR3 = load %>% select(c(4,5))
d_MR4 = load %>% select(c(2,5))


d_MR11$col = d_MR11$question_number
d_MR21$col = d_MR21$question_number
d_MR31$col = d_MR31$question_number
d_MR41$col = d_MR41$question_number

d_MR1$col = d_MR1$question_number
d_MR2$col = d_MR2$question_number
d_MR3$col = d_MR3$question_number
d_MR4$col = d_MR4$question_number

d_MR111$col = d_MR111$question_number
d_MR211$col = d_MR211$question_number
d_MR311$col = d_MR311$question_number
d_MR411$col = d_MR411$question_number
d_MR511$col = d_MR511$question_number

d_MR1111$col = d_MR1111$question_number
d_MR2111$col = d_MR2111$question_number
d_MR3111$col = d_MR3111$question_number
d_MR4111$col = d_MR4111$question_number
d_MR5111$col = d_MR5111$question_number

d_MR1$col = ifelse(d_MR1$col < 74, "PEST", d_MR1$col)
d_MR1$col = ifelse(d_MR1$question_number < 59, "Distress", d_MR1$col)
d_MR1$col = ifelse(d_MR1$question_number < 52, "Perspective Taking", d_MR1$col)
d_MR1$col = ifelse(d_MR1$question_number < 45, "Concern", d_MR1$col)
d_MR1$col = ifelse(d_MR1$question_number < 38, "Fantasy", d_MR1$col)
d_MR1$col = ifelse(d_MR1$question_number < 31, "BEES", d_MR1$col)

d_MR11$col = ifelse(d_MR11$col < 74, "PEST", d_MR11$col)
d_MR11$col = ifelse(d_MR11$question_number < 59, "Distress", d_MR11$col)
d_MR11$col = ifelse(d_MR11$question_number < 52, "Perspective Taking", d_MR11$col)
d_MR11$col = ifelse(d_MR11$question_number < 45, "Concern", d_MR11$col)
d_MR11$col = ifelse(d_MR11$question_number < 38, "Fantasy", d_MR11$col)
d_MR11$col = ifelse(d_MR11$question_number < 31, "BEES", d_MR11$col)

d_MR111$col = ifelse(d_MR111$col < 74, "PEST", d_MR111$col)
d_MR111$col = ifelse(d_MR111$question_number < 59, "Distress", d_MR111$col)
d_MR111$col = ifelse(d_MR111$question_number < 52, "Perspective Taking", d_MR111$col)
d_MR111$col = ifelse(d_MR111$question_number < 45, "Concern", d_MR111$col)
d_MR111$col = ifelse(d_MR111$question_number < 38, "Fantasy", d_MR111$col)
d_MR111$col = ifelse(d_MR111$question_number < 31, "BEES", d_MR111$col)

d_MR1111$col = ifelse(d_MR1111$col < 74, "PEST", d_MR1111$col)
d_MR1111$col = ifelse(d_MR1111$question_number < 59, "Distress", d_MR1111$col)
d_MR1111$col = ifelse(d_MR1111$question_number < 52, "Perspective Taking", d_MR1111$col)
d_MR1111$col = ifelse(d_MR1111$question_number < 45, "Concern", d_MR1111$col)
d_MR1111$col = ifelse(d_MR1111$question_number < 38, "Fantasy", d_MR1111$col)
d_MR1111$col = ifelse(d_MR1111$question_number < 31, "BEES", d_MR1111$col)


d_MR2$col = ifelse(d_MR2$col < 74, "PEST", d_MR2$col)
d_MR2$col = ifelse(d_MR2$question_number < 59, "Distress", d_MR2$col)
d_MR2$col = ifelse(d_MR2$question_number < 52, "Perspective Taking", d_MR2$col)
d_MR2$col = ifelse(d_MR2$question_number < 45, "Concern", d_MR2$col)
d_MR2$col = ifelse(d_MR2$question_number < 38, "Fantasy", d_MR2$col)
d_MR2$col = ifelse(d_MR2$question_number < 31, "BEES", d_MR2$col)

d_MR21$col = ifelse(d_MR21$col < 74, "PEST", d_MR21$col)
d_MR21$col = ifelse(d_MR21$question_number < 59, "Distress", d_MR21$col)
d_MR21$col = ifelse(d_MR21$question_number < 52, "Perspective Taking", d_MR21$col)
d_MR21$col = ifelse(d_MR21$question_number < 45, "Concern", d_MR21$col)
d_MR21$col = ifelse(d_MR21$question_number < 38, "Fantasy", d_MR21$col)
d_MR21$col = ifelse(d_MR21$question_number < 31, "BEES", d_MR21$col)

d_MR211$col = ifelse(d_MR211$col < 74, "PEST", d_MR211$col)
d_MR211$col = ifelse(d_MR211$question_number < 59, "Distress", d_MR211$col)
d_MR211$col = ifelse(d_MR211$question_number < 52, "Perspective Taking", d_MR211$col)
d_MR211$col = ifelse(d_MR211$question_number < 45, "Concern", d_MR211$col)
d_MR211$col = ifelse(d_MR211$question_number < 38, "Fantasy", d_MR211$col)
d_MR211$col = ifelse(d_MR211$question_number < 31, "BEES", d_MR211$col)

d_MR2111$col = ifelse(d_MR2111$col < 74, "PEST", d_MR2111$col)
d_MR2111$col = ifelse(d_MR2111$question_number < 59, "Distress", d_MR2111$col)
d_MR2111$col = ifelse(d_MR2111$question_number < 52, "Perspective Taking", d_MR2111$col)
d_MR2111$col = ifelse(d_MR2111$question_number < 45, "Concern", d_MR2111$col)
d_MR2111$col = ifelse(d_MR2111$question_number < 38, "Fantasy", d_MR2111$col)
d_MR2111$col = ifelse(d_MR2111$question_number < 31, "BEES", d_MR2111$col)



d_MR3$col = ifelse(d_MR3$col < 74, "PEST", d_MR3$col)
d_MR3$col = ifelse(d_MR3$question_number < 59, "Distress", d_MR3$col)
d_MR3$col = ifelse(d_MR3$question_number < 52, "Perspective Taking", d_MR3$col)
d_MR3$col = ifelse(d_MR3$question_number < 45, "Concern", d_MR3$col)
d_MR3$col = ifelse(d_MR3$question_number < 38, "Fantasy", d_MR3$col)
d_MR3$col = ifelse(d_MR3$question_number < 31, "BEES", d_MR3$col)

d_MR31$col = ifelse(d_MR31$col < 74, "PEST", d_MR31$col)
d_MR31$col = ifelse(d_MR31$question_number < 59, "Distress", d_MR31$col)
d_MR31$col = ifelse(d_MR31$question_number < 52, "Perspective Taking", d_MR31$col)
d_MR31$col = ifelse(d_MR31$question_number < 45, "Concern", d_MR31$col)
d_MR31$col = ifelse(d_MR31$question_number < 38, "Fantasy", d_MR31$col)
d_MR31$col = ifelse(d_MR31$question_number < 31, "BEES", d_MR31$col)

d_MR311$col = ifelse(d_MR311$col < 74, "PEST", d_MR311$col)
d_MR311$col = ifelse(d_MR311$question_number < 59, "Distress", d_MR311$col)
d_MR311$col = ifelse(d_MR311$question_number < 52, "Perspective Taking", d_MR311$col)
d_MR311$col = ifelse(d_MR311$question_number < 45, "Concern", d_MR311$col)
d_MR311$col = ifelse(d_MR311$question_number < 38, "Fantasy", d_MR311$col)
d_MR311$col = ifelse(d_MR311$question_number < 31, "BEES", d_MR311$col)

d_MR3111$col = ifelse(d_MR3111$col < 74, "PEST", d_MR3111$col)
d_MR3111$col = ifelse(d_MR3111$question_number < 59, "Distress", d_MR3111$col)
d_MR3111$col = ifelse(d_MR3111$question_number < 52, "Perspective Taking", d_MR3111$col)
d_MR3111$col = ifelse(d_MR3111$question_number < 45, "Concern", d_MR3111$col)
d_MR3111$col = ifelse(d_MR3111$question_number < 38, "Fantasy", d_MR3111$col)
d_MR3111$col = ifelse(d_MR3111$question_number < 31, "BEES", d_MR3111$col)


d_MR4$col = ifelse(d_MR4$col < 74, "PEST", d_MR4$col)
d_MR4$col = ifelse(d_MR4$question_number < 59, "Distress", d_MR4$col)
d_MR4$col = ifelse(d_MR4$question_number < 52, "Perspective Taking", d_MR4$col)
d_MR4$col = ifelse(d_MR4$question_number < 45, "Concern", d_MR4$col)
d_MR4$col = ifelse(d_MR4$question_number < 38, "Fantasy", d_MR4$col)
d_MR4$col = ifelse(d_MR4$question_number < 31, "BEES", d_MR4$col)


d_MR41$col = ifelse(d_MR41$col < 74, "PEST", d_MR41$col)
d_MR41$col = ifelse(d_MR41$question_number < 59, "Distress", d_MR41$col)
d_MR41$col = ifelse(d_MR41$question_number < 52, "Perspective Taking", d_MR41$col)
d_MR41$col = ifelse(d_MR41$question_number < 45, "Concern", d_MR41$col)
d_MR41$col = ifelse(d_MR41$question_number < 38, "Fantasy", d_MR41$col)
d_MR41$col = ifelse(d_MR41$question_number < 31, "BEES", d_MR41$col)

d_MR411$col = ifelse(d_MR411$col < 74, "PEST", d_MR411$col)
d_MR411$col = ifelse(d_MR411$question_number < 59, "Distress", d_MR411$col)
d_MR411$col = ifelse(d_MR411$question_number < 52, "Perspective Taking", d_MR411$col)
d_MR411$col = ifelse(d_MR411$question_number < 45, "Concern", d_MR411$col)
d_MR411$col = ifelse(d_MR411$question_number < 38, "Fantasy", d_MR411$col)
d_MR411$col = ifelse(d_MR411$question_number < 31, "BEES", d_MR411$col)

d_MR4111$col = ifelse(d_MR4111$col < 74, "PEST", d_MR4111$col)
d_MR4111$col = ifelse(d_MR4111$question_number < 59, "Distress", d_MR4111$col)
d_MR4111$col = ifelse(d_MR4111$question_number < 52, "Perspective Taking", d_MR411$col)
d_MR4111$col = ifelse(d_MR4111$question_number < 45, "Concern", d_MR4111$col)
d_MR4111$col = ifelse(d_MR4111$question_number < 38, "Fantasy", d_MR4111$col)
d_MR4111$col = ifelse(d_MR4111$question_number < 31, "BEES", d_MR4111$col)


d_MR511$col = ifelse(d_MR511$col < 74, "PEST", d_MR511$col)
d_MR511$col = ifelse(d_MR511$question_number < 59, "Distress", d_MR511$col)
d_MR511$col = ifelse(d_MR511$question_number < 52, "Perspective Taking", d_MR511$col)
d_MR511$col = ifelse(d_MR511$question_number < 45, "Concern", d_MR511$col)
d_MR511$col = ifelse(d_MR511$question_number < 38, "Fantasy", d_MR511$col)
d_MR511$col = ifelse(d_MR511$question_number < 31, "BEES", d_MR511$col)

d_MR5111$col = ifelse(d_MR5111$col < 74, "PEST", d_MR5111$col)
d_MR5111$col = ifelse(d_MR5111$question_number < 59, "Distress", d_MR5111$col)
d_MR5111$col = ifelse(d_MR5111$question_number < 52, "Perspective Taking", d_MR5111$col)
d_MR5111$col = ifelse(d_MR5111$question_number < 45, "Concern", d_MR5111$col)
d_MR5111$col = ifelse(d_MR5111$question_number < 38, "Fantasy", d_MR5111$col)
d_MR5111$col = ifelse(d_MR5111$question_number < 31, "BEES", d_MR5111$col)

#5 factors old
d_MR1111 <- d_MR1111 %>% arrange(MR1) %>% mutate(question_number = factor(question_number, levels = question_number))
ggplot(d_MR1111, aes(question_number, MR1, fill = d_MR1111$col, colour = d_MR1111$col)) + 
  geom_col(width = 0.6) + ylim(-2,2)+ggtitle("5 factors without PEST")+xlab("Question number")+theme(plot.title = element_text(hjust = 0.5))

ggsave("5 factors without PEST MR1.png", width = 15)

d_MR2111 <- d_MR2111 %>% arrange(MR2) %>% mutate(question_number = factor(question_number, levels = question_number))
ggplot(d_MR2111, aes(question_number, MR2, fill = d_MR2111$col, colour = d_MR2111$col)) + 
  geom_col(width = 0.6) + ylim(-2,2)+ggtitle("5 factors without PEST")+xlab("Question number")+theme(plot.title = element_text(hjust = 0.5))

ggsave("5 factors without PEST MR2.png", width = 15)


d_MR3111 <- d_MR3111 %>% arrange(MR3) %>% mutate(question_number = factor(question_number, levels = question_number))
ggplot(d_MR3111, aes(question_number, MR3, fill = d_MR3111$col, colour = d_MR3111$col)) + 
  geom_col(width = 0.6) + ylim(-2,2)+ggtitle("5 factors without PEST")+xlab("Question number")+theme(plot.title = element_text(hjust = 0.5))

ggsave("5 factors without PEST MR3.png", width = 15)


d_MR4111 <- d_MR4111 %>% arrange(MR4) %>% mutate(question_number = factor(question_number, levels = question_number))
ggplot(d_MR4111, aes(question_number, MR4, fill = d_MR4111$col, colour = d_MR4111$col)) + 
  geom_col(width = 0.6) + ylim(-2,2)+ggtitle("5 factors without PEST")+xlab("Question number")+theme(plot.title = element_text(hjust = 0.5))

ggsave("5 factors without PEST MR4.png", width = 15)


d_MR5111 <- d_MR5111 %>% arrange(MR5) %>% mutate(question_number = factor(question_number, levels = question_number))
ggplot(d_MR5111, aes(question_number, MR5, fill = d_MR5111$col, colour = d_MR5111$col)) + 
  geom_col(width = 0.6) + ylim(-2,2)+ggtitle("5 factors without PEST")+xlab("Question number")+theme(plot.title = element_text(hjust = 0.5))

ggsave("5 factors without PEST MR5.png", width = 15)


#5 factors new

d_MR111 <- d_MR111 %>% arrange(MR1) %>% mutate(question_number = factor(question_number, levels = question_number))
ggplot(d_MR111, aes(question_number, MR1, fill = d_MR111$col, colour = d_MR111$col)) + 
  geom_col(width = 0.6) + ylim(-2,2)+ggtitle("5 factors with PEST")+xlab("Question number")+theme(plot.title = element_text(hjust = 0.5))

ggsave("ggplot 5 factors with PEST MR1.png", width = 15)


d_MR211 <- d_MR211 %>% arrange(MR2) %>% mutate(question_number = factor(question_number, levels = question_number))
ggplot(d_MR211, aes(question_number, MR2, fill = d_MR211$col, colour = d_MR211$col)) + 
  geom_col(width = 0.5) + ylim(-2,2)+ggtitle("5 factors with PEST")+xlab("Question number")+
  rotate_x_text()+theme(plot.title = element_text(hjust = 0.5))

ggsave("ggplot 5 factors with PEST MR2.png", width = 15)


d_MR311 <- d_MR311 %>% arrange(MR3) %>% mutate(question_number = factor(question_number, levels = question_number))
ggplot(d_MR311, aes(question_number, MR3, fill = d_MR311$col, colour = d_MR311$col)) + 
  geom_col(width = 0.5) + ylim(-2,2)+ggtitle("5 factors with PEST")+xlab("Question number")+
  rotate_x_text()+theme(plot.title = element_text(hjust = 0.5))

ggsave("ggplot 5 factors with PEST MR3.png", width = 15)


d_MR411 <- d_MR411 %>% arrange(MR4) %>% mutate(question_number = factor(question_number, levels = question_number))
ggplot(d_MR411, aes(question_number, MR4, fill = d_MR411$col, colour = d_MR411$col)) + 
  geom_col(width = 0.5) + ylim(-2,2)+ggtitle("5 factors with PEST")+xlab("Question number")+
  rotate_x_text()+theme(plot.title = element_text(hjust = 0.5))

ggsave("ggplot 5 factors with PEST MR4.png", width = 15)


d_MR511 <- d_MR511 %>% arrange(MR5) %>% mutate(question_number = factor(question_number, levels = question_number))
ggplot(d_MR511, aes(question_number, MR5, fill = d_MR511$col, colour = d_MR511$col)) + 
  geom_col(width = 0.5) + ylim(-2,2)+ggtitle("5 factors with PEST")+xlab("Question number")+
  rotate_x_text()+theme(plot.title = element_text(hjust = 0.5))

ggsave("ggplot 5 factors with PEST MR5.png", width = 15)

#4 factors new

d_MR1 <- d_MR1 %>% arrange(MR1) %>% mutate(question_number = factor(question_number, levels = question_number))
ggplot(d_MR1, aes(question_number, MR1, fill = d_MR1$col, colour = d_MR1$col)) + 
  geom_col(width = 0.5) + ylim(-2,2)+ggtitle("4 factors with PEST")+xlab("Question number")+
  rotate_x_text()+theme(plot.title = element_text(hjust = 0.5))

ggsave("ggplot 4 factors with PEST MR1.png", width = 15)


d_MR2 <- d_MR2 %>% arrange(MR2) %>% mutate(question_number = factor(question_number, levels = question_number))
ggplot(d_MR2, aes(question_number, MR2, fill = d_MR2$col, colour = d_MR2$col)) + 
  geom_col(width = 0.5) + ylim(-2,2)+ggtitle("4 factors with PEST")+xlab("Question number")+
  rotate_x_text()+theme(plot.title = element_text(hjust = 0.5))

ggsave("ggplot 4 factors with PEST MR2.png", width = 15)


d_MR3 <- d_MR3 %>% arrange(MR3) %>% mutate(question_number = factor(question_number, levels = question_number))
ggplot(d_MR3, aes(question_number, MR3, fill = d_MR3$col, colour = d_MR3$col)) + 
  geom_col(width = 0.5) + ylim(-2,2)+ggtitle("4 factors with PEST")+xlab("Question number")+
  rotate_x_text()+theme(plot.title = element_text(hjust = 0.5))


ggsave("ggplot 4 factors with PEST MR3.png", width = 15)


d_MR4 <- d_MR4 %>% arrange(MR4) %>% mutate(question_number = factor(question_number, levels = question_number))
ggplot(d_MR4, aes(question_number, MR4, fill = d_MR4$col, colour = d_MR4$col)) + 
  geom_col(width = 0.5) + ylim(-2,2)+ggtitle("4 factors with PEST")+xlab("Question number")+
  rotate_x_text()+theme(plot.title = element_text(hjust = 0.5))


ggsave("ggplot 4 factors with PEST MR4.png", width = 15)




d_MR11 <- d_MR11 %>% arrange(MR3) %>% mutate(question_number = factor(question_number, levels = question_number))
ggplot(d_MR11, aes(question_number, MR3, fill = d_MR11$col, colour = d_MR11$col)) + 
  geom_col(width = 0.5) + ylim(-2,2)+ggtitle("4 factors without PEST")+xlab("Question number")+
  rotate_x_text()+theme(plot.title = element_text(hjust = 0.5))

ggsave("4 factors without PEST MR3.png", width = 15)


d_MR21 <- d_MR21 %>% arrange(MR1) %>% mutate(question_number = factor(question_number, levels = question_number))

ggplot(d_MR21, aes(question_number, MR1, fill = d_MR21$col, colour = d_MR21$col)) + 
  geom_col(width = 0.5) + ylim(-2,2)+ggtitle("4 factors without PEST")+xlab("Question number")+
  rotate_x_text()+theme(plot.title = element_text(hjust = 0.5))

ggsave("4 factors without PEST MR1.png", width = 15)


d_MR31 <- d_MR31 %>% arrange(MR2) %>% mutate(question_number = factor(question_number, levels = question_number))

ggplot(d_MR31, aes(question_number, MR2, fill = d_MR31$col, colour = d_MR31$col)) + 
  geom_col(width = 0.5) + ylim(-2,2)+ggtitle("4 factors without PEST")+xlab("Question number")+
  rotate_x_text()+theme(plot.title = element_text(hjust = 0.5))

ggsave("4 factors without PEST MR2.png", width = 15)


d_MR41 <- d_MR41 %>% arrange(MR4) %>% mutate(question_number = factor(question_number, levels = question_number))

ggplot(d_MR41, aes(question_number, MR4, fill = d_MR41$col, colour = d_MR41$col)) + 
  geom_col(width = 0.5) + ylim(-2,2)+ggtitle("4 factors without PEST")+xlab("Question number")+
  rotate_x_text()+theme(plot.title = element_text(hjust = 0.5))

ggsave("4 factors without PEST MR4.png", width = 15)



View(F1$scores)


p = as.data.frame(F1$scores)

result = which(p$MR1 == min(p$MR1), arr.ind = TRUE)

print(result)

library(tidyverse)


summarise_each(p, funs(max(., na.rm=TRUE)))




p1 = as.matrix(p)

heatmap(p1)


p$participant = seq(104)


install.packages("plotrix")
library(plotrix)






















print(result)

max(p$MR1)
max(p$MR2)
max(p$MR3)
max(p$MR4)
min(p$MR1)
min(p$MR2)
min(p$MR3)
min(p$MR4)












a2 = as.matrix(a)


ny = a[,59:73]


gammel = a[,1:58]



fagamel = fa.parallel(ny, fm = "minres", fa = "fa", sim = FALSE)

fany = fa.parallel(gammel, fm = "minres", fa = "fa", sim = FALSE)


fagamel = fa(r = gammel, nfactors = 4, rotate = "oblimin", fm = "minres")
fany = fa(r = ny, nfactors = 1, rotate = "oblimin", fm = "minres")



print(fagamel, cut = 0.4, order = TRUE)
print(fany, cut = 0.5, order = TRUE)

ss = fany$scores

library(psych)



gammelalt = fa33$scores

gammelalt = as.data.frame(gammelalt)


fa33 = fa(r = gammel, nfactors = 4, rotate = "oblimin", fm = "minres")

print(fa33, cut = 0.5, order = TRUE)





fa22 = fa(r = ny, nfactors = 1, rotate = "oblimin", fm = "minres")

d_MR1 = fa22$scores

d_MR1 = as.data.frame(d_MR1)

d_MR1$ID = seq(104)


d_MR1 <- d_MR1 %>% arrange(MR1) %>% mutate(ID = factor(ID, levels = ID))

ggplot(d_MR1, aes(ID, MR1, fill = MR1, color = MR1)) + 
  geom_col(stat = 'summary', fun.y = mean, width = 0.5) + ylim(-5,5)



max(d_MR1)


