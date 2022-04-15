library(haven)
library(dplyr)
library(tidyr)
library(gt)
library(Hmisc)

the_date <- as.character(Sys.Date())

# read the data
adsl <- read_sas("C:\\efficacy\\adsl.sas7bdat")
adrs <- read_sas("C:\\efficacy\\adrs.sas7bdat")

adsl$TRT01P <- gsub(" ", "", adsl$TRT01P)

adrs %>% 
  filter(RSEVAL=='Independent Central Review' & PARAMCD=='CBRSP')

adrs$TRT01P <- gsub(" ", "", adrs$TRT01P)

# get the big N in column headers from adsl
bign <- table(group=adsl$TRT01P)

# do the counting by TRT01P, AVALC

count0 <- 
  adrs %>%                   
  group_by(TRT01P, AVALC) %>%
  summarise(unique_subj = n_distinct(USUBJID))

# generate a frame data called comb
# it has all the treatment group and all response values CR, PR, SD, PD, NE

mat <- matrix(NA, nrow = 5, ncol = 3)
xy <- data.frame(mat)


xy[[1,1]] <- data.frame(X1='CR')
xy[[2,1]] <- data.frame(X1='PR')
xy[[3,1]] <- data.frame(X1='SD')
xy[[4,1]] <- data.frame(X1='PD')
xy[[5,1]] <- data.frame(X1='NE')

for (i in seq(1,5)){
  xy[[i,2]] <- data.frame(X2=i)
}

xy[[1,3]] <- data.frame(X3='Complete Response (CR)')
xy[[2,3]] <- data.frame(X3='Partial Response (PR)')
xy[[3,3]] <- data.frame(X3='Stable Disease (SD)')
xy[[4,3]] <- data.frame(X3='Progressive Disease (PD)')
xy[[5,3]] <- data.frame(X3='Not Evaluable (NE)')

X4 <- unique(count0$TRT01P)
comb <- xy %>%
  crossing(X4)

comb <- comb %>%
  mutate(AVALC=unlist(X1), TRT01P=X4, X2=unlist(X2), X3=unlist(X3))

# merge comb with count0 and calculate percentage 

m_count0 <- count0 %>%
  merge( comb, by=c("TRT01P", "AVALC"), all=TRUE) %>%
  mutate(denom=ifelse(TRT01P=='GroupA', bign[1], ifelse(TRT01P=='GroupB', bign[2], bign[3]))) %>%
  mutate(value=ifelse(is.na(unique_subj),"0", paste(unique_subj, "(", format(round(100*unique_subj/denom, 1), nsmall = 1), ")")))


# do the transpose 

a1 <- m_count0 %>%
 pivot_wider(id_cols=c(X3, X2, AVALC), names_from = TRT01P, values_from = value,
                  names_prefix = "") %>%
  mutate(block="Confirmed Best Overall Response", catlabel=X3)

a1 <- a1[order(a1$X2), ]


# do the ORR

mat <- matrix(NA, nrow = 1, ncol = 3)
z <- data.frame(mat)
z$GroupA <- z$X1
z$GroupB <- z$X2
z$GroupC <- z$X3

x <- table(group=adsl$TRT01P)
n <- table(group=adrs[(adrs$AVALC=='CR' | adrs$AVALC=='PR'),]$TRT01P)

# create a function to combine the proportion and confidence interval 
getci <- function (grp){
  grp <- {{grp}}
  n[grp] <- ifelse(is.na(n[grp]), 0, n[grp])
  t <- binconf(n[grp], x[grp], method="exact")
  ci <- paste0(round(100*t[1], digits=1), ' (', round(100*t[2], digits=1), ', ', round(100*t[3], digits=1),')')
  return(ci)
}

z$GroupA <- getci("GroupA")
z$GroupB <- getci("GroupB")
z$GroupC <- getci("GroupC")
z$catlabel <- "ORR (95% CI) [a]"
z$block <- "Objective Response Rate"

# do the CBR

mat <- matrix(NA, nrow = 1, ncol = 3)
y <- data.frame(mat)
y$GroupA <- y$X1
y$GroupB <- y$X2
y$GroupC <- y$X3

x <- table(group=adsl$TRT01P)
n <- table(group=adrs[(adrs$AVALC=='CR' | adrs$AVALC=='PR' | adrs$AVALC=="SD"),]$TRT01P)

# create a function to combine the proportion and confidence interval 
getci <- function (grp){
  grp <- {{grp}}
  n[grp] <- ifelse(is.na(n[grp]), 0, n[grp])
  t <- binconf(n[grp], x[grp], method="exact")
  ci <- paste0(round(100*t[1], digits=1), ' (', round(100*t[2], digits=1), ', ', round(100*t[3], digits=1),')')
  return(ci)
}

y$GroupA <- getci("GroupA")
y$GroupB <- getci("GroupB")
y$GroupC <- getci("GroupC")
y$catlabel <- "CBR (95% CI) [a]"
y$block <- "Clinical Benefit Rate"


a1 <- a1 %>%
  select(block,catlabel,GroupA, GroupB, GroupC)
z <- z %>%
  select(block,catlabel,GroupA, GroupB, GroupC)
y <- y %>%
  select(block,catlabel,GroupA, GroupB, GroupC)
  
df <- rbind(a1, z, y)

df %>% 
  gt(groupname_col="block") 


# use gt to do the reporting 
tab_html <- df %>% 
  gt(groupname_col="block") %>%
  
  tab_header(
    title = "Table 14.2.1 Confirmed Best Overall Response based on BICR assessment",
    subtitle = "ITT Population"
  ) %>%
  
  tab_source_note(
    source_note = "[a]: the confidence interval is based on Clopper-Pearson method."
  ) %>%
  
  
  tab_source_note(
    source_note = paste('Program Source: efficacy2.R            Executed: (Draft)',  the_date)
  ) %>%
  
  cols_label(
    catlabel= " ", 
    GroupA = paste0("Group A (N=", bign[1], ")"),  
    GroupB = paste0("Group B (N=", bign[2], ")"), 
    GroupC = paste0("Group C (N=", bign[3], ")")
    
  ) %>%
  
  tab_options(
    table.border.top.color = "white",
    heading.border.bottom.color = "black",
    table.border.bottom.color = "white",
    table_body.border.bottom.color = "black",
    table_body.hlines.color = "white", 
    row_group.border.bottom.color = "white", 
    row_group.border.top.color = "white", 
    column_labels.border.top.color = "black",
    column_labels.border.bottom.color = "black",
  ) %>%
  
  cols_align(
    align = "left",
    columns = c(catlabel)
  ) 

# output the HTML table

tab_html %>%
  gtsave("efficacy2.html", path = "C:\\efficacy" )
















