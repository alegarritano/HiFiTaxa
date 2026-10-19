#!/usr/bin/env Rscript
# =====================================================================
# Draft main Figure 2 (single results figure). 2 rows (16S, ITS) x 3 blocks:
#   mock detection (genus+species) | holdout per-rank accuracy | orphan backoff
# Emu/EMITS appear only on the mock (a holdout has no reads).
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
## ---------- mock detection data ----------
detc <- function(relab, exp_set, rank) {
  d <- relab %>% mutate(taxon=if(rank=="Species") norm_species(taxon) else norm_genus(taxon)) %>%
       filter(!is.na(taxon)) %>% group_by(taxon) %>% summarise(relab=sum(relab),.groups="drop")
  det <- length(intersect(exp_set, d$taxon[d$relab>0]))
  fo <- d %>% filter(!(taxon %in% exp_set), relab>FP_THRESH)
  fo <- if (rank=="Species") fo %>% filter(!is_placeholder_sp(taxon)) else
        fo %>% filter(!str_detect(taxon,"(ceae$| sp$|^Fungi|Incertae|_fam_|_gen_)"))
  tibble(rank=rank, Detected=det, Missed=length(exp_set)-det, `False pos.`=nrow(fo))
}
mock_detect <- function(marker) {
  if (marker=="16s"){dir<-file.path(BR,"results_atcc_16s");eg<-atcc_genus;es<-atcc_species;rl<-"Emu"}
  else              {dir<-file.path(BR,"results_fungi_its");eg<-fung_genus;es<-fung_species;rl<-"EMITS"}
  cnt<-read_feature_table(dir); bl<-read_blca(dir)$tab; b8<-read_blca_at80(dir); nb<-read_nb(dir)
  rc<-function(t,c) relab_from_counts(t$id,t[[c]],cnt)
  # share of reads (ASV read counts) with no assignment at this rank; read level profilers assign every read
  una<-function(t,c){ok<-t$id[!is.na(t[[c]]) & t[[c]]!=""]; 100*(1-sum(cnt$count[cnt$id %in% ok])/sum(cnt$count))}
  em<-function(rk) if(marker=="16s") read_emu(dir,tolower(rk)) else read_emits(dir,tolower(rk))
  bind_rows(
    detc(rc(bl,"Genus"),eg,"Genus")%>%mutate(cl="BLCA"), detc(rc(bl,"Species"),es,"Species")%>%mutate(cl="BLCA"),
    detc(rc(b8,"Genus"),eg,"Genus")%>%mutate(cl="BLCA@80"), detc(rc(b8,"Species"),es,"Species")%>%mutate(cl="BLCA@80"),
    detc(em("Genus"),eg,"Genus")%>%mutate(cl=rl), detc(em("Species"),es,"Species")%>%mutate(cl=rl),
    detc(rc(nb,"Genus"),eg,"Genus")%>%mutate(cl="NB"), detc(rc(nb,"Species"),es,"Species")%>%mutate(cl="NB")) %>%
    mutate(readlevel=rl) %>%
    mutate(unassigned=case_when(cl=="BLCA" & rank=="Genus"~una(bl,"Genus"), cl=="BLCA" & rank=="Species"~una(bl,"Species"),
                                cl=="BLCA@80" & rank=="Genus"~una(b8,"Genus"), cl=="BLCA@80" & rank=="Species"~una(b8,"Species"),
                                cl=="NB" & rank=="Genus"~una(nb,"Genus"), cl=="NB" & rank=="Species"~una(nb,"Species"),
                                TRUE~0))
}
acc <- read_csv(file.path(OUT,"tables","holdout_accuracy.csv"),show_col_types=FALSE,progress=FALSE)
bo  <- read_csv(file.path(OUT,"tables","holdout_orphan_backoff.csv"),show_col_types=FALSE,progress=FALSE)

## ---------- panel builders ----------
NUMSIZE <- 9/.pt       # bar-top number size (9 pt)
LEGSIZE <- 10
xt90 <- theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

p_detect <- function(dd, rl, leg, ceiling=NULL) {
  # dodged columns for Detected / Missed / False pos., coloured by classifier;
  # Expected shown as a grey reference column in the Detected group.
  nexp <- dd %>% group_by(rank) %>% summarise(nx=first(Detected+Missed),.groups="drop")
  expr <- nexp %>% transmute(rank, cl="Expected", Detected=nx, Missed=0, `False pos.`=0)
  d <- bind_rows(expr, dd %>% select(rank,cl,Detected,Missed,`False pos.`)) %>%
       pivot_longer(c(Detected,Missed,`False pos.`),names_to="cat",values_to="n")
  d$cat<-factor(d$cat,levels=c("Detected","Missed","False pos."))
  d$cl<-factor(d$cl,levels=c("Expected","BLCA","BLCA@80",rl,"NB")); d$rank<-factor(d$rank,levels=c("Genus","Species"))
  fmt_u<-function(x) ifelse(x<0.05 & x>0, "<0.1", formatC(x,format="f",digits=1))
  # unassigned reads (%) written inside the base of each Detected bar; Expected gets an empty label so the dodge aligns
  ul <- bind_rows(dd %>% transmute(rank, cl, lab=paste0(fmt_u(unassigned),"% unassigned")),
                  nexp %>% transmute(rank, cl="Expected", lab=""))
  ul$cat<-factor("Detected",levels=levels(d$cat)); ul$n<-0
  ul$cl<-factor(ul$cl,levels=levels(d$cl)); ul$rank<-factor(ul$rank,levels=c("Genus","Species"))
  g <- ggplot(d,aes(cat,n,fill=cl))+
    geom_col(position=position_dodge(0.8),width=0.8)+
    geom_text(data=ul,aes(x=cat,y=n,label=lab,group=cl),inherit.aes=FALSE,position=position_dodge(0.8),
              angle=90,hjust=-0.05,vjust=0.5,size=6.5/.pt,colour="white")+
    geom_text(aes(label=ifelse(n>0 | cat=="Detected",n,NA)),position=position_dodge(0.8),
              angle=90,hjust=-0.2,vjust=0.5,size=NUMSIZE)+
    facet_wrap(~rank,scales="free_y")+scale_fill_manual(values=CLF_COL,name=NULL)+
    scale_y_continuous(expand=expansion(mult=c(0,0.24)))+labs(x=NULL,y="# taxa")+
    theme_ms(11)+theme(legend.position=if(leg)"top" else "none",
                       legend.text=element_text(size=LEGSIZE))+xt90
  if (!is.null(ceiling)) {
    cd <- data.frame(rank=factor("Species",levels=c("Genus","Species")), y=ceiling)
    g <- g + geom_hline(data=cd,aes(yintercept=y),linetype="dashed",colour="grey30",linewidth=0.4) +
         geom_text(data=cd,aes(x=3.45,y=y,label=paste0(ceiling," nameable\nin UNITE")),inherit.aes=FALSE,
                   hjust=1,vjust=-0.25,size=8/.pt,colour="grey20",lineheight=0.9)
  }
  g
}
p_acc <- function(mk, leg) {
  # per-rank overall accuracy (% correct), dodged bars by classifier
  d<-acc%>%filter(marker==mk); d$Rank<-factor(d$Rank,levels=RANKS)
  if(mk=="its") levels(d$Rank)[levels(d$Rank)=="Domain"]<-"Kingdom"
  d$classifier<-factor(d$classifier,levels=c("BLCA","BLCA@80","NB"))
  ggplot(d,aes(Rank,accuracy,fill=classifier))+
    geom_col(position=position_dodge(0.8),width=0.75,colour="black",linewidth=0.15)+
    geom_text(aes(label=round(accuracy)),position=position_dodge(0.8),angle=90,hjust=-0.2,vjust=0.5,size=NUMSIZE)+
    scale_fill_manual(values=CLF_COL,name=NULL)+
    scale_y_continuous(limits=c(0,120),breaks=seq(0,100,25),oob=scales::oob_keep)+
    labs(x=NULL,y="Accuracy (% correct)")+
    theme_ms(11)+theme(legend.position=if(leg)"top" else "none",legend.text=element_text(size=LEGSIZE))+xt90
}
p_backoff <- function(mk, leg) {
  d<-bo%>%filter(marker==mk)%>%pivot_longer(c(Correct,Incorrect,Unassigned),names_to="cat",values_to="pct")
  d$cat<-factor(d$cat,levels=c("Unassigned","Incorrect","Correct"))
  d$classifier<-factor(d$classifier,levels=c("BLCA","BLCA@80","NB"))
  d$level<-factor(d$level,levels=c("Orphan species -> genus","Orphan genus -> family"),
                  labels=c("Species\n→ genus","Genus\n→ family"))
  ggplot(d,aes(classifier,pct,fill=cat))+geom_col(width=0.8)+facet_wrap(~level)+
    scale_fill_manual(values=CAT_COL,name=NULL,breaks=c("Correct","Incorrect","Unassigned"))+
    scale_y_continuous(expand=expansion(mult=c(0,0.02)))+labs(x=NULL,y="% of orphans")+
    theme_ms(11)+theme(legend.position=if(leg)"top" else "none",legend.text=element_text(size=LEGSIZE),
                       strip.text=element_text(size=8))+xt90
}

d16<-mock_detect("16s"); dIT<-mock_detect("its")
row16 <- p_detect(d16,"Emu",TRUE) | p_acc("16s",TRUE) | p_backoff("16s",TRUE)
rowIT <- p_detect(dIT,"EMITS",TRUE,ceiling=58) | p_acc("its",FALSE) | p_backoff("its",FALSE)
fig <- (row16 / rowIT) +
  patchwork::plot_layout(widths=c(1.3,1,1.15)) +
  patchwork::plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(face="bold", size=16))
ggsave(file.path(OUT,"figures","fig2_main.png"), fig, width=17, height=9, dpi=200)
ggsave(file.path(OUT,"figures","fig2_main.pdf"), fig, width=17, height=9)
write_csv(bind_rows(d16,dIT), file.path(OUT,"tables","mock_detection_unassigned.csv"))
cat("Wrote fig2_main.{png,pdf}\n")
