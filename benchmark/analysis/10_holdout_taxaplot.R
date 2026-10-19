#!/usr/bin/env Rscript
# =====================================================================
# Supplementary Figure S2: 16S holdout composition vs truth, genus AND species.
# The holdout HAS a truth (the GTDB reference lineage of each test sequence), so a
# composition taxaplot is meaningful. Each test sequence counts once; bars are the
# fraction of test sequences assigned to each taxon by Truth, BLCA, BLCA@80 and NB.
# Top-12 taxa by truth + Other + Unassigned. 16S only (the ITS panel was removed on
# revision; the holdout composition is shown for 16S at both genus and species rank).
# =====================================================================
.SCRIPT_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else normalizePath(getwd())
})
source(file.path(.SCRIPT_DIR, "00_common.R"))
comp1 <- function(taxa) {                 # composition (count=1 per seq), NA -> Unassigned
  t <- clean(taxa); t[is.na(t)] <- "Unassigned"
  as.data.frame(table(t)) %>% transmute(taxon=as.character(t), relab=Freq/sum(Freq))
}
holdout_comp <- function(marker, rank, ntop=12) {
  hd<-file.path(BR,sprintf("benchmark_holdout_%s_seed42",marker)); rd<-file.path(hd,"results")
  truth<-read_lineage_tsv(file.path(hd,"test_10.taxonomy"))
  bl<-read_blca(rd)$tab; b8<-read_blca_at80(rd); nb<-read_nb(rd)
  ids<-Reduce(intersect,list(truth$id,bl$id,b8$id,nb$id))
  ord<-function(x) x[match(ids,x$id),]
  T<-ord(truth); B<-ord(bl); B8<-ord(b8); N<-ord(nb)
  m<-list(Truth=comp1(T[[rank]]), BLCA=comp1(B[[rank]]), `BLCA@80`=comp1(B8[[rank]]), NB=comp1(N[[rank]]))
  all<-bind_rows(lapply(names(m),function(k) m[[k]] %>% mutate(method=k)))
  top<-m$Truth %>% filter(!taxon %in% c("Unassigned","Other")) %>% arrange(desc(relab)) %>% head(ntop) %>% pull(taxon)
  all<-all %>% mutate(grp=ifelse(taxon %in% top,taxon,ifelse(taxon=="Unassigned","Unassigned","Other"))) %>%
       group_by(method,grp) %>% summarise(relab=sum(relab),.groups="drop")
  lev<-c(top,"Other","Unassigned"); all$grp<-factor(all$grp,levels=lev)
  all$method<-factor(all$method,levels=c("Truth","BLCA","BLCA@80","NB"))
  cols<-setNames(c(TAXA_PAL[seq_along(top)],OTHER_COL,UNASSIGNED_COL),lev)
  list(df=all,cols=cols,marker=marker,rank=rank)
}
plot_h <- function(bc,title) ggplot(bc$df,aes(method,relab,fill=grp))+
  geom_col(width=0.72,colour="white",linewidth=0.2)+
  scale_fill_manual(values=bc$cols,name=NULL)+
  scale_y_continuous(labels=percent_format(accuracy=1),expand=expansion(mult=c(0,0.01)))+
  labs(x=NULL,y="Fraction of test sequences",title=title)+
  guides(fill=guide_legend(ncol=1,byrow=TRUE))+
  theme_ms(12)+theme(legend.position="right",legend.text=element_text(face="italic",size=8),
                     legend.key.size=unit(10,"pt"),axis.text.x=element_text(size=10))

g16<-holdout_comp("16s","Genus"); s16<-holdout_comp("16s","Species")
fig <- plot_h(g16,"a   16S holdout: genus") / plot_h(s16,"b   16S holdout: species")
ggsave(file.path(OUT,"figures","fig_holdout_taxaplot_16s.png"),fig,width=9,height=10,dpi=200)
ggsave(file.path(OUT,"figures","fig_holdout_taxaplot_16s.pdf"),fig,width=9,height=10)
cat("Wrote fig_holdout_taxaplot_16s.{png,pdf}\n")
