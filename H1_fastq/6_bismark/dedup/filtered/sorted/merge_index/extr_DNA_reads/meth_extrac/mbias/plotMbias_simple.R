#!/usr/bin/env Rscript

#setwd() where Mbias.txt files
#example input function:
#splitinfh("SRR10470169.DNA.M-bias.txt", "SRR10470169")

args <- commandArgs()
infh <- args[1]
sampleid <- args[2]
outdir <- args[3]

splitinfh <- function(infh, sampleid, outdir) {
  fh <- readLines(infh)
  totalen <- length(fh)
  fhsplits <- which(!nzchar(fh))
  nsplits <- length(fhsplits)
  #print(totalen)
  #print(fhsplits)
  #print(nsplits)
  
  #prepare dataframe
  dt <- data.frame(position=numeric(), cntmethylated=numeric(), cntunmethylated=numeric(), permethylation=numeric(), coverage=numeric(), context=character())
  #print(dt)
  
  #extract subheaders
  for (z in 1:nsplits) {
    subfh <- NULL
    if (z == 1) {
      startline <- 1
    } else {
      startline <- fhsplits[z-1]+1
    }
    endline <- fhsplits[z]-1
    
    subheader <- fh[startline]
    #print(subheader)
    
    #extract table from Mbias text file without top 2 lines (eg CpG context ===========) and change colnames
    subfh <- read.csv(text=fh[(startline+2):endline], sep="\t", header=TRUE)
    colnames(subfh) <- c("position", "cntmethylated", "cntunmethylated", "permethylation", "coverage")
    
    if (grepl("CpG", subheader)) {
      subfh$context = "CpG"
    } else if (grepl("CHG", subheader)) {
      subfh$context = "CHG"
    } else if (grepl("CHH", subheader)) {
      subfh$context = "CHH"
    }
    #put subfh for each context in premade, empty dt dataframe
    dt <- rbind(dt, subfh)

    outfig <- paste0(sampleid, ".plot.pdf")
    pdf(outdir/outfig, width=5.8, height=4)
    p1 <- ggplot(dt, aes(x=position, y=permethylation, color=context))+geom_line(aes(color=context))+theme_bw()+theme(plot.title=element_text(size=11),axis.text.x=element_text(color="black",size=11), axis.text.y=element_text(color="black",size=11),panel.grid.minor=element_blank())+scale_color_manual(values=c("red", "blue", "green"))+labs(x="R1 Base position", y="Percentage of methylation", title=sampleid)+scale_y_continuous(limits=c(0, 100), breaks=seq(0, 100, by = 10))
    print(p1)
    dev.off()
    
  } 
}



    
