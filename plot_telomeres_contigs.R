library(circlize)
rm(list=ls())
n = 10

yahsScaffolds = read.table("yahs.out_scaffolds_final.bed")

scaffolds = unique(yahsScaffolds[,1])
scaffolds

if(n > 0){
    scaffolds = scaffolds[1:n]
}

scaffolds

coords = as.data.frame(matrix(unlist(lapply(scaffolds, function(x){
    c(min(yahsScaffolds[ yahsScaffolds[,1] == x, 2:3]),
      max(yahsScaffolds[ yahsScaffolds[,1] == x, 2:3]))
})), ncol=2, byrow=TRUE))

colnames(coords) = c("start", "end")
coords$start

df = data.frame(name=scaffolds,
                start = coords$start,
                end=coords$end)
circos.genomicInitialize(df)
scafs1.10 = yahsScaffolds[ yahsScaffolds[,1] %in% scaffolds, ]
cols=sample(colors(), 10)
tels.data = read.table("telomeres/tidk_w_1000_telomeric_repeat_windows.1.bed")
colnames(tels.data) = c("scaffold", "start", "end", "n")
tels.data$n = tels.data$n/(max(tels.data$n))
tp_family = scafs1.10
names(tp_family) = c("scaffold", "start", "end", "contig")

conts = c()
i=1
for(i in 1:nrow(tels.data)){
    p = (tels.data$start[i] + tels.data$end[i])/2
    inds = tp_family[,1] == tels.data[i,1]
    print( c(p, tp_family[inds,"start"], tp_family[inds, "end"]) )
    truInds = p+500 >= tp_family[ inds, "start"] & p-500 <= tp_family[inds, "end"]
    print(c(sum(inds), sum(truInds)))
    conts = c(conts, (tp_family$contig[inds][which(truInds)[1]]))
}

conts



tels.data$c = conts

circos.genomicInitialize(scafs1.10)
circos.track(ylim = c(0, 1), 
    bg.col = cols,
    bg.border = NA, track.height = 0.05)

fdata = data.frame(scaffold=c(tp_family$scaffold, tels.data$scaffold),
                   start=c(tp_family$start, tels.data$start),
                   end=c(tp_family$end, tels.data$end),
                   contig=c(tp_family$contig, rep(NA, nrow(tels.data))),
                   n=c(rep(NA, nrow(tp_family)), tels.data$n),
                   col=c(rep(NA, nrow(tp_family)), as.factor(tels.data$c)))

n = max(tapply(tp_family$contig, tp_family$scaffold, function(x) length(unique(x))))
circos.genomicTrack(fdata, ylim = c(0.3, 1+0.4), 
                    panel.fun = function(region, value, ...) {
                        print(region)
                        print(value)
                        all_tx = unique(value$contig)
                        print(all_tx)
                        for(i in seq_along(all_tx)) {
                            l = value$contig == all_tx[i]
                                        # for each transcript
                            current_tx_start = min(region[l, 1])
                            current_tx_end = max(region[l, 2])
                            ##circos.lines(c(current_tx_start, current_tx_end), c(1,1),
                                         #c(n - i + 1, n - i + 1),
                            ##           col = "#CCCCCC")
                            circos.genomicRect(region[l, , drop = FALSE], ytop = 1 + 0.4, 
                                               ybottom = 1 - 0.7, col = "white", border = 'black')
                        }
                        allvals = which(!is.na(value$n))
                        circos.genomicPoints(region[allvals, ,drop=FALSE], value$n[allvals]+0.3, col=value$col[allvals], pch=value$col[allvals], cex=1.5, lw=2)
                        for(v in allvals){
                            print(v)
                        }
                    }, bg.border = "black", track.height = 0.3)



###########fd
