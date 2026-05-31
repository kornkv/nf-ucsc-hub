#!/usr/bin/env Rscript
suppressPackageStartupMessages({
library(data.table)
library(optparse)
})
##-------------------- read input ----------------
option_list <- list(
    make_option("--input", type="character", help="RepeatMasker .out file"),
    make_option("--outdir", type="character", default = "bed_by_repClass", help="Output directory [default: %default]"),
    make_option("--prefix", type="character", help="genome prefix")
    )
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$input), !is.null(opt$prefix))
rmsk  <- fread(
  cmd = sprintf(
    "awk 'BEGIN{OFS=\"\\t\"}{for(i=1;i<=15;i++)printf\"%%s%%s\",$i,(i==15?\"\\n\":OFS)}' %s",
    opt$input
  ),
  skip = 3,
  col.names = c(
    "SW_score","perc_div","perc_del","perc_ins",
    "seqname","start","end","remaining",
    "strand","repName","repClass","begin",
    "end_in_repeat","remaining_in_repeat","ID"
  )
)
# ----------------- Cleanup -----------------------
rmsk[strand == "C", strand := "-" ]
rmsk[,c("repClass_p1","repClass_p2"):=tstrsplit(repClass, "/")]
rmsk[repClass_p1 == "SINE?", repClass_p1 := "SINE"]
# ----------------- Convert to BED like ------------
bed <- rmsk[,.(
    chrom = seqname,
    start = as.integer(start - 1),
    end = as.integer(end),
    name = repName,
    score = 1000,
    strand = strand,
    thickStart = as.integer(start - 1),
    thickEnd = as.integer(end),
    itemRgb = "255,0,0",
    repClass_p1
)]
# ------------- priorities -----------------
preferred <- c("DNA","LINE","LTR","SINE","Unknown")
present  <- sort(unique(bed$repClass_p1))
preferred_present  <- preferred[preferred %in% present]
other_classes  <- setdiff(present, preferred_present)
priority_map  <-  c(
    setNames(seq_along(preferred_present), preferred_present),
    setNames(seq_along(other_classes) + length(preferred_present), other_classes)
)
# ------------ write BED and trackDb -----
dir.create(opt$outdir, showWarnings = FALSE)
trackdb_blocks <- list()
parent_block <-
"track RepeatMasker
compositeTrack on
shortLabel repeatMasker
longLabel Repeating Elements by RepeatMasker
group repeatMasker
priority 2
visibility dense
type bed 3 .
noInherit on"
bed[, {
      class <- unique(repClass_p1)
      priority <- as.numeric(priority_map[[class]])
      bed_file <- sprintf("%s/%s.%s.bed", opt$outdir, opt$prefix, class)
      fwrite(.SD, bed_file, sep = "\t", col.names = FALSE, quote = FALSE)
       track_text <- sprintf(
"
    track repeatMasker_%s
    parent RepeatMasker
    shortLabel %s
    longLabel %s Repeating Elements by RepeatMasker
    priority %d
    spectrum on
    maxWindowToDraw 10000000
    colorByStrand 50,50,150 150,50,50
    type bigBed 9
    bigDataUrl %s.%s.bb",
    class, class, class, priority, opt$prefix, class
  )
  trackdb_blocks[[class]] <<- list(
    priority = priority,
    text = track_text
  )
 
 }, by = repClass_p1]
priorities <- sapply(trackdb_blocks, `[[`, "priority")
ord <- order(priorities)
trackdb_blocks_sorted <- trackdb_blocks[ord]
child_texts <- sapply(trackdb_blocks_sorted, `[[`, "text")
final_trackdb <- c(
  parent_block,
  child_texts
)
writeLines(
  final_trackdb,
  file.path(opt$outdir, "paste2trackdb.txt")
)