library(ggplot2)
library(ggrepel)
library(Cairo)

args = commandArgs(trailingOnly=TRUE)

segment_csv_path <- args[1]
segment_svg_path <- args[2]
len_target_seq <- as.numeric(args[3])
max_len_no_split <- as.numeric(args[4])

cat("Load data\n")

dat <- data.frame(length=c(0,len_target_seq), heigth=c(0,0))
segment_data <- read.csv(segment_csv_path, quote = "'")
segment_data$color <- factor(segment_data$color, levels = c(3,2,1,-1,-2,-3))

if(len_target_seq > max_len_no_split){
  lower_bound <- -1.1
}else{
  lower_bound <- -0.7
}

g <- (
  ggplot(dat, aes(x=length, y=heigth, group=1)) +
    geom_line(arrow=arrow(angle = 90, ends="both")) +
    scale_y_continuous(limits=c(lower_bound,0.7), n.breaks=20) +
    geom_segment(data=segment_data, aes(x=start, y=height_arrow, xend=end,
                                         yend=height_arrow, color=color),
                 lineend = "butt", linejoin = "mitre", size=2, show.legend = FALSE,
                 arrow=arrow(angle=45, length=unit(0.3,"cm"),
                             type='open', ends=segment_data$arrow_end)) +
    # geom_line is only for legend the + 1000 makes sure that the lines are not visible
    geom_line(data=segment_data, aes(x=start, y=height_arrow + 1000, color=color), size = 3,) +
    scale_colour_manual(values = c("red3", "royalblue", "lawngreen", "purple", "orange", "turquoise"), drop = FALSE) +
    guides(color=guide_legend(title="Frame")) +
    theme_void()
)

if(nrow(segment_data) < 12){
 g <- g + geom_text_repel(data=segment_data, aes(x=start, y = height_text, label=label, segment.size = 0), point.size = NA, direction = "x")
}

if(len_target_seq > max_len_no_split){
  start_child <- seq(0, len_target_seq, max_len_no_split/2)
  end_child <- seq(max_len_no_split, len_target_seq, max_len_no_split/2)
  if(end_child[length(end_child)] != len_target_seq){
    end_child <- c(end_child, len_target_seq)
  }
  start_child <- start_child[1:length(end_child)]
  child_iterator <- seq(1,length(start_child))
  height <- rep(c(-0.9,-1.1), length(start_child))[1:length(start_child)]
  child_dat <- data.frame(start=start_child, end=end_child, height=height, child_iterator=child_iterator)

  g <- (g + geom_segment(data=child_dat, aes(x=start, y=height, xend=end,
                                         yend=height),
                 lineend = "butt", linejoin = "mitre", show.legend = FALSE,
                 arrow=arrow(angle=90, length=unit(0.2,"cm"), 
                             type='open', ends="both"))
 )
 if (length(start_child) < 15){
    g <- (g +
          geom_text(data=child_dat, aes(x=(end_child+start_child)/2, y=height+0.05, label=child_iterator)) +
          annotate(geom="text", label="Child Sections", x=len_target_seq/2,  y=-0.75)
    )
 } else{
    g <- g + annotate(geom="text", label="Child Sections", x=len_target_seq/2,  y=-0.75)
 }
}

CairoSVG(segment_svg_path, width=7, height=4)
g
