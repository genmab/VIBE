#' Plot a boxplot pre and post treatment per gene per indication with Kruskal Wallis test, Fold change and sample counts
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import ggsignif
#'
#'
#' @param df dataframe harmonized by VIBE's harmonize_df function with paired pre/post info
#' @param dataset name of dataset to be mentioned for captioning and sorting plot-axis.
#' @param symbol_list list of genes for creating the plot
#' @param method the summary stat method to use, either "mean" or "median"
#' @param plot_groups Character vector with values you wish to include in the plot, threshold will be calculated for target1 and target2 on all data in df. plot_groups values must be present in the grouping_var column.
#' @param treatment_levels pre and post-treatment levels
#'
#' @return boxplot comparing the statistical differences of the sample distribution per gene, per indication
#'
#' @examples
#' # prepare data
#'data("VIBE_data")
#'df_harmonized = harmonize_df(df = VIBE_data,
#'                             col_patientid = "patient_no",
#'                           col_sampleid = "analysis",
#'                           col_indication = "tumor",
#'                           col_treatment = "treatment_flag",
#'                           col_symbol = "gene",
#'                           col_value = "log2_tpm",
#'                           unit = "tpm",
#'                           keep_extra_cols = c("database"))
#'
#'
#' # example1:
#'plot_box_pre_post(df = df_harmonized,
#'                  dataset = "VIBE dummy data (not real data!)",
#'                  symbol_list = c("Immune target", "Tumor target", "Gene 1"),
#'                  method = "median")
#' # example2: limit indications plotted to selection
#'plot_box_pre_post(df = df_harmonized,
#'                  dataset = "VIBE dummy data (not real data!)",
#'                  symbol_list = c("Immune target", "Tumor target", "Gene 1"),
#'                  method = "median", plot_groups = c("CERV", "NSCLC", "SARC"))

#' # example3: calculate fold change based on mean values instead of median
#'plot_box_pre_post(df = df_harmonized,
#'                  dataset = "VIBE dummy data (not real data!)",
#'                  symbol_list = c("Immune target", "Tumor target", "Gene 1"),
#'                  method = "mean", plot_groups = c("CERV", "NSCLC", "SARC"))

plot_box_pre_post <- function (df,
                                   dataset = NA, # was datatype
                                   symbol_list,
                               method = c("mean", "median"),
                               plot_groups = NA,
                               treatment_levels = c("pre-treatment", "post-treatment")){
  # filter data so only groups of interest are plotted
  if(!is.na(plot_groups[1])) {
    df <- filter(df, indication %in% plot_groups)
  }

  ifelse(!is.na(dataset), dataset <- dataset, dataset <- "Dataset")

  # extract the rows with corresponding gene names
  df_sel <- ungroup(df) %>%
    dplyr::filter(symbol %in% symbol_list) %>%
    mutate(treatment = factor(treatment, levels = treatment_levels)) %>%
    arrange(symbol, treatment, patientid)

  # Calclate p values based on Kruskal test and posthoc test
  kw <- df_sel %>%
    group_by(indication,symbol) %>%
    summarise(p=round(kruskal.test(value ~ treatment)$p.value,2),
              y=max(value),
              x=1.5) %>%
    ungroup() %>%
    mutate(y=max(y))

  # Calculate fold change
  fc <- df_sel %>%
    arrange(treatment) %>%
    group_by(indication,symbol, treatment)  %>%
    summarise(stat = if(method == "mean"){mean(value)} else if (method == "median"){median(value)},
              y = max(value)) %>%
    group_by(indication, symbol) %>%
    summarise(fold_change = stat[2]/stat[1],
              y = max(y),
              x = 1.5) %>%
    ungroup() %>%
    mutate(y = max(y))

   ## Calculate count for the groups (changed as code above counts patients double; 1x pre and 1x post treatmnet)
  count_samples <- df_sel %>%
    arrange(treatment) %>%
    group_by(indication,symbol, treatment)  %>%
    summarize(count=n(),
              y=min(value),
              x=1.5)%>%
    ungroup() %>%
    mutate(y=min(y)) %>%
    select(-treatment) %>%
    distinct()

  ### Caption for the plot
  caption = paste0("Dataset: ", dataset, "\n",
                   "KW: ", " P vales are calculated based on the Kruskal Wallis test\n",
                   "FC: ", " Fold change  between pre-treatment and post-treatment condition for each gene per indication \n",
                   "n: ", " Number of samples in pre-treatment and post-treatment condition for each gene per indication")

  title = paste("Expression of genes pre- and post-treatment in", dataset, "dataset")
  counttype = df_sel$unit[1]

  # 6. Plot with two facets
  boxplot_pre_post <-
    ggplot(df_sel, aes(x=treatment, y=value)) +
    geom_boxplot((aes(color=treatment, fill=treatment)), alpha=0.5) +
    geom_line(aes(group=patientid), linetype="11",linewidth=0.9, alpha=0.5, colour = gmb_palette[["grey"]])+
    geom_point(size=0.5, alpha=0.8,colour = gmb_palette[["grey"]])+
    scale_color_manual(values=c(gmb_palette[["lightgreen"]],gmb_palette[["yellow"]]))+
    scale_fill_manual(values=c(gmb_palette[["lightgreen"]], gmb_palette[["yellow"]]))+
    facet_grid(symbol~indication, switch = "y") +
    labs(title = title, caption = caption)+
    ylim(min(df_sel$value)-2,max(fc$y)*1.2)+
    theme(
      # axis test, line and titles
      axis.text.x = element_text(color="black",size=12, angle=30, hjust=1),
      axis.text.y=element_text(color="black",size=12),
      axis.line = element_line(color="grey",linewidth=0.5),
      axis.title=element_text(size=20, face="bold", color="black"),
      legend.position = "None",
      #panel background, margins and title
      panel.background = element_blank(),
      plot.margin = unit(c(1,3,1,1), "lines"),
      plot.title = element_text(size = 20, face = "bold"),
      #facet text and background setting
      strip.background = element_rect(linewidth = 0.5, color="black",fill="white"),
      strip.placement = "outside",
      strip.text.x = element_text(size=16, face="bold",color="black"),
      strip.text.y = element_text(angle=0,size=16, face="bold",color="black"),
      plot.caption = element_text(size=12))+
    #X and Y axis titles
    xlab("Indication and treatment groups") +
    ylab(paste0("Expression log2(",counttype,")"))+
    # Geom text from Kruskal Wallis test and Fold change
    geom_text(data=kw,aes(x=x, y=y*1.04, label=paste0("KW p=",p))) +
    geom_text(data=fc,aes(x=x, y=y*1.15, label=paste0("FC=",round(fold_change,2))))+
    geom_text(data=count_samples, aes(x=x, y=y-0.2, label=paste0("n=",count)))

  return (boxplot_pre_post)
}

