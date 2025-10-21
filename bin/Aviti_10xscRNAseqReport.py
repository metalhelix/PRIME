#! /n/ngs/tools/primary_smk/config/primary_smk/bin/python

import subprocess as sb
import requests as req
import pandas as pd
import glob
import os
import argparse
import re
import urllib3
import numpy as np
import sys

req.packages.urllib3.util.ssl_.DEFAULT_CIPHERS = "ALL:@SECLEVEL=1"
pd.set_option("display.max_columns", None)

parser = argparse.ArgumentParser(
    description="Combine multiple non-FLEX 10x reports into a single report"
)
parser.add_argument(
    "--secondary_path",
    help="Path to secondary analysis results folder for the order of interest",
)
args = parser.parse_args()

##get path and order ID
order_path = args.secondary_path
os.chdir(args.secondary_path)

order_id = re.search("MOLNG-\\d\\d\\d\\d", order_path).group(0)
print(order_id)

if not order_id.startswith("MOLNG-"):
    print(
        "Could not find Order ID, make sure you run the script inside the order folder"
    )
    exit()

myToken = "493d945ab7b64bcd1bc05996c8371e26"
head = {"x-zan-apitoken": myToken, "Accept": "application/json"}

shiny_server = f"http://shiny:3431/Bioinformatics/{str(args.secondary_path).split('Bioinformatics/secondary/')[-1]}"

# get Flowcell Names from an order


def flowcell_name(order_id):
    machine_flowcells = []
    url = (
        "http://lims.stowers.org/zanmodules/molecular-biology/ngs/requests/%s/flowcells"
        % order_id
    )
    response = req.get(url, headers=head)
    if response.status_code == 200:
        data = response.json()
        machinetype = data[0]["machineMode"]["label"]
        machine_flowcells.append(str(machinetype))
        # print(data[0]['machineMode']['label'])
        for f in data:
            machine_flowcells.append(str(f["ngsCellId"]))
    return machine_flowcells


flowcell = flowcell_name(order_id)[1:]


# get sample Names and info into a dataframe
def get_flowcellSamples(flowcells):
    # list_1 =[]
    flowcell_data = []

    for f in flowcells:
        url = (
            "http://lims.stowers.org/zanmodules/molecular-biology/ngs/flowcells/%s/samples"
            % f
        )
        response = req.get(url, headers=head)
        if response.status_code == 200:
            data = response.json()
            for samples in data["samples"]:
                df = pd.DataFrame.from_dict(
                    dict([(k, pd.Series(v)) for k, v in samples.items()])
                )
                df["FlowcellName"] = f
                df = df[df["prnOrderNo"] == order_id]
                flowcell_data.append(df)
    df_concat = pd.concat(flowcell_data, ignore_index=True)
    return df_concat


df_samples = get_flowcellSamples(flowcell)

# Pull Lims Info
genome = df_samples["speciesName"][0]
ordertype = df_samples["orderType"][0]
objective = df_samples["analysisGoals"][0]
machinetype = flowcell_name(order_id)[0]
link = "https://lims.stowers.org/#/molecular-biology/ngs/requests/" + order_id

# Extract only necesaary columns, drop duplicates, sort by Sample, Create a directory called scripts and write data to targets.csv
df_samples = df_samples[["prnOrderNo", "sampleName", "FlowcellName", "libID"]]
df_samples["FlowcellName"] = df_samples["FlowcellName"].astype(str)
df_samples = df_samples.drop_duplicates()
df_samples = df_samples.sort_values(by=["libID"])
sb.call("mkdir -p scripts/", shell=True)
df_samples.to_csv("scripts/targets.tsv", sep="\t", header=True, index=False)

# Read web_summmary csv files and extract info

allFiles = glob.glob(order_path + "/*/outs/metrics_summary.csv")
list_summary = []
for file_ in allFiles:
    filename_split = str(file_).strip().split("/")
    df_temp = pd.read_csv(
        file_,
        index_col=None,
        header=0,
        usecols=[
            "Estimated Number of Cells",
            "Mean Reads per Cell",
            "Median Genes per Cell",
            "Number of Reads",
            "Valid Barcodes",
            "Sequencing Saturation",
            "Reads Mapped to Genome",
            "Fraction Reads in Cells",
            "Median UMI Counts per Cell",
        ],
    )
    df_temp["name"] = filename_split[-3]
    df_temp["10X web summary"] = "https://webfs" + file_.replace(
        "metrics_summary.csv", "web_summary.html"
    )
    df_temp["CellBender Summary"] = (
        "https://webfs" + file_.split("outs")[0] + "cellbender_outs/cellbender_output_report.html"
    )
    df_temp["PRIME Seurat summary"] = (
        "https://webfs" + file_.split("outs")[0] + "PRIME_SC_out/SeuratReport.html"
    )
    list_summary.append(df_temp)

# concate df's
frame_summary = pd.concat(list_summary, axis=0, ignore_index=True)

# Add sample names to summary file
df_sub_samples = df_samples[["sampleName", "libID"]].drop_duplicates()
merged_summary = pd.merge(
    frame_summary, df_sub_samples, how="left", left_on="name", right_on="libID"
).drop_duplicates()
df_summary = merged_summary[
    [
        "libID",
        "sampleName",
        "10X web summary",
        "PRIME Seurat summary",
        "CellBender Summary",
        "Median UMI Counts per Cell",
        "Estimated Number of Cells",
        "Mean Reads per Cell",
        "Median Genes per Cell",
        "Number of Reads",
        "Reads Mapped to Genome",
        "Fraction Reads in Cells",
    ]
]
df_summary = df_summary.sort_values(by=["libID"])
df_summary.dropna(subset=["libID"], inplace=True)
# print(len(df_summary))
df_summary.to_csv("scripts/summary.csv", sep=",", header=True, index=False)

if len(df_summary) > 0:
    # Write R markdown script
    script = """
---
title: "10x scRNA-seq Analysis Report"
output:
  html_document:
    code_folding: hide
    collapsed: yes
    highlight: tango
    number_sections: no
    theme: cerulean
    toc: yes
    toc_float: yes
fontsize: 30pt
---

``` {r setup, echo=FALSE, message=FALSE, results="hide", warning=FALSE}
library(rmdformats)
library(prettydoc)
library(hrbrthemes)
library(tint)
library(tufte)

library(ggplot2)
library(RColorBrewer)
library(pander)
library(knitr)
library(scales)
library("gridExtra")
library("cowplot")
theme_x <- theme_bw() + theme(text = element_text(size=14))

panderOptions('table.style','rmarkdown')
panderOptions('table.split.table',Inf)


options(width=50)
opts_chunk$set(error=FALSE)
opts_chunk$set(echo=FALSE)
opts_chunk$set(tidy=FALSE)
opts_chunk$set(warnings=FALSE)

```


## [%s](%s)

### Machine Type
%s

%s


### Method
* [Cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) count version 8.0.1 was used in the initial single cell analysis of the samples. </br>
* [CellBender](https://www.nature.com/articles/s41592-023-01943-7) version 0.3.0 were used on h5 files output by Cellranger for further denoising. </br>
* [Seurat](https://satijalab.org/seurat/) version 5.1.0 was used for data filteration, normalization, and single cell analysis. </br>

### Objective
%s


```{r table, echo=TRUE}
df_targets <- read.table('scripts/targets.tsv', sep = '\\t',header = TRUE, as.is = TRUE)
colnames(df_targets) <- c('Order', 'Sample Name', 'Flowcell', 'Library')
Counts <- row.names(df_targets)
df_targets <- cbind(Counts, df_targets)
df_targets$Flowcell <- as.character(df_targets$Flowcell)
pander(df_targets)
```

## ShinyApp
You can explore your data via a preconstructed [ShinyApp](%s)

## Web Summary 

### Summary Table
```{r web_summary, echo=TRUE}
df_summary <- read.csv('scripts/summary.csv', sep = ',',header = TRUE, as.is = TRUE)
colnames(df_summary) <- c("Library", "Sample Name", "10X web summary", 'PRIME Seurat summary', 'CellBender Summary',"Median UMI Counts per Cell", "Estimated Number of Cells" , "Mean Reads per Cell","Median Genes per Cell", "Number of Reads" ,  "Reads Mapped to Genome" , "Fraction Reads in Cells"  )

df_summary$`Estimated Number of Cells` = as.numeric(gsub(",", "", df_summary$`Estimated Number of Cells`))
df_summary$`Number of Reads` = as.numeric(gsub(",", "", df_summary$`Number of Reads`))
df_summary$`Mean Reads per Cell` = as.numeric(gsub(",", "", df_summary$`Mean Reads per Cell`))
df_summary$`Median Genes per Cell` = as.numeric(gsub(",", "",df_summary$`Median Genes per Cell`))

df_summary_overall <- df_summary[, c("Library", "Sample Name", "10X web summary", 'CellBender Summary', 'PRIME Seurat summary')]

df_summary_overall$`10X web summary` <- paste( "<a href=" , df_summary_overall$`10X web summary`, ' target=\\"_blank\\"> ', "10X_Summary  </a>", sep ="")
df_summary_overall$`PRIME Seurat summary` <- paste( "<a href=" , df_summary_overall$`PRIME Seurat summary`, ' target=\\"_blank\\"> ', "Seurat_Summary  </a>", sep ="")
df_summary_overall$`CellBender Summary` <- paste( "<a href=" , df_summary_overall$`CellBender Summary`, ' target=\\"_blank\\"> ', "CellBender_Summary  </a>", sep ="")

pander(df_summary_overall)
```
</br>

### 10X Stats Summary Table
```{r web_summary_10X, echo=TRUE}
df_summary_10x <- df_summary[, c("Library", "Sample Name", "Median UMI Counts per Cell", "Estimated Number of Cells" , "Mean Reads per Cell","Median Genes per Cell", "Number of Reads" , "Reads Mapped to Genome" , "Fraction Reads in Cells")]
pander(df_summary_10x)
```
</br>

### Initial Seurat Stats Summary Table
```{r web_summary_seurat, echo=TRUE}
seurat_obj_list <- df_summary$`PRIME Seurat summary` %%>%%
    gsub("SeuratReport.html", "SeuratObj.rds", .) %%>%%
    gsub("https://webfs", "", .)

Total_number_of_cells_list <- list()
Mean_reads_per_cell_list <- list()
Median_genes_per_cell_list <- list()

for (i in 1:length(seurat_obj_list)) {
    seurat_obj <- readRDS(file=seurat_obj_list[i])
    Mean_reads_per_cell_list[i] <- mean(seurat_obj@meta.data$nCount_RNA)
    Total_number_of_cells_list[i] <- length(seurat_obj@meta.data$nCount_RNA)
    Median_genes_per_cell_list[i] <- median(seurat_obj@meta.data$nFeature_RNA)
    }

df_seurat <- data.frame(matrix(ncol = 0, nrow = length(seurat_obj_list)))
df_seurat$libID <- df_summary$libID
df_seurat$`Sample Name` <- df_summary$`Sample Name`
df_seurat <- cbind(df_seurat, `Estimated Number of Cells`= unlist(Total_number_of_cells_list))
df_seurat <- cbind(df_seurat, `Mean Reads per Cell`= unlist(Mean_reads_per_cell_list))
df_seurat <- cbind(df_seurat, `Median Genes per Cell`= unlist(Median_genes_per_cell_list))
pander(df_seurat)
```
<br/>

### Summary Plots
```{r number_cells, echo=TRUE, warning=FALSE, message=FALSE, fig.width=12, fig.height=12}
df_summary_10x$source <- "10X"
df_seurat$source <- "Cellbender + Seurat"
df_sum_cell_number <- rbind(df_summary_10x[,c("Sample Name", "Estimated Number of Cells", "source")], df_seurat[,c("Sample Name", "Estimated Number of Cells", "source")])

p1 <- ggplot(df_sum_cell_number, aes(x = `Sample Name`, y = `Estimated Number of Cells`, fill = source)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  scale_y_continuous(name="Estimated Number of Cells", 
                     labels=comma) + 
  xlab("Sample Name") +
  ggtitle("Estimated Number of Cells") +
  theme_x +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5), legend.position = "top") +
  guides(fill=guide_legend(title ="Pipelines", title.position = "top"))+
  theme(plot.margin=unit(c(0,0.5,0.5,0.2),"cm"))


df_sum_median_gene <- rbind(df_summary_10x[,c("Sample Name", "Median Genes per Cell", "source")], df_seurat[,c("Sample Name", "Median Genes per Cell", "source")])

p2 <- ggplot(df_sum_median_gene, aes(x = `Sample Name`, y = `Median Genes per Cell`, fill = source)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  scale_y_continuous(name="Median Genes per Cell", 
                     labels=comma) + 
  xlab("Sample Name") +
  ggtitle("Median Genes per Cell") +
  theme_x +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5), legend.position = "top") +
  guides(fill=guide_legend(title ="Pipelines", title.position = "top"))+
  theme(plot.margin=unit(c(0,0.5,0.5,0.2),"cm"))


df_sum_mean_reads <- rbind(df_summary_10x[,c("Sample Name", "Mean Reads per Cell", "source")], df_seurat[,c("Sample Name", "Mean Reads per Cell", "source")])

p3 <- ggplot(df_sum_mean_reads, aes(x = `Sample Name`, y = `Mean Reads per Cell`, fill = source)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  scale_y_continuous(name="Mean Reads per Cell", 
                     labels=comma) + 
  xlab("Sample Name") +
  ggtitle("Mean Reads per Cell") +
  theme_x +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5), legend.position = "top") +
  guides(fill=guide_legend(title ="Pipelines", title.position = "top"))+
  theme(plot.margin=unit(c(0,0.5,0.5,0.2),"cm"))
  
  
#Function to Extract legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#Extract legend from plot one(p1)
com_legend <- get_legend(p1)

#Remove the legend from plot one (p1)
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")

#Matrix of plots
grid.arrange( com_legend, p1,p2,p3,  ncol=2, nrow = 3, 
              layout_matrix = rbind(c(1,1), c(2,3), c(4,5)),
              widths = c(7,7), heights = c(1.5,5, 5))
             

```

## Individual datasets
The following table provides direct links to directories of Cellranger & Cellbender output files for each dataset.
```{r}
dat <- df_summary[,c('Sample Name', '10X web summary', 'PRIME Seurat summary')]
dat$`10X web summary` <- gsub("/web_summary.html", "", dat$`10X web summary`)
dat$`10X web summary` <- gsub("10X_Summary", "Cellranger Folder", dat$`10X web summary`)
dat$`PRIME Seurat summary` <- gsub("/PRIME_SC_out/SeuratReport.html", "/cellbender_outs/", dat$`PRIME Seurat summary`)
dat$`PRIME Seurat summary` <- gsub("Seurat_Summary", "Cellbender Folder", dat$`PRIME Seurat summary`)

colnames(dat) <- c("Sample ID", "10X Results Directory", "Cellbender Results Directory")
row.names(dat) <- NULL
pander(dat)

```


</br>
(Contact Bowen of CompBio if link does not work)

###  Session information

**Contact:** Please email Bowen (by2747@stowers.org) if you have questions or suggestions

**Generated:** `r format(Sys.time(), "%%a %%b %%d %%Y, %%I:%%M %%p")`

For reproducibility, this analysis was performed with Cell Ranger, R and Python session:

```{r version, echo=FALSE, comment=NA}

message("Cell Ranger v7.1.0")

# genome_version

message("Genomes: %s")
#genome
```

``` {r session_info, echo=FALSE, comment=NA}
sessionInfo()
```

""" % (order_id, link, machinetype, ordertype, objective, shiny_server, genome)

    file = open("analysis_10xscRNASeq.Rmd", "w")
    file.write(script)
    file.close()

    cmd_rmarkdown = "R -e \"rmarkdown::render('analysis_10xscRNASeq.Rmd',output_file='analysis_10xscRNASeq.html')\""
    sb.call(cmd_rmarkdown, shell=True)

else:
    sys.exit(
        "summary table is empty, check 'sample names' and 'libID' fetched from LIMS"
    )
