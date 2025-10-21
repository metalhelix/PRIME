#!/usr/bin/env python

import sys
import argparse
import subprocess
import os
import shutil
import pandas as pd


parser = argparse.ArgumentParser()

parser.add_argument('--app_location', default="/n/core/Bioinformatics/ShinyApps")
parser.add_argument('--second_folder')
# parser.add_argument('--molngID')
parser.add_argument('--wwwPath')
args=parser.parse_args()

destination_folder = f"{args.app_location}/{str(args.second_folder).split('Bioinformatics/secondary/')[-1]}"
os.makedirs(destination_folder, exist_ok=True)

df_summary = pd.read_csv(f"{args.second_folder}/scripts/summary.csv")
rds_list=[]
markerGene_list = []
for index, row in df_summary.iterrows():
    rds_list.append(f"\"{row['sampleName']}\"=readRDS(\"{args.second_folder}/{row['libID']}/PRIME_SC_out/SeuratObj.rds\")")
    markerGene_list.append(f"{args.second_folder}/{row['libID']}/PRIME_SC_out/Seurat_Report_gen_plots/AllMarkers.csv")
    
samples = ",".join([f"\"{v}\"" for v in df_summary["sampleName"].to_list()])
display_sample = f"\"{df_summary['sampleName'].to_list()[0]}\""
rds_files = ", ".join(rds_list) 

df_marker_genes = pd.read_csv(markerGene_list[0]) #The first library will be auto displayed in the app, so we select for marker genes of the first library to auto display
display_marker_gene = [f"\"{v}\"" for v in df_marker_genes["gene"].to_list() if "ENS" not in v][0]

shiny_app = f"""library(shiny)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(shinybusy)
library(bslib)

options(shiny.maxRequestSize=50*1024^2) # Do we still need this? No file to upload here

ui <- fluidPage(
	add_busy_spinner(spin = "fading-circle"),
	titlePanel("PRIME Seurat Visualization"),
	theme = bslib::bs_theme(version = 4, bootswatch = "minty" ),

	tags$head(
    tags$style(HTML("
      .logo {{
        position: absolute;
        top: 10px;
        right: 10px;
        z-index: 10;
      }}
      .title-container {{
        position: relative;
        padding-top: 50px;  /* Adjust padding to avoid overlap with logo */
      }}
    "))
  ),
  
  # Logo in the top right corner
	tags$div(
    class = "logo",
    tags$img(src = "logo.png", height = "70px")  # Adjust height as needed
	),

	fluidRow(
		column(3,
			selectInput(inputId="dataset", label="Dataset", choices = c({samples}), selected = {display_sample}),
			textInput("gene", label="Gene Symbol", value={display_marker_gene}), 
			selectInput(inputId="clust", label="Cluster", choices = c("all",0:30), selected = "all"),
            actionButton("submit", "Submit", class = "btn-sm btn-success", width = "100px",style='height:50px')  
		),
		column(8,
			tags$h4("Expression violin plot"),
			plotOutput("violinPlot")
		)
	),

	fluidRow(
		column(6,
			tags$h4("UMAP plot"),
			plotOutput("umapPlot",width="800px",height="700px")
		),
		column(5,
			tags$h4("Feature plot"),
			plotOutput("featurePlot",width="800px",height="700px")
		)
	),

	tags$div(
		"ShinyApp stored under: /n/core/Bioinformatics/ShinyApps/",
    	style = "
			position: fixed;
			bottom: 10px;
			right: 10px;
			color: red;
			font-size: 20px;
			font-weight: bold;
			"
		)
)

server<-function(input, output){{
    
    selected_dataset <- reactiveVal({display_sample})
	observeEvent(input$submit, {{selected_dataset(input$dataset)}})

	selected_gene <- reactiveVal({display_marker_gene})
	observeEvent(input$submit, {{selected_gene(input$gene)}})
	
	selected_clust <- reactiveVal("all")
	observeEvent(input$submit, {{selected_clust(input$clust)}})

	dataInput <- reactive({{
		switch(selected_dataset(),
		{rds_files}
		)
	}})


	checkgene <- reactive({{
    obj=dataInput()
	validate(need(selected_gene() %in% rownames(obj@assays$SCT),"Gene Symbol/ID not yet found in data set or has been filtered out"))
    selected_gene()
	}})

	checkclust <- reactive({{
		obj=dataInput()
		validate(need((selected_clust() == "all" | selected_clust() %in% unique(obj@meta.data$seurat_clusters)),"cluster not yet found in data set"))
		selected_clust()
	}})

	output$umapPlot<- renderPlot({{
		obj=dataInput()
		checkclust()
		dp <- DimPlot(object = obj,reduction="umap")
		orig_xlim<-ggplot_build(dp)$layout$panel_scales_x[[1]]$range$range
		orig_ylim<-ggplot_build(dp)$layout$panel_scales_y[[1]]$range$range
		b<-ggplot2::ggplot_build(dp)
		cols<-unique(b[[1]][[1]][,c("group","colour")])
		cols<-cols[order(cols[,1]),]
		my_clust <- selected_clust()

	if(my_clust =="all")
	{{
		p <- DimPlot(object = obj,reduction="umap",label=T, label.size = 6,pt.size=.05)
		p<- p + xlim(orig_xlim) + ylim(orig_ylim)
	}}else
	{{
		p <- DimPlot(object = obj,reduction="umap",cols=cols[cols[,1]==(as.numeric(my_clust)+1),2],cells=colnames(obj)[obj@meta.data$seurat_clusters==my_clust],label=T,label.size = 6,pt.size=.05)
		p<- p + xlim(orig_xlim) + ylim(orig_ylim)
	}}
	print(p)
	}})

	output$featurePlot <- renderPlot({{
		obj=dataInput()
		checkgene()
		# mygene <- getensid(input$gene)
		mygene <- selected_gene()
		dp <- DimPlot(object = obj,reduction="umap")
		orig_xlim<-ggplot_build(dp)$layout$panel_scales_x[[1]]$range$range
		orig_ylim<-ggplot_build(dp)$layout$panel_scales_y[[1]]$range$range
		if(mygene %in% rownames(obj@assays$SCT)){{
			p2 <- FeaturePlot(obj, reduction="umap",features=mygene, cols = c("#BEBEBE32","blue"), pt.size = .1)
			p2 + labs(title=mygene,sep='') + xlim(orig_xlim) + ylim(orig_ylim)
		}}
	}})

	output$violinPlot <- renderPlot({{
        obj=dataInput()
        checkgene()
		mygene <- selected_gene()
        mytitle <- mygene
        if(mygene %in% rownames(obj@assays$SCT)){{
          p <- VlnPlot(obj, mygene,pt.size=.05)
		  p + labs(title=mytitle,sep='')
        }}
	}})

	output$umapPlot_orig.ident<- renderPlot({{
		obj=dataInput()
		p <- DimPlot(object = obj,reduction="umap",pch.use=20,pt.size=.1,group.by="orig.ident")
		print(p)
	}})
}}

shinyApp(ui=ui, server=server)

"""

newfile = open(os.path.join(destination_folder,"app.R"), "w")
newfile.write(shiny_app)
newfile.close()

shutil.copytree(args.wwwPath, f"{destination_folder}/www", dirs_exist_ok=True)