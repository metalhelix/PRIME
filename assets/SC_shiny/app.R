library(shiny)
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
      .logo {
        position: absolute;
        top: 10px;
        right: 10px;
        z-index: 10;
      }
      .title-container {
        position: relative;
        padding-top: 50px;  /* Adjust padding to avoid overlap with logo */
      }
    "))
  ),
  
  # Logo in the top right corner
	tags$div(
    class = "logo",
    tags$img(src = "logo.png", height = "70px")  # Adjust height as needed
	),

	fluidRow(
		column(3,
			selectInput(inputId="dataset", label="Dataset", choices = c("enh6_bud-3ss","enh6_8-10ss","enh4_bud-4ss","enh4_9-12ss"), selected = "enh6_bud-3ss"),
			textInput("gene", label="Gene Symbol", value="aldob"), 
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

server<-function(input, output){
  
selected_dataset <- eventReactive(input$submit, {input$dataset})

  dataInput <- reactive({
    switch(selected_dataset(),
	"enh6_bud-3ss"=readRDS("/n/core/Bioinformatics/secondary//Sauka-Spengler/sa2722/MOLNG-3950.GRCz11.Ens_110/L66525/PRIME_SC_out/SeuratObj.rds"), "enh6_8-10ss"=readRDS("/n/core/Bioinformatics/secondary//Sauka-Spengler/sa2722/MOLNG-3950.GRCz11.Ens_110/L66526/PRIME_SC_out/SeuratObj.rds"), "enh4_bud-4ss"=readRDS("/n/core/Bioinformatics/secondary//Sauka-Spengler/sa2722/MOLNG-3950.GRCz11.Ens_110/L66527/PRIME_SC_out/SeuratObj.rds"), "enh4_9-12ss"=readRDS("/n/core/Bioinformatics/secondary//Sauka-Spengler/sa2722/MOLNG-3950.GRCz11.Ens_110/L66528/PRIME_SC_out/SeuratObj.rds")
    )
  })

  selected_gene <- eventReactive(input$submit, {input$gene})

  checkgene <- reactive({
  obj <- dataInput()
  
  # Perform the validation
  validate(
    need(selected_gene() %in% rownames(obj@assays$SCT), "Gene Symbol/ID not yet found in data set or has been filtered out")
  )
  selected_gene()  # Return the selected gene if valid
})

  selected_clust <- eventReactive(input$submit, {input$clust})

  checkclust <- reactive({
    obj=dataInput()
    validate(need((selected_clust() == "all" | selected_clust() %in% unique(obj@meta.data$seurat_clusters)),"Cluster not yet found in data set"))
	selected_clust()
  })

  output$umapPlot<- renderPlot({
	obj=dataInput()
	checkclust()
	dp <- DimPlot(object = obj,reduction="umap")
	orig_xlim<-ggplot_build(dp)$layout$panel_scales_x[[1]]$range$range
	orig_ylim<-ggplot_build(dp)$layout$panel_scales_y[[1]]$range$range
	b<-ggplot2::ggplot_build(dp)
	cols<-unique(b[[1]][[1]][,c("group","colour")])
	cols<-cols[order(cols[,1]),]

	if(input$clust =="all")
	{
		p <- DimPlot(object = obj,reduction="umap",label=T, label.size = 6,pt.size=.05)
		p<- p + xlim(orig_xlim) + ylim(orig_ylim)
	}else
	{
		p <- DimPlot(object = obj,reduction="umap",cols=cols[cols[,1]==(as.numeric(input$clust)+1),2],cells=colnames(obj)[obj@meta.data$seurat_clusters==input$clust],label=T,label.size = 6,pt.size=.05)
		p<- p + xlim(orig_xlim) + ylim(orig_ylim)
	}
	print(p)
  })

  output$featurePlot <- renderPlot({
		obj=dataInput()
		checkgene()
		# mygene <- getensid(input$gene)
		mygene <- input$gene
		dp <- DimPlot(object = obj,reduction="umap")
		orig_xlim<-ggplot_build(dp)$layout$panel_scales_x[[1]]$range$range
		orig_ylim<-ggplot_build(dp)$layout$panel_scales_y[[1]]$range$range
		if(mygene %in% rownames(obj@assays$SCT)){
			p2 <- FeaturePlot(obj, reduction="umap",features=mygene, cols = c("#BEBEBE32","blue"), pt.size = .1)
			p2 + labs(title=mygene,sep='') + xlim(orig_xlim) + ylim(orig_ylim)
		}
  })

  output$violinPlot <- renderPlot({
        obj=dataInput()
        checkgene()
        # mygene <- getensid(input$gene)
		mygene <- input$gene
		# mysym <-getsym(mygene)
		mysym <-mygene
		mytitle<-paste(mygene," (",mysym,")",sep='')
        if(mygene %in% rownames(obj@assays$SCT)){
          p <- VlnPlot(obj, mygene,pt.size=.05)
		  p + labs(title=mytitle,sep='')
        }
  })

  output$umapPlot_orig.ident<- renderPlot({
	obj=dataInput()
	p <- DimPlot(object = obj,reduction="umap",pch.use=20,pt.size=.1,group.by="orig.ident")
	print(p)
  })
}

shinyApp(ui=ui, server=server)

