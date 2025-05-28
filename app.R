library(shiny)
library(rentrez)
library(ggplot2)
library(plotly)
library(DT)
library(dplyr)
library(stringr)
library(htmltools)
library(rmarkdown)

ui <- fluidPage(
  titlePanel("Shiny Genomics Explorer"),
  sidebarLayout(
    sidebarPanel(
      # Nuevo panel para carga de archivos
      h4("Cargar datos desde archivo"),
      fileInput("file_input", "Seleccionar archivo FASTA o CSV",
                accept = c(".fasta", ".fa", ".csv"),
                buttonLabel = "Examinar..."),
      actionButton("load_file_btn", "Cargar Datos desde Archivo"),
      helpText("Formatos aceptados: FASTA (.fasta, .fa) o CSV con columnas 'Accession' y 'seq'"),
      hr(),
      
      textInput("search_term", "Término de búsqueda (gen, especie, etc.):", value = "BRCA1 Homo sapiens"),
      actionButton("search_btn", "Buscar Accession IDs"),
      br(), br(),
      DTOutput("acc_table"),
      hr(),
      actionButton("load_btn", "Cargar Datos de Secuencias"),
      hr(),
      sliderInput("min_len", "Longitud Mínima de Secuencia:", min = 0, max = 50000, value = 0, step = 100),
      selectInput("plot_type", "Tipo de Gráfico:",
                  choices = c("Distribución de Longitudes", "Contenido GC", "Frecuencia de nucleótidos")),
      numericInput("kmer_size", "Tamaño de k-mer:", value = 3, min = 1, max = 6),  # Reducido el máximo para mejor rendimiento
      hr(),
      h4("Visualización Avanzada"),
      sliderInput("window_size", "Tamaño de ventana GC:", min = 10, max = 1000, value = 100, step = 10),
      sliderInput("window_step", "Paso ventana GC:", min = 1, max = 500, value = 50, step = 1),
      selectInput("density_type", "Densidad:", choices = c("Longitud", "% GC")),
      selectInput("kmer_vis", "Visualizar k-mer:", choices = c("Tabla", "PCA"), selected = "Tabla"),
      hr(),
      h4("Exportación y Reportes"),
      downloadButton("download_summary", "Descargar Resumen CSV"),
      downloadButton("download_kmer", "Descargar k-mer CSV"),
      downloadButton("download_plot", "Descargar Gráfico PNG"),
      selectInput("report_format", "Formato reporte:", choices = c("HTML", "PDF", "Word")),
      downloadButton("download_report", "Generar Reporte"),
      hr(),
      helpText("Nota: Para grandes conjuntos de datos, el análisis puede tardar varios minutos.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Secuencias",       DTOutput("seq_table")),
        tabPanel("Visor Secuencias", uiOutput("seq_viewer")),
        tabPanel("Análisis k-mer",   DTOutput("kmer_table")),
        tabPanel("Gráficos",         
                 tabsetPanel(
                   tabPanel("Histograma",     plotlyOutput("plot")),
                   tabPanel("Densidad",       plotlyOutput("density_plot")),
                   tabPanel("Heatmap GC",     plotlyOutput("heatmap_gc")),
                   tabPanel("PCA k-mer",      plotlyOutput("pca_plot"))
                 )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # 1. Buscar Accession IDs con manejo de errores
  acc_df <- reactiveVal(data.frame())
  observeEvent(input$search_btn, {
    req(input$search_term)
    tryCatch({
      withProgress(message = "Buscando en NCBI...", value = 0.3, {
        res <- entrez_search(db = "nucleotide", term = input$search_term, retmax = 50)  # Reducido a 50 resultados
        incProgress(0.6)
        if (length(res$ids) == 0) {
          showNotification("No se encontraron resultados.", type = "warning")
          acc_df(data.frame())
          return()
        }
        sums <- entrez_summary(db = "nucleotide", id = res$ids)
        df <- data.frame(
          Accession = sapply(sums, function(x) x$accessionversion),
          Descripción = sapply(sums, function(x) x$title),
          stringsAsFactors = FALSE
        ) %>% arrange(Accession)
        acc_df(df)
        incProgress(1)
      })
    }, error = function(e) {
      showNotification(paste("Error al buscar en NCBI:", e$message), type = "error")
      acc_df(data.frame())
    })
  })
  
  # 2. Selección de Accessions
  output$acc_table <- renderDT({
    df <- acc_df()
    req(nrow(df) > 0)
    datatable(df, selection = "multiple", rownames = FALSE,
              options = list(pageLength = 5, autoWidth = TRUE, scrollX = TRUE))
  })
  
  # 3. Carga de secuencias con manejo de errores
  seqs_info <- reactive({
    input$load_btn
    isolate({
      sel <- input$acc_table_rows_selected
      req(length(sel) > 0)
      ids <- acc_df()[sel, "Accession"]
      
      tryCatch({
        withProgress(message = "Descargando secuencias...", value = 0, {
          fasta <- entrez_fetch(db = "nucleotide", id = ids, rettype = "fasta")
          incProgress(0.5)
          recs <- unlist(strsplit(fasta, "\n>"))
          
          seq_list <- lapply(1:length(recs), function(i) {
            incProgress(i/length(recs)/2)
            lines <- strsplit(recs[i], "\n")[[1]]
            header <- gsub("^>", "", lines[1])
            acc <- strsplit(header, " ")[[1]][1]
            seq <- paste(lines[-1], collapse = "")
            data.frame(Accession = acc, seq = seq, stringsAsFactors = FALSE)
          })
          
          do.call(rbind, seq_list)
        })
      }, error = function(e) {
        showNotification(paste("Error al cargar secuencias:", e$message), type = "error")
        return(data.frame())
      })
    })
  })
  
  # 4. Estadísticas básicas con validación
  summary_df <- reactive({
    df <- seqs_info()
    req(nrow(df) > 0)
    
    tryCatch({
      df %>% mutate(
        Longitud = nchar(seq),
        GC = round(str_count(toupper(seq), "[GC]") / nchar(seq) * 100, 2)
      ) %>% filter(Longitud >= input$min_len)
    }, error = function(e) {
      showNotification("Error al calcular estadísticas", type = "error")
      return(data.frame())
    })
  })
  
  output$seq_table <- renderDT({
    df <- summary_df()
    req(nrow(df) > 0)
    datatable(df %>% select(Accession, Longitud, GC), 
              options = list(dom = 't', scrollX = TRUE))
  })
  
  # 5. Análisis k-mer optimizado
  kmer_df <- reactive({
    df <- seqs_info()
    req(nrow(df) > 0)
    k <- as.integer(input$kmer_size)
    req(k > 0)
    
    tryCatch({
      s <- toupper(paste0(df$seq, collapse = ""))
      n <- nchar(s)
      if (k > n) return(data.frame())
      
      # Limitar a las primeras 50,000 bases para rendimiento
      if (n > 50000) {
        s <- substr(s, 1, 50000)
        n <- 50000
        showNotification("Secuencia truncada a 50,000 bases para análisis k-mer", type = "warning")
      }
      
      kmers <- substring(s, 1:(n - k + 1), k:n)
      tab <- sort(table(kmers), decreasing = TRUE)
      data.frame(kmer = names(tab), Frecuencia = as.integer(tab), stringsAsFactors = FALSE)
    }, error = function(e) {
      showNotification("Error en análisis k-mer", type = "error")
      return(data.frame())
    })
  })
  
  output$kmer_table <- renderDT({
    df <- kmer_df()
    req(nrow(df) > 0)
    datatable(df, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # 6. Gráfico principal mejorado
  plot_obj <- reactive({
    df <- summary_df()
    req(nrow(df) > 0)
    
    tryCatch({
      switch(input$plot_type,
             "Distribución de Longitudes" = {
               ggplot(df, aes(Longitud)) + 
                 geom_histogram(binwidth = 100, fill = "steelblue", alpha = 0.7) + 
                 labs(title = "Distribución de Longitudes", x = "Longitud", y = "Frecuencia")
             },
             "Contenido GC" = {
               ggplot(df, aes(x = reorder(Accession, GC), y = GC)) + 
                 geom_bar(stat = "identity", fill = "darkgreen") + 
                 labs(title = "Contenido GC (%)", x = "Accession", y = "% GC") + 
                 coord_flip() + 
                 theme(axis.text.y = element_text(size = 8))
             },
             "Frecuencia de nucleótidos" = {
               seq_all <- toupper(paste0(df$seq, collapse = ""))
               ft <- as.data.frame(table(strsplit(seq_all, "")[[1]]), stringsAsFactors = FALSE)
               colnames(ft) <- c("Nucleótido", "Frecuencia")
               ggplot(ft, aes(Nucleótido, Frecuencia, fill = Nucleótido)) + 
                 geom_bar(stat = "identity") + 
                 labs(title = "Frecuencia de nucleótidos") +
                 scale_fill_manual(values = c(A = "red", T = "blue", C = "green", G = "orange"))
             })
    }, error = function(e) {
      showNotification("Error al generar gráfico", type = "error")
      return(ggplot() + labs(title = "Error al generar gráfico"))
    })
  })
  
  output$plot <- renderPlotly({ 
    p <- plot_obj()
    ggplotly(p) %>% config(displayModeBar = TRUE)
  })
  
  # 7. Gráficos avanzados
  
  # 7a. Densidad
  output$density_plot <- renderPlotly({
    df <- summary_df()
    req(nrow(df) > 0)
    
    tryCatch({
      if (input$density_type == "Longitud") {
        p <- ggplot(df, aes(x = Longitud)) + 
          geom_density(fill = "steelblue", alpha = 0.5) + 
          labs(title = "Densidad de Longitudes", x = "Longitud")
      } else {
        p <- ggplot(df, aes(x = GC)) + 
          geom_density(fill = "darkgreen", alpha = 0.5) + 
          labs(title = "Densidad de % GC", x = "% GC")
      }
      ggplotly(p)
    }, error = function(e) {
      showNotification("Error en gráfico de densidad", type = "error")
      return(plotly_empty())
    })
  })
  
  # 7b. Heatmap GC optimizado
  output$heatmap_gc <- renderPlotly({
    df <- seqs_info()
    req(nrow(df) > 0)
    
    ws <- input$window_size
    step <- input$window_step
    max_windows <- 100  # Más reducido para mejor rendimiento
    
    tryCatch({
      withProgress(message = "Calculando contenido GC...", value = 0, {
        heat_df <- do.call(rbind, lapply(1:min(nrow(df), 10), function(i) {  # Limitar a 10 secuencias
          incProgress(1/min(nrow(df), 10), detail = paste("Procesando", df$Accession[i]))
          
          seq <- toupper(df$seq[i])
          len <- nchar(seq)
          if(len < ws) return(NULL)
          
          starts <- seq(1, len - ws + 1, by = step)
          if(length(starts) > max_windows) {
            starts <- round(seq(1, max(starts), length.out = max_windows))
          }
          
          gc <- sapply(starts, function(pos) {
            frag <- substr(seq, pos, pos + ws - 1)
            round(str_count(frag, "[GC]") / ws * 100, 2)
          })
          
          data.frame(
            Accession = df$Accession[i],
            Position = starts + floor(ws/2),
            GC = gc,
            stringsAsFactors = FALSE
          )
        }))
      })
      
      if(is.null(heat_df) || nrow(heat_df) == 0) return(NULL)
      
      plot_ly(
        data = heat_df,
        x = ~Position,
        y = ~Accession,
        z = ~GC,
        type = "heatmap",
        colors = colorRamp(c("blue", "white", "red")),
        colorbar = list(title = "% GC"),
        hoverinfo = "text",
        text = ~paste("Accession:", Accession,
                      "<br>Position:", Position,
                      "<br>GC%:", GC)
      ) %>%
        layout(
          title = "Contenido GC por Ventana Deslizante",
          xaxis = list(title = "Posición Genómica"),
          yaxis = list(title = "Accession", categoryorder = "category ascending"),
          margin = list(l = 150)
        )
    }, error = function(e) {
      showNotification("Error en heatmap GC", type = "error")
      return(plotly_empty())
    })
  })
  
  # 7c. PCA optimizado para k-mers
  output$pca_plot <- renderPlotly({
    df <- seqs_info()
    req(nrow(df) > 0)
    k <- as.integer(input$kmer_size)
    req(k > 0 && k <= 6)  # Limitado a k=6 para rendimiento
    
    tryCatch({
      # Limitar a las primeras 5000 bases por secuencia para PCA
      seqs_short <- lapply(df$seq, function(s) substr(toupper(s), 1, 5000))
      
      # Calcular frecuencias de k-mer para cada secuencia (limitado a k-mers encontrados)
      kmers_list <- lapply(seqs_short, function(s) {
        n <- nchar(s)
        if(n < k) return(table(character(0)))
        kmers <- substring(s, 1:(n - k + 1), k:n)
        table(kmers)
      })
      
      # Obtener todos los k-mers únicos
      all_kmers <- unique(unlist(lapply(kmers_list, names)))
      
      # Crear matriz de frecuencias
      mat <- do.call(rbind, lapply(kmers_list, function(x) {
        counts <- rep(0, length(all_kmers))
        names(counts) <- all_kmers
        counts[names(x)] <- x
        counts
      }))
      
      # Filtrar k-mers muy raros
      col_sums <- colSums(mat)
      mat <- mat[, col_sums > 1 & nchar(colnames(mat)) == k]
      
      if(ncol(mat) < 2) {
        showNotification("No hay suficientes k-mers para PCA", type = "warning")
        return(plotly_empty())
      }
      
      # PCA
      pca <- prcomp(mat, scale. = TRUE)
      scores <- as.data.frame(pca$x[,1:2])
      scores$Accession <- df$Accession
      
      p <- ggplot(scores, aes(x = PC1, y = PC2, text = Accession)) +
        geom_point(color = "steelblue", size = 3) +
        labs(title = paste("PCA de", k, "-mers"), x = "PC1", y = "PC2")
      
      ggplotly(p, tooltip = "text") %>% config(displayModeBar = TRUE)
    }, error = function(e) {
      showNotification("Error en análisis PCA", type = "error")
      return(plotly_empty())
    })
  })
  
  # 8. Visor de secuencias coloreadas optimizado
  output$seq_viewer <- renderUI({
    df <- seqs_info()
    req(nrow(df) > 0)
    
    # Limitar a 5 secuencias para visualización
    if(nrow(df) > 5) {
      showNotification("Mostrando solo las primeras 5 secuencias", type = "warning")
      df <- df[1:5, ]
    }
    
    tryCatch({
      seq_divs <- lapply(1:nrow(df), function(i) {
        # Limitar a 200 caracteres por secuencia para visualización
        seq <- strsplit(toupper(substr(df$seq[i], 1, 200)), "")[[1]]
        spans <- lapply(seq, function(b) {
          color <- switch(b, A = "red", T = "blue", C = "green", G = "orange", "black")
          tags$span(style = paste0("color:", color), b)
        })
        tags$div(
          tags$b(df$Accession[i], " (primeros 200 nucleótidos)"),
          tags$div(style = "overflow-x:auto; white-space: pre; font-family: monospace; margin-bottom:1em;",
                   do.call(tagList, spans))
        )
      })
      do.call(tagList, seq_divs)
    }, error = function(e) {
      tags$div("Error al visualizar secuencias", style = "color:red;")
    })
  })
  
  # 9. Descargas
  
  output$download_summary <- downloadHandler(
    filename = function() paste0("resumen_genomica_", Sys.Date(), ".csv"),
    content = function(file) {
      tryCatch({
        write.csv(summary_df(), file, row.names = FALSE)
      }, error = function(e) {
        showNotification("Error al generar archivo CSV", type = "error")
      })
    }
  )
  
  output$download_kmer <- downloadHandler(
    filename = function() paste0("kmer_frecuencias_", Sys.Date(), ".csv"),
    content = function(file) {
      tryCatch({
        write.csv(kmer_df(), file, row.names = FALSE)
      }, error = function(e) {
        showNotification("Error al generar archivo k-mer", type = "error")
      })
    }
  )
  
  output$download_plot <- downloadHandler(
    filename = function() paste0("grafico_genomica_", Sys.Date(), ".png"),
    content = function(file) {
      tryCatch({
        p <- plot_obj()
        ggsave(file, plot = p, device = "png", width = 10, height = 7, dpi = 300)
      }, error = function(e) {
        showNotification("Error al guardar gráfico", type = "error")
      })
    }
  )
  
  # 10. Generación de reporte con plantilla integrada
  output$download_report <- downloadHandler(
    filename = function() {
      paste0("reporte_genomica_", Sys.Date(), 
             switch(input$report_format,
                    "HTML" = ".html",
                    "PDF" = ".pdf",
                    "Word" = ".docx"))
    },
    content = function(file) {
      tryCatch({
        # Crear reporte temporal con plantilla integrada
        tempReport <- file.path(tempdir(), "report.Rmd")
        file.copy("report_template.Rmd", tempReport, overwrite = TRUE)
        
        # Si no existe el template, crear uno básico
        if(!file.exists(tempReport)) {
          writeLines(c(
            "---",
            "title: 'Reporte de Análisis Genómico'",
            "output:", 
            "  html_document: default",
            "  pdf_document: default",
            "  word_document: default",
            "params:",
            "  resumen: NULL",
            "  kmers: NULL",
            "  ws: 100",
            "  step: 50",
            "---",
            "",
            "```{r setup, include=FALSE}",
            "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
            "```",
            "",
            "# Resumen de Secuencias",
            "",
            "```{r}",
            "if(!is.null(params$resumen) && nrow(params$resumen) > 0) {",
            "  knitr::kable(head(params$resumen), caption = 'Resumen de secuencias')",
            "} else {",
            "  cat('No hay datos de secuencias disponibles')",
            "}",
            "```",
            "",
            "# Análisis de k-mers",
            "",
            "```{r}",
            "if(!is.null(params$kmers) && nrow(params$kmers) > 0) {",
            "  knitr::kable(head(params$kmers), caption = 'Frecuencias de k-mers')",
            "} else {",
            "  cat('No hay datos de k-mers disponibles')",
            "}",
            "```",
            "",
            "# Parámetros de análisis",
            "",
            "- Tamaño de ventana GC: `r params$ws`",
            "- Paso de ventana GC: `r params$step`",
            "",
            "Reporte generado el `r Sys.Date()`"
          ), tempReport)
        }
        
        # Parámetros para el reporte
        params <- list(
          resumen = summary_df(),
          kmers = kmer_df(),
          ws = input$window_size,
          step = input$window_step
        )
        
        # Renderizar el reporte
        render(tempReport, output_file = file,
               params = params,
               envir = new.env(parent = globalenv()),
               output_format = switch(input$report_format,
                                      "HTML" = html_document(),
                                      "PDF" = pdf_document(),
                                      "Word" = word_document())
        )
      }, error = function(e) {
        showNotification(paste("Error al generar reporte:", e$message), type = "error")
      })
    }
  )
}

shinyApp(ui, server)