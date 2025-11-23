

library(shiny)
library(leaflet)
library(leaflet.extras)
library(sf)
library(rgee)
library(htmltools)
library(jsonlite)
library(shinyjs)




# Initialize rgee (assumes you have already run ee_Initialize() interactively at least once)
# If you need to run initialization programmatically: ee_Initialize()
rgee::ee_Initialize(service_account = "gee-service@pdr-jnu.iam.gserviceaccount.com",
              credentials= "service_account.json")

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Select area → get B11, LULC, NDVI → Compute PDR → Download"),
  sidebarLayout(
    sidebarPanel(
      h4("Instructions"),
      tags$ol(
        tags$li("Draw a polygon on the map (top-left drawing tools)."),
        tags$li("Click 'Run PDR' to fetch images and compute PDR."),
        tags$li("View result, download PNG or export GeoTIFF.")
      ),
      actionButton("run_pdr", "Run PDR", class = "btn-primary"),
      br(), br(),
      downloadButton("download_png", "Download PNG"),
      br(), br(),
      p("GeoTIFF export (Earth Engine): after computation you can export the GeoTIFF to your Google Drive. Click the 'Export to Drive' button below."),
      actionButton("export_drive", "Export GeoTIFF to Google Drive")
    ),
    mainPanel(
      leafletOutput("map", height = 600),
      hr(),
      h4("PDR image preview"),
      uiOutput("preview_ui"),
      verbatimTextOutput("status")
    )
  )
)

server <- function(input, output, session) {
  # initialize map with draw toolbar
  output$map <- renderLeaflet({
    leaflet() %>%
      addProviderTiles("Esri.WorldImagery") %>%
      addDrawToolbar(
        targetGroup = "draw",
        polylineOptions = FALSE,
        circleOptions = FALSE,
        editOptions = editToolbarOptions(selectedPathOptions = selectedPathOptions()),
        polygonOptions = drawPolygonOptions(shapeOptions = drawShapeOptions()),
        rectangleOptions = drawRectangleOptions(shapeOptions = drawShapeOptions())
      ) %>%
      addLayersControl(overlayGroups = c("draw"), options = layersControlOptions(collapsed = FALSE))
  })
  
  # store last drawn feature as sf
  drawn_sf <- reactiveVal(NULL)
  observeEvent(input$map_draw_new_feature, {
    feat <- input$map_draw_new_feature
    geojson <- jsonlite::toJSON(feat$geometry, auto_unbox = TRUE)
    sfobj <- sf::st_read(geojson, quiet = TRUE)
    drawn_sf(sfobj)
    leafletProxy("map") %>% clearGroup("draw") %>% addPolygons(data = sfobj, group = "draw")
  })
  
  # remove/edit events
  observeEvent(input$map_draw_edited_features, {
    # grab first feature (if multiple)
    feats <- input$map_draw_edited_features$features
    if(length(feats) >= 1) {
      geojson <- jsonlite::toJSON(feats[[1]]$geometry, auto_unbox = TRUE)
      sfobj <- sf::st_read(geojson, quiet = TRUE)
      drawn_sf(sfobj)
      leafletProxy("map") %>% clearGroup("draw") %>% addPolygons(data = sfobj, group = "draw")
    }
  })
  
  observeEvent(input$map_draw_deleted_features, {
    drawn_sf(NULL)
    leafletProxy("map") %>% clearGroup("draw")
  })
  
  # helper: compute PDR in Earth Engine for a geometry
  compute_pdr_ee <- function(ee_geom) {
    # Parameters (example weights) - you can tune these
    a_b11 <- 0.5
    b_ndvi <- 0.4
    c_lulc <- 0.1
    
    # Sentinel-2 surface reflectance collection
    s2 <- ee$ImageCollection("COPERNICUS/S2_SR")$
      filterBounds(ee_geom)$
      filterDate("2021-01-01", "2021-12-31")$   # example: year 2021; you can param.
      filter(ee$Filter$lt("CLOUDY_PIXEL_PERCENTAGE", 30))$
      map(function(img){
        # Scale factors for SR: bands are integer scaled by 1e4 (in S2_SR)
        img$select(c("B4","B8","B11"))$multiply(0.0001)$copyProperties(img, img$propertyNames())
      })$
      median()  # median composite
    
    # compute NDVI = (B8 - B4) / (B8 + B4)
    ndvi <- s2$normalizedDifference(c("B8","B4"))$rename("NDVI")
    
    # B11 band (already scaled)
    b11 <- s2$select("B11")$rename("B11")
    
    # LULC: ESA WorldCover v100 (2020) has 'Map' band
    lulc <- ee$Image("ESA/WorldCover/v100/2020")$select("Map")$rename("LULC")
    
    
    # Map LULC classes to factor values (example mapping)
    # We'll create an image lulc_factor where each class maps to a numeric factor
    # Example mapping (simplified):
    # 10 Tree cover -> 0.1, 20 Shrubland -> 0.2, 30 Grassland -> 0.3, 40 Cropland -> 0.8, 50 Built-up -> 0.1, 80 Water -> 0.0, others -> 0.2
    lulc_factor <- lulc$remap(
      from = ee$List(c(10,20,30,40,50,60,70,80,90,100)),
      to   = ee$List(c(0.1,0.2,0.3,0.8,0.1,0.2,0.2,0.0,0.2,0.2))
    )$rename("LULC_F")
    # Then optionally, set missing values to 0.2 like this:
    lulc_factor <- lulc_factor$where(lulc_factor$eq(0), 0.2)
    
    
    # Normalize B11 to 0..1 by simple min/max clip (example clip values)
    # NOTE: you might want to adjust min/max depending on local ranges; here we clip between 0 and 0.5
    b11_norm <- b11$clamp(0, 0.5)$divide(0.5)$rename("B11_NORM")
    
    # Compose PDR as weighted sum (example)
    pdr <- b11_norm$multiply(a_b11)$
      add(ndvi$multiply(b_ndvi))$
      add(lulc_factor$multiply(c_lulc))$
      rename("PDR_raw")
    
    # Optionally scale PDR to 0-1
    pdr_norm <- pdr$unitScale(0, 1)$rename("PDR")
    
    # return an Image with layers for visualization & download
    return(list(pdr = pdr_norm, vis_bands = list(b11 = b11, ndvi = ndvi, lulc = lulc)))
  }
  
  # reactive to run PDR when button clicked
  pdr_result <- eventReactive(input$run_pdr, {
    sfobj <- drawn_sf()
    if (is.null(sfobj)) {
      return(list(error = "Please draw a polygon on the map first."))
    }
    
    # convert sf to ee geometry
    # ensure CRS is EPSG:4326
    sf_ll <- sf::st_transform(sfobj, 4326)
    ee_geom <- sf_as_ee(sf_ll$geometry)  # rgee helper
    
    # call compute
    res <- tryCatch({
      computed <- compute_pdr_ee(ee_geom)
      list(ee_geom = ee_geom, pdr = computed$pdr, vis = computed$vis_bands, error = NULL)
    }, error = function(e) {
      list(error = paste0("Error computing PDR: ", e$message))
    })
    
    return(res)
  })
  
  # show status & preview using a thumbnail URL from EE
  output$status <- renderPrint({
    res <- pdr_result()
    if (is.null(res)) return("Waiting...")
    if (!is.null(res$error)) return(res$error)
    "PDR computed successfully. Preview and download ready."
    
  })
  


  # create a thumbnail URL for preview
  preview_url <- reactive({
    res <- pdr_result()
    if (is.null(res) || !is.null(res$error)) return(NULL)
    
    # Correct visualize call
    ee_img <- res$pdr$visualize(
      bands = c("PDR"),
      min = 0,
      max = 1,
      palette = c( "blue",          # low values
                   "white",
                   
                   "lavender",
                   
                   "darkgreen",
                   "lightyellow",
                   
                   "yellow",
                   "orange",
                   "maroon",
                   "red" ,
                   
                   "orangered",
                   "darkred",
                   "black" )
    )
    # build a URL tile / thumbnail
    thumb_params <- list(
      region = res$ee_geom,
      
      dimensions = 1024,
      format = "png"
    )
    url <- ee$Image$getThumbURL(ee_img, thumb_params)
    return(url)
  })
  
  output$preview_ui <- renderUI({
    url <- preview_url()
    if (is.null(url)) {
      tags$p("No preview yet.")
    } else {
      tags$div(
        tags$img(src = url, style = "max-width:100%;height:auto;border:1px solid #ddd;"),
        br(),
        tags$p("Preview generated from median Sentinel-2 (2021), NDVI, LULC. Use 'Download PNG' or 'Export to Drive' for GeoTIFF.")
      )
    }
  })
  
  # PNG download: we'll request the same thumbnail URL and serve it to the user
  output$download_png <- downloadHandler(
    filename = function() { paste0("PDR_preview_", Sys.Date(), ".png") },
    content = function(file) {
      url <- preview_url()
      if (is.null(url)) stop("No preview available.")
      # download the image to file
      httr::GET(url, httr::write_disk(file, overwrite = TRUE))
    }
  )
  
  # Export GeoTIFF to Google Drive on button click
  observeEvent(input$export_drive, {
    res <- pdr_result()
    if (is.null(res) || !is.null(res$error)) {
      showNotification("No computed PDR to export. Draw area and click Run PDR first.", type = "error")
      return()
    }
    
    # create an export task to Drive
    task_img <- res$pdr$multiply(10000)$toInt()  # scale and convert to int for export
    task <- ee$batch$Export$image$toDrive(
      image = task_img,
      description = paste0("PDR_export_", as.integer(Sys.time())),
      folder = "earthengine_exports",
      fileNamePrefix = paste0("PDR_", format(Sys.Date(), "%Y%m%d")),
      region = res$ee_geom$geometry(),
      scale = 30,
      maxPixels = 1e13
    )
    task$start()
    showNotification("Export started to your Google Drive -> folder 'earthengine_exports'. Check Earth Engine Tasks console or your Drive.", type = "message", duration = 10)
  })
  
}

shinyApp(ui, server)


