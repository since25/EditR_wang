# === 在文件最开头加上 ===
options(shiny.sanitize.errors = FALSE)   # 页面直接显示报错
options(shiny.fullstacktrace = TRUE)     # 打印完整报错堆栈

server <- function(input, output, session) {
  # 初始化结果矩阵
  results_matrix <- reactiveVal(NULL)
  
  # 处理数据
  observeEvent(input$process, {
    req(input$sequence, input$ab1_files, input$target_base)
    
    # 获取用户输入的序列和目标碱基
    sequence <- toupper(input$sequence)
    target_base <- toupper(input$target_base)
    sliced_sequence <- strsplit(sequence, "")[[1]]
    row_indices <- which(sliced_sequence == target_base)
    column_indices <- c(7, 8, 9, 10)
    if (target_base %in% c("C", "A")) {
        row_indices_name <- row_indices
    } else if (target_base %in% c("G", "T")) {
        row_indices_name <- 21 - rev(row_indices)
    }
    
    # ------------------------------
    # 关键修改1: 在并行任务前生成列名
    # ------------------------------
    column_names <- c()
    for (row_index in row_indices_name) {
      for (column_index in column_indices) {
        base_perc <- switch(
          as.character(column_index),
          "7" = "A",
          "8" = "C",
          "9" = "G",
          "10" = "T",
          "Unknown"
        )
        column_names <- c(column_names, paste0("N", row_index, "_", base_perc))
      }
    }
    column_names <- c(column_names, "Table_Name")
    ncol_n <- length(column_names)
    
    # 获取文件路径和名称
    file_paths <- input$ab1_files$datapath
    file_names <- input$ab1_files$name
    
    # 设置并行计算
    # 更稳妥的并行核心设置
    #max_cores <- 4  # 建议不超过分配给虚拟机的核心数
    cl <- parallel::makeCluster(6)
    #cl <- parallel::makeCluster(max_cluster_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    
    # ------------------------------
    # 关键修改2: 显式导出必要变量
    # ------------------------------
    #parallel::clusterExport(cl, varlist = c("column_names", "row_indices", "column_indices"), envir = environment())
    #parallel::clusterExport(cl, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
    
    parallel::clusterEvalQ(cl, {
      source("global.R")
    })
    
    withProgress(message = paste0("处理文件中... ", format(Sys.time(), "%H:%M:%S")), value = 0, {
      total_files <- length(file_paths)
      progress_step <- 1 / total_files
      
      # 并行处理文件
      results_list <- foreach::foreach(
        i = 1:total_files,
        .combine =  function(...) as.data.frame(do.call(rbind, list(...))),
        .packages = c("sangerseqR", "Biostrings", "magrittr", "dplyr")
      ) %dopar% {
        file_path <- file_paths[i]
        output_file <- tempfile(pattern = paste0("result_", i, "_"), fileext = ".csv")
        
        # 处理文件
        getEditing(file_path, sequence, outputF = output_file)
        
        # 提取数据
        row_results <- rep(NA, ncol_n)
        col_counter <- 1
        for (row_index in row_indices) {
          for (column_index in column_indices) {
            cell_value <- tryCatch({
              extract_cell(output_file, row_index, column_index)
            }, error = function(e) NA)
            row_results[col_counter] <- cell_value
            col_counter <- col_counter + 1
          }
        }
        row_results[length(column_names)] <- file_names[i]
        
        if (file.exists(output_file)) file.remove(output_file)
        
        # 返回带列名的向量
        stats::setNames(row_results, column_names)
      }
      
      # 主线程更新进度
      for (i in 1:total_files) {
        incProgress(progress_step, detail = paste("已完成文件:", file_names[i]))
      }
    })
        # 👇 加上这段修正结构
    if (is.null(dim(results_list))) {
      results_list <- as.data.frame(t(results_list), stringsAsFactors = FALSE)
    } else {
      results_list <- as.data.frame(results_list, stringsAsFactors = FALSE)
    }    
    # ------------------------------
    # 关键修改3: 确保列名正确设置
    # ------------------------------
    if (target_base %in% c("G", "T")) {
       n <- ncol(results_list)
       results_list <- cbind(results_list[, (n-1):1], results_list[, n, drop = FALSE])
    }
    results <- as.data.frame(results_list)
    colnames(results) <- column_names
    
    # 更新结果矩阵
    results_matrix(results)
    #output$status <- renderText("数据处理完成！", format(Sys.time(), "%H:%M:%S"))
    output$status <- renderText(paste0("数据处理完成！时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
   })
  # 显示结果表格
  output$results <- renderTable({
    req(results_matrix())
    results_matrix()
  })
  
  # 下载结果
  output$download <- downloadHandler(
    filename = function() {
      "output.csv"
    },
    content = function(file) {
      write.csv(results_matrix(), file, row.names = FALSE)
    }
  )
}
