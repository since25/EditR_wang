
#完美版本1



server <- function(input, output, session) {
  # 初始化结果矩阵
  results_matrix <- reactiveVal(NULL)
  
  # 处理数据
  observeEvent(input$process, {
    req(input$sequence)  # 确保序列已输入
    req(input$ab1_files)  # 确保文件已上传
    req(input$target_base)  # 确保目标碱基已输入
    
    # 获取用户输入的序列和目标碱基
    sequence <- input$sequence
    target_base <- toupper(input$target_base)  # 强制转换为大写
    sliced_sequence <- strsplit(sequence, "")[[1]]
    row_indices <- c()
    column_indices <- c(7, 8, 9, 10)
    
    # 找到所有目标碱基的位置
    for (i in 1:length(sliced_sequence)) {
      if (sliced_sequence[i] == target_base) {
        row_indices <- c(row_indices, i)
      }
    }
    
    # 初始化结果矩阵
    ncol_n <- length(row_indices) * length(column_indices) + 1
    results <- matrix(ncol = ncol_n, nrow = nrow(input$ab1_files))
    
    # 设置列名
    column_names <- c()
    for (row_index in row_indices) {
      for (column_index in column_indices) {
        # 将列索引映射为对应的碱基百分比名称
        base_perc <- switch(
          as.character(column_index),
          "7" = "A.perc",
          "8" = "C.perc",
          "9" = "G.perc",
          "10" = "T.perc",
          "Unknown"  # 如果列索引不在 7,8,9,10 中，标记为 Unknown
        )
        column_names <- c(column_names, paste0("Row_", row_index, "_", base_perc))
      }
    }
    column_names <- c(column_names, "Table_Name")
    colnames(results) <- column_names
    
    # 获取上传的文件路径
    file_paths <- input$ab1_files$datapath
    
    # 使用 withProgress 显示任务进度
    withProgress(message = "处理文件中...", value = 0, {
      # 顺序处理文件
      for (i in 1:length(file_paths)) {
        file_path <- file_paths[i]
        
        # 更新当前任务进度
        incProgress(1 / length(file_paths), detail = paste("正在处理文件:", input$ab1_files$name[i]))
        
        # 打印当前正在处理的文件路径（用于调试）
        print(paste("Processing file:", file_path))
        
        # 调用 getEditing 函数处理 .ab1 文件
        output_file <- tempfile(fileext = ".csv")  # 临时文件用于存储结果
        getEditing(file_path, sequence, outputF = output_file)
        
        # 提取指定单元格内容
        col_counter <- 1
        for (row_index in row_indices) {
          for (column_index in column_indices) {
            results[i, col_counter] <- extract_cell(output_file, row_index, column_index)
            col_counter <- col_counter + 1
          }
        }
        
        # 存储文件名
        results[i, ncol_n] <- input$ab1_files$name[i]
      }
    })
    
    # 更新结果矩阵
    results_matrix(results)
    output$status <- renderText("数据处理完成！")
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


