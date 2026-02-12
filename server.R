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
    
    # 设置并行计算 - 动态检测可用核心数
    # 获取系统可用核心数，预留1个核心给系统
    sys_cores <- parallel::detectCores()
    max_cores <- min(sys_cores - 1, 6)  # 不超过6个核心，预留1个给系统
    max_cores <- max(1, max_cores)  # 至少1个核心

    # 如果核心数不足3个，切换到顺序处理以提高稳定性
    use_parallel <- max_cores >= 2

    if (use_parallel) {
      cl <- parallel::makeCluster(max_cores)
      doParallel::registerDoParallel(cl)
      on.exit({
        tryCatch({
          parallel::stopCluster(cl)
        }, error = function(e) {})
      })
    }

    # 获取脚本所在目录，确保集群节点能正确加载global.R
    script_dir <- getwd()
    global_r_path <- file.path(script_dir, "global.R")

    # 为每个并行进程创建独立的工作目录，避免临时文件冲突
    temp_dir <- file.path(tempdir(), paste0("editr_process_", Sys.getpid()))
    dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)

    withProgress(message = paste0("处理文件中... ", format(Sys.time(), "%H:%M:%S")), value = 0, {
      total_files <- length(file_paths)
      progress_step <- 1 / total_files

      if (use_parallel) {
        # 在集群节点上加载global.R，使用完整路径
        parallel::clusterEvalQ(cl, {
          source(global_r_path)
        })

        # 并行处理文件
        results_list <- foreach::foreach(
          i = 1:total_files,
          .combine = function(...) as.data.frame(do.call(rbind, list(...))),
          .packages = c("sangerseqR", "Biostrings", "magrittr", "dplyr")
        ) %dopar% {
          # 每个进程使用独立的临时目录
          proc_temp_dir <- file.path(temp_dir, paste0("worker_", Sys.getpid()))
          dir.create(proc_temp_dir, showWarnings = FALSE, recursive = TRUE)

          file_path <- file_paths[i]
          output_file <- file.path(proc_temp_dir, paste0("result_", i, "_", Sys.getpid(), ".csv"))

          tryCatch({
            # 处理文件
            getEditing(file_path, sequence, outputF = output_file)

            # 提取数据
            row_results <- rep(NA, ncol_n)
            col_counter <- 1
            for (row_index in row_indices) {
              for (column_index in column_indices) {
                cell_value <- tryCatch({
                  extract_cell(output_file, row_index, column_index)
                }, error = function(e) {
                  warning(paste("提取单元格失败:", row_index, column_index, e$message))
                  NA
                })
                row_results[col_counter] <- cell_value
                col_counter <- col_counter + 1
              }
            }
            row_results[length(column_names)] <- file_names[i]

            # 清理临时文件
            if (file.exists(output_file)) {
              tryCatch(file.remove(output_file), error = function(e) {})
            }

            # 返回带列名的向量
            stats::setNames(row_results, column_names)
          }, error = function(e) {
            warning(paste("处理文件失败:", file_path, e$message))
            # 返回NA向量
            stats::setNames(rep(NA, ncol_n), column_names)
          })
        }
      } else {
        # 顺序处理文件（更稳定）
        results_list <- NULL
        for (i in 1:total_files) {
          file_path <- file_paths[i]
          output_file <- file.path(temp_dir, paste0("result_", i, ".csv"))

          tryCatch({
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

            # 清理临时文件
            if (file.exists(output_file)) {
              tryCatch(file.remove(output_file), error = function(e) {})
            }

            row_result_df <- as.data.frame(t(stats::setNames(row_results, column_names)),
                                          stringsAsFactors = FALSE)

            if (is.null(results_list)) {
              results_list <- row_result_df
            } else {
              results_list <- rbind(results_list, row_result_df)
            }
          }, error = function(e) {
            warning(paste("处理文件失败:", file_path, e$message))
            na_row <- as.data.frame(t(rep(NA, ncol_n)), stringsAsFactors = FALSE)
            colnames(na_row) <- column_names
            na_row[1, length(column_names)] <- file_names[i]

            if (is.null(results_list)) {
              results_list <- na_row
            } else {
              results_list <- rbind(results_list, na_row)
            }
          })

          # 更新进度
          incProgress(progress_step, detail = paste("已完成文件:", file_names[i]))
        }
      }

      # 确保results_list是数据框
      if (use_parallel) {
        if (is.null(dim(results_list))) {
          results_list <- as.data.frame(t(results_list), stringsAsFactors = FALSE)
        } else {
          results_list <- as.data.frame(results_list, stringsAsFactors = FALSE)
        }
      }

      # 确保列名正确设置
      if (target_base %in% c("G", "T")) {
         n <- ncol(results_list)
         if (n > 1) {
           results_list <- cbind(results_list[, (n-1):1], results_list[, n, drop = FALSE])
         }
      }
      results <- as.data.frame(results_list)
      colnames(results) <- column_names

      # 更新结果矩阵
      results_matrix(results)
      output$status <- renderText(paste0("数据处理完成！时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    })
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
