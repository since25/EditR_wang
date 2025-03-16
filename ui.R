ui <- fluidPage(
  titlePanel("AB1 文件分析工具"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("可用于tBE的sanger测序分析，请注意输入序列必须与测序方向一致，即5-3方向序列，若运行失败，请尝试复制互补序列。分析结果与ab1测序质量相关，结果与预期不符时，请导入snapgene查看峰图是否异常"),
      textInput("sequence", "输入序列", value = ""),  # 序列输入框
      textInput("target_base", "输入目标碱基（A/T/G/C）", value = "C"),  # 目标碱基输入框
      fileInput("ab1_files", "上传 AB1 文件", multiple = TRUE, accept = ".ab1"),  # 文件上传按钮
      actionButton("process", "处理数据")  # 处理数据按钮
    ),
    
    mainPanel(
      textOutput("status"),  # 状态信息
      verbatimTextOutput("task_progress"),  # 当前任务进度
      tableOutput("results"),  # 结果表格
      downloadButton("download", "下载结果")  # 下载按钮
    )
  )
)
