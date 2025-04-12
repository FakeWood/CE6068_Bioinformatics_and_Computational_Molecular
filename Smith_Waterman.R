library(ggplot2)
library(reshape2)
library(grid)

smith_waterman_plot <- function(seq1, seq2, match = 2, mismatch = -1, gap = -2) {
  n <- nchar(seq1)
  m <- nchar(seq2)
  s1 <- unlist(strsplit(seq1, ""))
  s2 <- unlist(strsplit(seq2, ""))
  
  score <- matrix(0, n + 1, m + 1)
  trace <- matrix("", n + 1, m + 1)
  max_score <- 0
  max_pos <- c(1, 1)
  
  for (i in 1:(n + 1)) {
    score[i, 1] <- 0
    trace[i, 1] <- "U"
  }
  for (j in 1:(m + 1)) {
    score[1, j] <- 0
    trace[1, j] <- "L"
  }
  trace[1, 1] <- "0"
  
  for (i in 2:(n + 1)) {
    for (j in 2:(m + 1)) {
      match_score <- ifelse(s1[i - 1] == s2[j - 1], match, mismatch)
      diag <- score[i - 1, j - 1] + match_score
      up <- score[i - 1, j] + gap
      left <- score[i, j - 1] + gap
      
      score[i, j] <- max(diag, up, left, 0)
      trace[i, j] <- switch(
        which.max(c(diag, up, left, 0)),
        "D", "U", "L", "0"
      )
      
      if (score[i, j] > max_score) {
        max_score <- score[i, j]
        max_pos <- c(i,j)
      }
    }
  }
  
  # 回溯
  align1 <- ""
  align2 <- ""
  i <- max_pos[1]
  j <- max_pos[2]
  path <- list(c(i, j))  # 記錄回溯路徑
  
  while (i > 1 || j > 1) {
    if (score[i, j] <= 0) {
      break
    }
    
    dir <- trace[i, j]
    if (dir == "D") {
      align1 <- paste0(s1[i - 1], align1)
      align2 <- paste0(s2[j - 1], align2)
      i <- i - 1
      j <- j - 1
    } else if (dir == "U") {
      align1 <- paste0(s1[i - 1], align1)
      align2 <- paste0("-", align2)
      i <- i - 1
    } else {
      align1 <- paste0("-", align1)
      align2 <- paste0(s2[j - 1], align2)
      j <- j - 1
    }
    path <- append(path, list(c(i, j)))
  }
  
  cat("Alignment:\n")
  cat(align1, "\n")
  cat(align2, "\n")
  cat("\nFinal Score:", score[n + 1, m + 1], "\n")
  
  # reshape for ggplot
  df <- melt(score)
  colnames(df) <- c("i", "j", "score")
  
  # 回溯箭頭
  arrows <- do.call(rbind, lapply(seq(length(path) - 1), function(k) {
    from <- path[[k + 1]]
    to <- path[[k]]
    data.frame(
      x = from[2], y = -from[1],
      xend = to[2], yend = -to[1]
    )
  }))
  
  # 顯示矩陣 + 箭頭 + 分數
  p <- ggplot(df, aes(x = j, y = -i)) +
    geom_tile(aes(fill = score), color = "white") +
    geom_text(aes(label = score), size = 4) +
    geom_segment(data = arrows,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.15, "inches")),
                 color = "red", linewidth = 1) +
    scale_fill_gradient(low = "white", high = "skyblue") +
    coord_fixed() +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  # 加上序列標籤（字母）
  # seq2 -> x 軸
  p <- p + annotate("text", x = 2:(m + 1), y = 0, label = s2, size = 5, fontface = "bold") +
    annotate("text", x = 0, y = -(2:(n + 1)), label = s1, size = 5, fontface = "bold")
  
  print(p)
}

# 測試
smith_waterman_plot("ACGTGACCAGTT", "CGTGACCGGTT")
