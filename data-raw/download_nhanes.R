dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

nhanes_files <- data.frame(
  file = c("DEMO_L.XPT", "BMX_L.XPT", "DIQ_L.XPT", "GHB_L.XPT", "MCQ_L.XPT"),
  url = paste0(
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/",
    c("DEMO_L.XPT", "BMX_L.XPT", "DIQ_L.XPT", "GHB_L.XPT", "MCQ_L.XPT")
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(nhanes_files))) {
  destination <- file.path("data/raw", nhanes_files$file[i])
  if (!file.exists(destination)) {
    message("Downloading ", nhanes_files$file[i])
    download.file(
      nhanes_files$url[i], destination,
      mode = "wb", quiet = FALSE
    )
  }
}

write.csv(
  nhanes_files,
  "data/derived/nhanes_download_manifest.csv",
  row.names = FALSE
)

