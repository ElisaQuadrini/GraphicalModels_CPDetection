# ==============================================================================
# Extract and Process Kenneth R. French 10 Industry Portfolios (Daily)
# EXTENDED HORIZON: 2005 - 2025 to Capture Multiple Change-Points
# Run from the repository root: source("application/1_eda_realdata.R")
# ==============================================================================

library(lubridate)
library(readr)
library(GGally)
library(ggplot2)
library(tidyr)
library(dplyr)
library(corrplot)

# Ensure output directories exist (data and figures are not tracked as empty
# folders by git, so this call is idempotent and safe to keep in the script)
dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/application", showWarnings = FALSE, recursive = TRUE)


# --- Step 1: Download and extract raw archive ---------------------------------
dataset_url <- "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/10_Industry_Portfolios_daily_CSV.zip"
temp_zip    <- tempfile(fileext = ".zip")
temp_dir    <- tempfile()

download.file(dataset_url, destfile = temp_zip, mode = "wb")
unzip(temp_zip, exdir = temp_dir)

csv_filepath <- list.files(temp_dir, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE)
raw_text_lines <- read_lines(csv_filepath)

# --- Step 2: Locate and parse Value-Weighted Daily Returns -------------------
header_index <- which(grepl("NoDur", raw_text_lines))[1]
data_start_index <- header_index

# Find the first blank line after the data block to locate table boundary
blank_indices <- which(trimws(raw_text_lines) == "")
data_end_index <- blank_indices[blank_indices > data_start_index][1] - 1

# Define target industry sector codes
target_industries <- c("NoDur", "Durbl", "Manuf", "Enrgy", "HiTec", 
                       "Telcm", "Shops", "Hlth", "Utils")

# Read comma-separated values using I() to avoid literal string warnings
daily_returns_raw <- read_csv(
  I(paste(raw_text_lines[data_start_index:data_end_index], collapse = "\n")),
  col_types = cols(.default = col_double(), ...1 = col_character()),
  show_col_types = FALSE
)

# Rename the first column (date column in French's standard layout)
colnames(daily_returns_raw)[1] <- "date_raw"

# --- Step 3: Clean dates and standard missing-value codes ---------------------
daily_clean <- daily_returns_raw %>%
  mutate(date = ymd(date_raw)) %>%
  filter(!is.na(date)) %>%
  dplyr::select(date, all_of(target_industries)) %>%  # <--- Risolto il conflitto qui con dplyr::
  mutate(across(all_of(target_industries), ~ na_if(.x, -99.99))) %>%
  mutate(across(all_of(target_industries), ~ na_if(.x, -99.990))) %>%
  mutate(across(all_of(target_industries), ~ na_if(.x, -999)))

# --- Step 4: Filter target time frame (EXPANDED: Jan 2005 - Dec 2025) --------
start_date <- as.Date("2005-01-01")
end_date   <- as.Date("2025-12-31")

daily_filtered <- daily_clean %>%
  filter(date >= start_date, date <= end_date)

# --- Step 5: Aggregate daily returns into weekly log-returns -----------------
weekly_returns <- daily_filtered %>%
  mutate(week = floor_date(date, unit = "week", week_start = 1)) %>%  # Monday start
  group_by(week) %>%
  summarise(
    across(all_of(target_industries), ~ log(prod(1 + .x / 100))),
    trading_days = n(),
    .groups = "drop"
  ) %>%
  filter(trading_days >= 3) %>%  # Exclude incomplete holiday/boundary weeks
  dplyr::select(-trading_days)   

# --- Step 6: Construct data matrix Y for model input ------------------------
Y <- as.matrix(weekly_returns[, target_industries])
rownames(Y) <- as.character(weekly_returns$week)

# --- Diagnostics check --------------------------------------------------------
cat("=== EXTENDED MATRIX SUMMARY ===\n")
cat("Matrix Dimensions (T x p):", dim(Y), "\n")
cat("Total Missing (NA) values:", sum(is.na(Y)), "\n\n")

# save the data matrix for use by application/2_application_realdata.R
saveRDS(Y, file = "data/industry_portfolios_2005_2025.rds")

# ==============================================================================
# EXPLORATORY DATA ANALYSIS (EDA) - EXACT GGPLOT REPLICATION
# ==============================================================================

# --- Format Data for ggplot2 --------------------------------------------------
df_y <- as.data.frame(Y)
df_y$Date <- as.Date(rownames(Y))

df_long <- df_y %>%
  pivot_longer(cols = -Date, names_to = "Sector", values_to = "Log_Return")

# --- PLOT 1: Weekly Log-Returns -----------------------------------------------
# Lines styled exactly like screenshot 1
p1 <- ggplot(df_long, aes(x = Date, y = Log_Return, color = Sector)) +
  geom_line(alpha = 0.8, linewidth = 0.6) +
  geom_vline(xintercept = as.Date("2020-03-09"), linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_y_continuous(breaks = seq(-0.2, 0.2, by = 0.1), limits = c(min(df_long$Log_Return, -0.2), max(df_long$Log_Return, 0.2))) +
  theme_minimal() +
  labs(
    title = paste0("Weekly Log-Returns across 9 Sectors (", format(min(df_long$Date), "%Y"), " - ", format(max(df_long$Date), "%Y"), ")"),
    subtitle = "Red dashed line indicates COVID-19 market crash (March 2020)",
    x = "Date",
    y = "Log Return",
    color = "Sector"
  ) +
  theme(
    plot.title = element_text(size = 13, face = "plain", hjust = 0),
    plot.subtitle = element_text(size = 11, hjust = 0),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p1)


# --- PLOT 2: 4-Week Rolling Volatility ----------------------------------------
rolling_sd <- apply(Y, 2, function(col) {
  zoo::rollapply(col, width = 4, FUN = sd, fill = NA, align = "right")
})
df_roll <- as.data.frame(rolling_sd)
df_roll$Date <- as.Date(rownames(Y))

df_roll_long <- df_roll %>%
  pivot_longer(cols = -Date, names_to = "Industry", values_to = "Rolling_SD")

p2 <- ggplot(df_roll_long, aes(x = Date, y = Rolling_SD, color = Industry)) +
  geom_line(alpha = 0.8, linewidth = 0.6) +
  theme_minimal() +
  labs(
    title = "4-Week Rolling Volatility (Standard Deviation)",
    subtitle = "Spikes in rolling SD suggest potential change-points in covariance matrix",
    x = "Date",
    y = "Rolling Std Dev",
    color = "Industry"
  ) +
  theme(
    plot.title = element_text(size = 13, face = "plain", hjust = 0),
    plot.subtitle = element_text(size = 11, hjust = 0),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p2)


# --- PLOT 3: Sector Correlation Matrix Heatmap --------------------------------
p3_scatter <- ggpairs(
  as.data.frame(Y),
  lower = list(continuous = wrap("points", alpha = 0.3, size = 0.8, color = "steelblue")),
  diag  = list(continuous = wrap("densityDiag", fill = "lightblue", alpha = 0.5)),
  upper = list(continuous = wrap("cor", size = 4, color = "black"))
) +
  theme_minimal() +
  labs(
    title = paste0("Scatterplot & Correlation Matrix (", format(min(weekly_returns$week), "%Y"), "-", format(max(weekly_returns$week), "%Y"), ")"),
    subtitle = "Lower: Scatter plots | Diag: Densities | Upper: Correlation coefficients"
  ) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 9, face = "bold") # Etichette dei settori
  )

print(p3_scatter)


# ==============================================================================
# SAVE PLOTS TO PDF 
# ==============================================================================

# --- Save Plot 1: Log-Returns Time Series -------------------------------------
ggsave(
  filename = "figures/application/realdata_weekly_log_returns.pdf",
  plot = p1,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300
)

# --- Save Plot 2: 4-Week Rolling Volatility ----------------------------------
ggsave(
  filename = "figures/application/realdata_rolling_volatility.pdf",
  plot = p2,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300
)

# --- Save Plot 3: Scatterplot & Correlation Matrix ---------------------------
ggsave(
  filename = "figures/application/realdata_scatterplot_correlation_matrix.pdf",
  plot = p3_scatter,
  width = 12,
  height = 12,
  units = "in",
  dpi = 300
)

