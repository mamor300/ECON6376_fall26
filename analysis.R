#########
# Matt Amor
# ECON 6376 - Fall 2026
# Midterm
#########

# Vintage datasets used for this project
## are available at https://github.com/mamor300/ECON6376_fall26

pacman::p_load(
  tidyverse,
  ggplot2,
  forecast,
  patchwork,
  urca,
  gt,
  readxl,
  stringr
)

# Project Start and Frequency
Begin <- c(1981,01)
Begin_diff <- c(1980,04)
Period <- 4

###################
# Avg vs Med Home Prices
###################
# Loading Data
{
  # U.S. Census Bureau and U.S. Department of Housing and Urban Development, 
  # Average Sales Price of Houses Sold for the United States [ASPUS], 
  # retrieved from FRED, Federal Reserve Bank of St. Louis; 
  # https://fred.stlouisfed.org/series/ASPUS, October 6, 2026.
  ASPUS0 <- read_csv("ASPUS.csv")
  ASPUS <- ASPUS0[73:228,]
  
  # U.S. Census Bureau and U.S. Department of Housing and Urban Development, 
  # Median Sales Price of New Houses Sold for the United States [MSPUS], 
  # retrieved from FRED, Federal Reserve Bank of St. Louis; 
  # https://fred.stlouisfed.org/series/MSPUS, October 6, 2026.
  MSPUS0 <- read_csv("MSPUS.csv")
  MSPUS <- MSPUS0[73:228,]
  
  ## Prices are standardized without a deflator
  # U.S. Bureau of Labor Statistics, 
  # Consumer Price Index for All Urban Consumers: All Items in U.S. City Average [CPIAUCSL], 
  # retrieved from FRED, Federal Reserve Bank of St. Louis; 
  # https://fred.stlouisfed.org/series/CPIAUCSL, October 6, 2026.
  CPI0 <- read_csv("CPIAUCSL.csv")
  CPI <- CPI0[137:292,]
  
  # Combining ASPUS and MSPUS into one data frame
  SPUS <- cbind(ASPUS,MSPUS$MSPUS,CPI$CPIAUCSL) |>
    rename(
      AverageHomePrice = ASPUS,
      MedianHomePrice = "MSPUS$MSPUS",
      CPI = `CPI$CPIAUCSL`
    ) |>
    mutate(
      AverageHomePrice = round(AverageHomePrice/CPI*332.985, digits = 2),
      MedianHomePrice = round(MedianHomePrice/CPI*332.985, digits = 2),
      Spread = AverageHomePrice - MedianHomePrice,
      # Change_avg = (AverageHomePrice - lag(AverageHomePrice)) / lag(AverageHomePrice),
      # Change_med = (MedianHomePrice - lag(MedianHomePrice)) / lag(MedianHomePrice),
      # Change_diff = (Spread - lag(Spread)) / lag(Spread),
      Spread_Med = (Spread/MedianHomePrice)*100
    )
}

# Average Home Prices
{
  average <- ts(SPUS$AverageHomePrice, start = Begin, frequency = Period)
  
  # Plotting the average
  p1_average <- autoplot(average) + 
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 0.4,
                linetype = "dashed")+
    theme_bw() +
    labs(title = expression(paste("FIGURE 1 - Average Home Prices 1981 - 2019")),
         subtitle = "", 
         y = "Price in 2026 USD")
  p2_average <- ggAcf(average, lag.max = 24) + theme_bw() + ggtitle("ACF")
  p3_average <- ggPacf(average, lag.max = 24) + theme_bw() + ggtitle("PACF")
  p_average <- p1_average / (p2_average + p3_average)
  p_average
  # ggsave("Figure01.png", p_average)
  
  # Testing average Stationarity
  # Strong Evidence - Non-Stationary
  dftest_average <- ur.df(average)
  summary(dftest_average)
  kpss_average <- ur.kpss(average)
  summary(kpss_average)
  
  # Detrending average
  time_average <- 1:length(average)
  trend_average <- lm(as.numeric(average)~time_average)
  de_average <- ts(residuals(trend_average),start = Begin, frequency = Period)
  
  # Plotting detrended average
  p1_de_average <- autoplot(de_average) + 
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 0.4,
                linetype = "dashed")+
    theme_bw() +
    labs(title = expression(paste("FIGURE 2 - Detreneded Avg. Home Prices")),
         subtitle = "", 
         y = "Price (2026 USD)")
  p2_de_average <- ggAcf(de_average, lag.max = 24) + theme_bw() + ggtitle("ACF")
  p3_de_average <- ggPacf(de_average, lag.max = 24) + theme_bw() + ggtitle("PACF")
  p_de_average <- p1_de_average / (p2_de_average + p3_de_average)
  p_de_average
  # ggsave("Figure02.png", p_de_average)
  
  # Testing Detrended average Stationarity
  # Strong Evidence - df:Stationary ; kpss:Stationary
  dftest_de_average <- ur.df(de_average)
  summary(dftest_de_average)
  kpss_de_average <- ur.kpss(de_average)
  summary(kpss_de_average)
  
  
  # Differencing Average Home Price
  diff_average <- diff(average)
  
  # Plotting detrended average
  p1_diff_average <- autoplot(diff_average) + 
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 0.4,
                linetype = "dashed")+
    theme_bw() +
    labs(title = expression(paste("FIGURE 3 - Differenced Avg. Home Prices")),
         subtitle = "", 
         y = "Price (2026 USD)")
  p2_diff_average <- ggAcf(diff_average, lag.max = 24) + theme_bw() + ggtitle("ACF")
  p3_diff_average <- ggPacf(diff_average, lag.max = 24) + theme_bw() + ggtitle("PACF")
  p_diff_average <- p1_diff_average / (p2_diff_average + p3_diff_average)
  p_diff_average
  # ggsave("Figure03.png", p_diff_average)
  
  # Testing Detrended average Stationarity
  # Strong Evidence - df:Stationary ; kpss:Stationary
  dftest_diff_average <- ur.df(diff_average)
  summary(dftest_diff_average)
  kpss_diff_average <- ur.kpss(diff_average)
  summary(kpss_diff_average)
  
  
  # Creating a Table for DF and KPSS tests
  unit_root_rows <- function(y, name, df_type = "drift", kpss_type = "mu") {
    df <- ur.df(y, type = df_type)
    kp <- ur.kpss(y, type = kpss_type)
    
    tibble(
      series    = name,
      test      = c("ADF", "KPSS"),
      null      = c("Unit root", "Stationary"),
      statistic = c(df@teststat[1, 1], kp@teststat),
      cv1       = c(df@cval[1, "1pct"],  kp@cval[1, "1pct"]),
      cv5       = c(df@cval[1, "5pct"],  kp@cval[1, "5pct"]),
      cv10      = c(df@cval[1, "10pct"], kp@cval[1, "10pct"])
    ) |>
      mutate(decision = case_when(
        test == "ADF"  & statistic < cv5 ~ "Reject null",
        test == "KPSS" & statistic > cv5 ~ "Reject null",
        TRUE                             ~ "Fail to reject"))
  }
  
  tab_avg <- bind_rows(
    unit_root_rows(average, "Average Home Prices"),
    unit_root_rows(de_average, name = "Detrended Avg. Home Prices", df_type = "none"),
    unit_root_rows(diff_average, "Differenced Avg. Home Prices")
  )
  
  Table2 <- tab_avg |>
    gt(groupname_col = "series") |>
    tab_header(title = "TABLE 2 - Unit-root and stationarity tests") |>
    tab_spanner(label = "Critical values", columns = c(cv1, cv5, cv10)) |>
    cols_label(test = "Test", null = "Null", statistic = "Statistic",
               cv1 = "1%", cv5 = "5%", cv10 = "10%", decision = "Decision at 5%") |>
    fmt_number(columns = c(statistic, cv1, cv5, cv10), decimals = 3)|>
    cols_width(
      test ~ px(80),
      null ~ px(80),
      statistic ~ px(110),
      cv1 ~ px(80),
      cv5 ~ px(80),
      cv10 ~ px(80),
      decision ~ px(130)
    ) |>
    cols_align(align = "center")
  Table2
  # gtsave(Table2, "Table02.png", expand = 50)
  
}

# Median Home Prices
{
  median <- ts(SPUS$MedianHomePrice, start = Begin, frequency = Period)
  
  # Plotting the median
  p1_median <- autoplot(median) + 
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 0.4,
                linetype = "dashed")+
    theme_bw() +
    labs(title = expression(paste("FIGURE 4 - Median Home Prices 1981 - 2019")),
         subtitle = "", 
         y = "Price in 2026 USD")
  p2_median <- ggAcf(median, lag.max = 24) + theme_bw() + ggtitle("ACF")
  p3_median <- ggPacf(median, lag.max = 24) + theme_bw() + ggtitle("PACF")
  p_median <- p1_median / (p2_median + p3_median)
  p_median
  # ggsave("Figure04.png", p_median)
  
  # Testing median Stationarity
  # Strong Evidence - Non-Stationary
  dftest_median <- ur.df(median)
  summary(dftest_median)
  kpss_median <- ur.kpss(median)
  summary(kpss_median)
  
  # Detrending median
  time_median <- 1:length(median)
  trend_median <- lm(as.numeric(median)~time_median)
  de_median <- ts(residuals(trend_median),start = Begin, frequency = Period)
  
  # Plotting detrended median
  p1_de_median <- autoplot(de_median) + 
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 0.4,
                linetype = "dashed")+
    theme_bw() +
    labs(title = expression(paste("FIGURE 5 - Detrended Med. Home Prices")),
         subtitle = "", 
         y = "Price (2026 USD)")
  p2_de_median <- ggAcf(de_median, lag.max = 24) + theme_bw() + ggtitle("ACF")
  p3_de_median <- ggPacf(de_median, lag.max = 24) + theme_bw() + ggtitle("PACF")
  p_de_median <- p1_de_median / (p2_de_median + p3_de_median)
  p_de_median
  # ggsave("Figure05.png", p_de_median)
  
  # Testing Detrended median Stationarity
  # Strong Evidence - df:Stationary ; kpss:Stationary
  dftest_de_median <- ur.df(de_median)
  summary(dftest_de_median)
  kpss_de_median <- ur.kpss(de_median)
  summary(kpss_de_median)
  
  # Differencing median Home Price
  diff_median <- diff(median)
  
  # Plotting differenced median
  p1_diff_median <- autoplot(diff_median) + 
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 0.4,
                linetype = "dashed")+
    theme_bw() +
    labs(title = expression(paste("FIGURE 6 - Differenced Median Home Prices 1981 - 2019")),
         subtitle = "", 
         y = "Prices (2026 USD)")
  p2_diff_median <- ggAcf(diff_median, lag.max = 24) + theme_bw() + ggtitle("ACF")
  p3_diff_median <- ggPacf(diff_median, lag.max = 24) + theme_bw() + ggtitle("PACF")
  p_diff_median <- p1_diff_median / (p2_diff_median + p3_diff_median)
  p_diff_median
  # ggsave("Figure06.png", p_diff_median)
  
  # Testing Detrended median Stationarity
  # Strong Evidence - df:Stationary ; kpss:Stationary
  dftest_diff_median <- ur.df(diff_median)
  summary(dftest_diff_median)
  kpss_diff_median <- ur.kpss(diff_median)
  summary(kpss_diff_median)
  
  # Generating a table for median home prices
  tab_med <- bind_rows(
    unit_root_rows(median, "Median Home Prices"),
    unit_root_rows(de_median, name = "Detrended Med. Home Prices", df_type = "none"),
    unit_root_rows(diff_median, "Differenced Med. Home Prices")
  )
  
  t_tab_med <- tab_med |>
    gt(groupname_col = "series") |>
    tab_header(title = "TABLE 3 - ADF and KPSS Test Results of Median Home Prices") |>
    tab_spanner(label = "Critical values", columns = c(cv1, cv5, cv10)) |>
    cols_label(
      test = "Test", 
      null = "Null", 
      statistic = "Statistic",
      cv1 = "1%", 
      cv5 = "5%", 
      cv10 = "10%", 
      decision = "Decision at 5%") |>
    fmt_number(columns = c(statistic, cv1, cv5, cv10), decimals = 3) |>
    cols_width(
      test ~ px(80),
      null ~ px(80),
      statistic ~ px(110),
      cv1 ~ px(80),
      cv5 ~ px(80),
      cv10 ~ px(80),
      decision ~ px(130)
    ) |>
    cols_align(align = "center")
  t_tab_med
  # gtsave(t_tab_med, "Table03.png", expand = 50)
}

# Avg-Med Spread as a Percent of the Median Price
{
  spread_med <- ts(SPUS$Spread_Med, start = Begin, frequency = Period)
  
  # Plotting the spread_med
  p1_spread_med <- autoplot(spread_med) + 
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 0.4,
                linetype = "dashed")+
    theme_bw() +
    labs(title = expression(paste("FIGURE 7 - Skew of Home Price Distribution 1981 - 2019")),
         subtitle = "", 
         y = "% of Median Price")
  p2_spread_med <- ggAcf(spread_med, lag.max = 24) + theme_bw() + ggtitle("ACF")
  p3_spread_med <- ggPacf(spread_med, lag.max = 24) + theme_bw() + ggtitle("PACF")
  p_spread_med <- p1_spread_med / (p2_spread_med + p3_spread_med)
  p_spread_med
  # ggsave("Figure07.png", p_spread_med)
  
  # Testing spread_med Stationarity
  # Mixed Evidence - df:Non-Stationary ; kpss:stationary
  dftest_spread_med <- ur.df(spread_med)
  summary(dftest_spread_med)
  kpss_spread_med <- ur.kpss(spread_med)
  summary(kpss_spread_med)
  
  # Detrending spread_med
  time_spread_med <- 1:length(spread_med)
  trend_spread_med <- lm(as.numeric(spread_med)~time_spread_med)
  de_spread_med <- ts(residuals(trend_spread_med),start = Begin, frequency = Period)
  
  # Plotting detrended spread_med
  p1_de_spread_med <- autoplot(de_spread_med) + 
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 0.4,
                linetype = "dashed")+
    theme_bw() +
    labs(title = expression(paste("FIGURE 8 - Detrended Skew of Home Prices 1981-2019")),
         subtitle = "", 
         y = "% of Median")
  p2_de_spread_med <- ggAcf(de_spread_med, lag.max = 24) + theme_bw() + ggtitle("ACF")
  p3_de_spread_med <- ggPacf(de_spread_med, lag.max = 24) + theme_bw() + ggtitle("PACF")
  p_de_spread_med <- p1_de_spread_med / (p2_de_spread_med + p3_de_spread_med)
  p_de_spread_med
  # ggsave("Figure08.png", p_de_spread_med)
  
  # Testing Detrended spread_med Stationarity
  # Strong Evidence - df:Stationary ; kpss:Stationary
  dftest_de_spread_med <- ur.df(de_spread_med)
  summary(dftest_de_spread_med)
  kpss_de_spread_med <- ur.kpss(de_spread_med)
  summary(kpss_de_spread_med)
  
  # Generating a table for home prices skew
  tab_spread <- bind_rows(
    unit_root_rows(spread_med, "Skew of Home Prices", df_type = "drift"),
    unit_root_rows(de_spread_med, name = "Detrended Skew of Home Prices", df_type = "none")
  )
  
  t_tab_spread <- tab_spread |>
    gt(groupname_col = "series") |>
    tab_header(title = "TABLE 4 - ADF and KPSS Test Results of Skew of Home Prices") |>
    tab_spanner(label = "Critical values", columns = c(cv1, cv5, cv10)) |>
    cols_label(test = "Test", null = "Null", statistic = "Statistic",
               cv1 = "1%", cv5 = "5%", cv10 = "10%", decision = "Decision at 5%") |>
    fmt_number(columns = c(statistic, cv1, cv5, cv10), decimals = 3)|>
    cols_width(
      test ~ px(80),
      null ~ px(80),
      statistic ~ px(110),
      cv1 ~ px(80),
      cv5 ~ px(80),
      cv10 ~ px(80),
      decision ~ px(130)
    ) |>
    cols_align(align = "center")
  t_tab_spread
  # gtsave(t_tab_spread, "Table04.png", expand = 50)
  
}

###################
# Income Inequality
###################

# Loading Realtime Inequality Data
{
  # Monthly Inequality data from realtimeinequality.org
  Inequality <- read_xlsx("Full Database.xlsx")
  # Measuring inequality by share of income
  ShareInc <- Inequality |>
    mutate(unit = as.factor(unit)) |>
    filter(unit %in% "adult_households") |>
    select(year, quarter, month, group, factor_income_share) |>
    pivot_wider(
      names_from = "group",
      values_from = "factor_income_share"
    )
  ShareIncQ0 <- ShareInc |>
    mutate(date = paste0(year,"_Q",quarter))|>
    select(date, everything(), -c(year, quarter, month))|>
    group_by(date) |>
    summarise(across(.cols = everything(),\(x) mean(x,na.rm = TRUE)*100), .groups = "drop")
}  

# Raw Inequality Data to 1981 - 2019
{
  ShareIncQ <- ShareIncQ0[20:176,]
  
  ## Creating TS of Top 10%
  Top10 <- ts(ShareIncQ$"Top 10%",start = Begin_diff, frequency = Period)
  
  ## Plotting Top 10
  p1_Top10 <- autoplot(Top10) + 
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 0.4,
                linetype = "dashed")+
    theme_bw() +
    labs(title = expression(paste("FIGIURE 9 - Income Share of the Top 10% of Income Earners")),
         subtitle = "", 
         y = "% of Total Income")
  p2_Top10 <- ggAcf(Top10, lag.max = 24) + theme_bw() + ggtitle("ACF")
  p3_Top10 <- ggPacf(Top10, lag.max = 24) + theme_bw() + ggtitle("PACF")
  p_Top10 <- p1_Top10 / (p2_Top10 + p3_Top10)
  p_Top10
  # ggsave("Figure09.png", p_Top10)
  
  # Testing Stationarity of top 10% 
  # Strong Evidence - Non-Stationary
  dftest_Top10 <- ur.df(Top10)
  summary(dftest_Top10)
  kpss_Top10 <- ur.kpss(Top10)
  summary(kpss_Top10)
}

# Detrending Inequality
{
time_Top10 <- 1:length(Top10)
trend_Top10 <- lm(as.numeric(Top10)~time_Top10)
de_Top10 <- ts(residuals(trend_Top10), start = Begin, frequency = Period)

# Plotting Detrended Income Inequality
p1_de_Top10 <- autoplot(de_Top10) + 
  geom_smooth(method = "lm",
              se = FALSE,
              linewidth = 0.4,
              linetype = "dashed")+
  theme_bw() +
  labs(title = expression(paste("FIGURE 10 - Detrended Income Share of the Top 10% of Income Earners")),
       subtitle = "", 
       y = "% of Income")
p2_de_Top10 <- ggAcf(de_Top10, lag.max = 24) + theme_bw() + ggtitle("ACF")
p3_de_Top10 <- ggPacf(de_Top10, lag.max = 24) + theme_bw() + ggtitle("PACF")
p_de_Top10 <- p1_de_Top10 / (p2_de_Top10 + p3_de_Top10)
p_de_Top10
# ggsave("Figure10.png", p_de_Top10)

# Testing Stationarity of Detrended Inequality 
# Mixed Evidence - df: stationary at 5% ; kpss: stationary at 5%
dftest_de_Top10 <- ur.df(de_Top10)
summary(dftest_de_Top10)
kpss_de_Top10 <- ur.kpss(de_Top10)
summary(kpss_de_Top10)
}

# Differencing Top10
{
diff_Top10 <- diff(Top10)
p1_diff_Top10 <- autoplot(diff_Top10) + 
  geom_smooth(method = "lm",
              se = FALSE,
              linewidth = 0.4,
              linetype = "dashed")+
  theme_bw() +
  labs(title = expression(paste("FIGURE 11 - First Differencing of Inequality")),
       subtitle = "", 
       y = "% of Income")
p2_diff_Top10 <- ggAcf(diff_Top10, lag.max = 24) + theme_bw() + ggtitle("ACF")
p3_diff_Top10 <- ggPacf(diff_Top10, lag.max = 24) + theme_bw() + ggtitle("PACF")
p_diff_Top10 <- p1_diff_Top10 / (p2_diff_Top10 + p3_diff_Top10)
p_diff_Top10
# ggsave("Figure11.png", p_diff_Top10)

# Testing Stationarity of Differenced Top10
# Strong Evidence - df & kpss support stationarity
dftest_diff_Top10 <- ur.df(diff_Top10)
summary(dftest_diff_Top10)
kpss_diff_Top10 <- ur.kpss(diff_Top10)
summary(kpss_diff_Top10)
}

# Generating a table for Inequality
{
tab_ineq <- bind_rows(
  unit_root_rows(Top10, "Inequality", df_type = "drift"),
  unit_root_rows(de_Top10, name = "Detrended Inequality", df_type = "none"),
  unit_root_rows(diff_Top10, name = "Differenced Inequality")
)

t_tab_ineq <- tab_ineq |>
  gt(groupname_col = "series") |>
  tab_header(title = "TABLE 5 - ADF and KPSS Test Results of Inequality") |>
  tab_spanner(label = "Critical values", columns = c(cv1, cv5, cv10)) |>
  cols_label(test = "Test", null = "Null", statistic = "Statistic",
             cv1 = "1%", cv5 = "5%", cv10 = "10%", decision = "Decision at 5%") |>
  fmt_number(columns = c(statistic, cv1, cv5, cv10), decimals = 3)|>
  cols_width(
    test ~ px(80),
    null ~ px(80),
    statistic ~ px(110),
    cv1 ~ px(80),
    cv5 ~ px(80),
    cv10 ~ px(80),
    decision ~ px(130)
  ) |>
  cols_align(align = "center")
t_tab_ineq
# gtsave(t_tab_ineq, "Table05.png", expand = 50)
}

###################
# Housing Inventory
###################

# Loading Census Inventory Data
{
  Inventory_rough <- read_excel(
    "histtab8.xlsx", 
    col_names = FALSE,
    col_types = "text")
  ForSale <- which(str_detect(Inventory_rough[[1]], "For sale only") & !is.na(Inventory_rough[[1]]))
  get_block <- function(i) {
    # quarter-header row: nearest row above containing "Qtr"
    j <- i - 1
    while (!any(str_detect(unlist(Inventory_rough[j, 2:5]), "Qtr"), na.rm = TRUE)) j <- j - 1
    # year row: nearest non-empty row above that
    k <- j - 1
    while (all(is.na(Inventory_rough[k, ]))) k <- k - 1
    yr_cell <- na.omit(unlist(Inventory_rough[k, ]))[1]       # e.g. "19651", "1979r1", "2009"
    
    tibble(
      year          = as.integer(str_sub(yr_cell, 1, 4)),  # drops footnote suffix
      revised       = str_detect(yr_cell, "r\\d"),
      quarter       = 1:4,
      for_sale_only = as.numeric(unlist(Inventory_rough[i, 2:5]))      # thousands of units
    )
  }
  
  Inventory0 <- bind_rows(lapply(ForSale, get_block)) |>
    filter(!is.na(for_sale_only)) |>
    mutate(date = as.Date(sprintf("%d-%02d-01", year, 3 * quarter - 2))) |>
    arrange(year, quarter, desc(revised)) |>
    distinct(year, quarter, .keep_all = TRUE) |>
    select(date, for_sale_only)
  
  Inventory <- Inventory0[64:220,]
}

# Generating TS for Housing inventory
{
  inventory <- ts(Inventory$for_sale_only, start = Begin_diff, frequency = Period)
  
  # Plotting the inventory
  p1_inventory <- autoplot(inventory) + 
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 0.4,
                linetype = "dashed")+
    theme_bw() +
    labs(title = expression(paste("FIGURE 12 - Housing Inventory 1981 - 2019")),
         subtitle = "", 
         y = "thousands of units")
  p2_inventory <- ggAcf(inventory, lag.max = 24) + theme_bw() + ggtitle("ACF")
  p3_inventory <- ggPacf(inventory, lag.max = 24) + theme_bw() + ggtitle("PACF")
  p_inventory <- p1_inventory / (p2_inventory + p3_inventory)
  p_inventory
  # ggsave("Figure12.png", p_inventory)
  
  # Testing inventory Stationarity
  # Strong Evidence - Non-Stationary
  dftest_inventory <- ur.df(inventory)
  summary(dftest_inventory)
  kpss_inventory <- ur.kpss(inventory)
  summary(kpss_inventory)
  
  # Differencing inventory
  diff_inventory <- diff(inventory)
  
  # Plotting differenced inventory
  p1_diff_inventory <- autoplot(diff_inventory) + 
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 0.4,
                linetype = "dashed")+
    theme_bw() +
    labs(title = expression(paste("FIGURE 13 - Differenced Housing Inventory")),
         subtitle = "", 
         y = "thousands of Units")
  p2_diff_inventory <- ggAcf(diff_inventory, lag.max = 24) + theme_bw() + ggtitle("ACF")
  p3_diff_inventory <- ggPacf(diff_inventory, lag.max = 24) + theme_bw() + ggtitle("PACF")
  p_diff_inventory <- p1_diff_inventory / (p2_diff_inventory + p3_diff_inventory)
  p_diff_inventory
  # ggsave("Figure13.png", p_diff_inventory)
  
  # Testing Differenced inventory Stationarity
  # Strong Evidence - df:Stationary ; kpss:Stationary
  dftest_diff_inventory <- ur.df(diff_inventory)
  summary(dftest_diff_inventory)
  kpss_diff_inventory <- ur.kpss(diff_inventory)
  summary(kpss_diff_inventory)
  
  # Generating Detrended Inventory TS
  time_inventory <- 1:length(inventory)
  trend_inventory <- lm(as.numeric(inventory)~time_inventory)
  de_inventory <- ts(residuals(trend_inventory),start = Begin_diff, frequency = Period)
  
  # Plotting detrended inventory
  p1_de_inventory <- autoplot(de_inventory) + 
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 0.4,
                linetype = "dashed")+
    theme_bw() +
    labs(title = expression(paste("FIGURE 14 - Detrended Housing Inventory")),
         subtitle = "", 
         y = "thousands of Units")
  p2_de_inventory <- ggAcf(de_inventory, lag.max = 24) + theme_bw() + ggtitle("ACF")
  p3_de_inventory <- ggPacf(de_inventory, lag.max = 24) + theme_bw() + ggtitle("PACF")
  p_de_inventory <- p1_de_inventory / (p2_de_inventory + p3_de_inventory)
  p_de_inventory
  # ggsave("Figure14.png", p_de_inventory)
  
  # Generating a table for Inventory
  tab_inv <- bind_rows(
    unit_root_rows(inventory, "Housing Inventory"),
    unit_root_rows(de_inventory, name = "Detrended Inventory", df_type = "none"),
    unit_root_rows(diff_inventory, name = "Differenced Inventory")
  )
  
  t_tab_inv <- tab_inv |>
    gt(groupname_col = "series") |>
    tab_header(title = "TABLE 6 - ADF and KPSS Test Results of Inventory") |>
    tab_spanner(label = "Critical values", columns = c(cv1, cv5, cv10)) |>
    cols_label(test = "Test", null = "Null", statistic = "Statistic",
               cv1 = "1%", cv5 = "5%", cv10 = "10%", decision = "Decision at 5%") |>
    fmt_number(columns = c(statistic, cv1, cv5, cv10), decimals = 3)|>
    cols_width(
      test ~ px(80),
      null ~ px(80),
      statistic ~ px(110),
      cv1 ~ px(80),
      cv5 ~ px(80),
      cv10 ~ px(80),
      decision ~ px(130)
    ) |>
    cols_align(align = "center")
  t_tab_inv
  # gtsave(t_tab_inv, "Table06.png", expand = 50)
  
}

###################
# Data Table
###################

  # Generating information to populate data summary table 
{
  series_list <- list(
    "Average Home Prices" = ASPUS,
    "Median Home Prices" = MSPUS,
    "Average - Median Home Prices" = SPUS,
    "Top 10% Income Earners" = ShareIncQ,
    "Housing Inventory" = Inventory
  )
  # Individual series ASPUS and MSPUS were retrieved from FRED,
  ## but originated from US Census Bureau [https://www.census.gov/construction/nrs/data/series.html]
  source_list <- tibble(
    source = c(
      "ASPUS - US Census/HUD via FRED",
      "MSPUS - US Census/HUD via FRED",
      "Constructed from ASPUS and MSPUS",
      "Realtime Inequality",
      "CPS/HVS: Table 8 - US Census Bureau" 
    ),
    url = c(
      "https://fred.stlouisfed.org/series/ASPUS",
      "https://fred.stlouisfed.org/series/MSPUS",
      "https://fred.stlouisfed.org/series/ASPUS",
      "realtimeinequality.org",
      "https://www.census.gov/housing/hvs/data/histtabs.html"
    )
  )
  N_list <- list(
    length(ASPUS$ASPUS),
    length(MSPUS$MSPUS),
    length(SPUS$Spread_Med),
    length(ShareIncQ$`Top 10%`),
    length(Inventory$for_sale_only)
  )
  Units_list <- list(
    "2017 USD",
    "2017 USD",
    "Share of Median Price (%)",
    "Share of Income Earned (%)",
    "Housing Units (thousands)"
  )
  Seasonal_list <- list(
    "Yes",
    "Yes",
    "Yes",
    "No",
    "Yes"
  )
  Trans_list <- list(
    "Detrended",
    "Detrended",
    "Detrended",
    "Diff. (1)",
    "Diff. (1)"
  )
  Acc_list <- list(
    "10/06/2026",
    "10/06/2026",
    "10/06/2026",
    "10/07/2026",
    "10/21/2026"
  )
}
  # Building Data Summary
  data_table <- tibble(
    Series = names(series_list),
    `Source & Series ID` = source_list$url,
    Frequency = "Quarterly",
    Span = "1981Q1 - 2019Q4",
    N = N_list,
    Units = Units_list,
    SA = Seasonal_list,
    Transformation = Trans_list,
    Accessed = Acc_list
  )
  
  # Generating Data Summary Table
  Table1 <- data_table |>
    gt() |>
    fmt_url(columns = `Source & Series ID`, label = source_list$source ) |>
    tab_header(
      title = "TABLE 1 - Data Sources and Summary"
    ) 
  Table1
  # gtsave(filename = "Table01.html",Table1)
  
  # Converting TS objects into data frames
  data_SPUS <- as.data.frame(de_spread_med) |>
    rename( "Avg-Med Spread" = x)
  data_Top10 <- as.data.frame(diff_Top10) |>
    rename( "Top 10% Income" = x)
  data_Inv <- as.data.frame(diff_inventory) |>
    rename( "Housing Inventory" = x)
  data_block <- as.data.frame(c(
    rep("train", times = 132),
    rep("weights", times = 12),
    rep("holdout", times = 12)))|>
    rename(block = 1)
  
  # Collating series into one csv
  data_csv <- ASPUS |>
    cbind(
      MSPUS$MSPUS,
      data_SPUS,
      data_Top10,
      data_Inv,
      data_block
    ) |>
    rename(
      `Average Home Price` = ASPUS,
      `Median Home Price` = `MSPUS$MSPUS`)
  
  # Saving csv
  # write_csv(data_csv, "data.csv")
###################
# Section 5
###################

spread <- de_spread_med
plot <- p_de_spread_med
# rm(list = setdiff(ls(), c("spread", "plot")))

# Dividing TS into Train--Weights--Holdout
T_origin <- 132

t_dates <- 1:T_origin           # train window:   1981Q1 - 2013Q4
w_dates <- T_origin + 1:12      # weights window: 2014Q1 - 2016Q4
h_dates <- T_origin + 13:24     # holdout:        2017Q1 - 2019Q4

spread_tr <- ts(spread[t_dates], start = c(1981,01),frequency = 4)
spread_wt <- ts(spread[w_dates], start = c(2014,01),frequency = 4)
spread_ho <- ts(spread[h_dates], start = c(2017,01),frequency = 4)

# Testing ARIMA (2,0,0) model and alts
fit_100 <- Arima(spread_tr, order = c(1,0,0), include.mean = FALSE)
fit_200 <- Arima(spread_tr, order = c(2,0,0), include.mean = FALSE)
fit_300 <- Arima(spread_tr, order = c(3,0,0), include.mean = FALSE)
fit_101 <- Arima(spread_tr, order = c(1,0,1), include.mean = FALSE)
fit_201 <- Arima(spread_tr, order = c(2,0,1), include.mean = FALSE)

# Without "stationary = TRUE" auto.arima wants to pick ARIMA(0,1,1) w/ AIC=-622
## It tests for stationarity with KPSS (testing with ADF resolves the issue). 
## This TS (spread_tr; not spread_med tested in section 4) fails stationarity at 5% w/ KPSS, 
## so auto.arima does not entertain any model without a d term.
## The function only agrees with actual best fit when it's forced to look
fit_auto <- auto.arima(
  spread_tr,
  # stationary = TRUE,
  seasonal = FALSE,
  stepwise = FALSE,
  allowmean = FALSE,
  trace = TRUE,
  # test = "adf",
  ic = "aic"
)

fit_list <- list(
  "ARIMA(1,0,0)" = fit_100,
  "ARIMA(2,0,0)" = fit_200,
  "ARIMA(3,0,0)" = fit_300,
  "ARIMA(1,0,1)" = fit_101,
  "ARIMA(2,0,1)" = fit_201
)

fit_table <- tibble(
  Model = names(fit_list),
  k = vapply(fit_list,function(f) length(f$coef) + 1L, integer(1)),
  AIC = vapply(fit_list, function(f) f$aic, numeric(1)),
  AICc = vapply(fit_list, function(f) f$aicc, numeric(1)),
  BIC = vapply(fit_list, function(f) f$bic, numeric(1))
)

t_fit <- fit_table |>
  gt() |>
  tab_header(
    title = "TABLE 7 - Comparison of Home Price Skew Models"
  ) |>
  cols_width(
    Model~px(150),
    k~px(50),
    AIC~px(130),
    AICc~px(130),
    BIC~px(130)
  )
t_fit
# gtsave(filename = "Table07.png",t_fit, expand = 20)

###################
# Section 6
###################

# Testing #1 model ARIMA(1,0,1)
checkresiduals(fit_101, lag = 10)

# Testing #2 model - ARIMA(2,0,1)
checkresiduals(fit_201, lag = 10)

# Testing #3 model - ARIMA(3,0,0)
checkresiduals(fit_300, lag = 10)

###################
# Section 7 
###################

fits <- list(
  arima_1 = fit_101, 
  arima_2 = fit_201, 
  arima_3 = fit_300
)

upto <- function(x, n) {
  # First n observations of a ts, KEEPING its frequency and dates.
  window(x, end = time(x)[n])
}
kavg <- function(x, k, dates) {
  x <- as.numeric(x)
  sapply(dates, function(t) {
    if (is.infinite(k)) mean(x[1:(t - 1)]) else mean(x[(t - k):(t - 1)])
  })
}
# Function that generates the 6 tracks
tracks_on <- function(x, fits, dates, k = 3) {
  xs <- upto(x, max(dates))
  os <- function(fit) as.numeric(fitted(Arima(xs, model = fit)))[dates]
  cbind(
    arima_1 = os(fits[[1]]),
    arima_2 = os(fits[[2]]),
    arima_3 = os(fits[[3]]),
    locf    = kavg(x, 1,   dates),
    avg_k   = kavg(x, k,   dates),
    meanf   = kavg(x, Inf, dates)
  )
}

# Defining the Dates
d0 <- as.Date(sprintf("%d-%02d-01", start(spread)[1], start(spread)[2]))
date_vec <- seq(d0, by = "quarter", length.out = length(spread))

# Running Predictions for each track using the weights sample
F_w <- tracks_on(spread, fits, w_dates) 
all.equal(F_w[, "locf"], as.numeric(spread)[w_dates - 1])
all.equal(unname(F_w[1, "avg_k"]), mean(as.numeric(spread)[(w_dates[1] - 3):(w_dates[1] - 1)]))
y_w <- as.numeric(spread)[w_dates]
y_t <- as.numeric(spread)
forecasts_csv <- data.frame(date = date_vec[w_dates], round(F_w, 3))
# write_csv(forecasts_csv, "forecasts.csv")

# Plotting 6 track predictions through weights window 
## Function will work with holdouts by changing "dates_all"
dates_all <- c(w_dates)
tr_df <- data.frame(
  time  = rep(as.numeric(time(spread))[dates_all], 6),
  value = c(rbind(F_w)),
  track = rep(colnames(F_w), each = 12))
p_forecasts <- ggplot() +
  geom_line(data = subset(tr_df), #, track != "meanf"
            aes(time, value, colour = track), linewidth = 0.6) +
  geom_line(aes(x = as.numeric(time(spread))[dates_all],
                y = as.numeric(spread)[dates_all]),
            colour = "black", linewidth = 1.1) +
  geom_vline(xintercept = mean(as.numeric(time(spread))[T_origin + 1]),
             linetype = "dashed") +
  # annotate("text", x = 2018.4, y = 4.25, label = "weights window") +
  # annotate("text", x = 2019.4, y = 4.25, label = "holdout") +
  labs(
    x = "Year", 
    y = "% of Median Price",
    title = "FIGURE 15 - Actuals (black) and Six One-Step Tracks",
    # subtitle = "(meanf is relatively constant at 0.0039 - 0.0025)",
    color = "Track"
  ) +
  theme_bw() + theme(legend.position = "bottom")
p_forecasts
# ggsave("Figure15.png",p_forecasts)

###################
# Section 8
###################

# Creating a data frame containing the residuals from the models
err <- forecasts_csv |>
  mutate(
    actuals = spread[w_dates],
    arima_1 = actuals - arima_1,
    arima_2 = actuals - arima_2,
    arima_3 = actuals - arima_3,
    avg_k = actuals - avg_k,
    locf = actuals - locf,
    meanf = actuals - meanf
  ) |>
  select(-date, -actuals)

# Building a data frame reporting on the error stats
err_stats <- tibble(
  Model = colnames(err),
  MSE = colMeans(err^2),
  MAE = abs(colMeans(err))
)|>
  mutate(
    RMSE = sqrt(MSE)
  ) |>
  select(Model, MSE, RMSE, MAE)

# Creating a table from the error stats
t_err <- err_stats |>
  gt()|>
  tab_header(
    title = "TABLE 8 - Comparison of Models' Loss Functions"
  ) |>
  fmt_number(c(MSE, MAE, RMSE), decimals = 5) |>
  cols_width(
    Model~px(150),
    MSE~px(130),
    MAE~px(130),
    RMSE~px(130)
  ) 
t_err
# gtsave(t_err, "Table08.png", expand = 20)

