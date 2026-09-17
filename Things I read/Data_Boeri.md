# Data Requirements 

This document provides instructions for obtaining all datasets required to replicate "A Wartime Labor Market: The Case of Ukraine."

Replication requires both **public** and **restricted** datasets. The former are included in the replication package, the latter are not. Researchers should independently obtain restricted datasets and place them in `/data/raw/restricted/` before running the replication.

One note on public data from the State Statistics Service of Ukraine (SSSU). These were collected from the official website (https://ukrstat.gov.ua/). However, on November 5, 2025, the SSSU launched a new site (https://stat.gov.ua/) and stopped updating the old one. The new site is gradually migrating historical data into a new database format, and not all datasets have been transferred yet. The SSSU has not announced how long the old site will remain online. Therefore, if a dataset is not available via the provided link, use the file title referenced in this document to search for it on the new site.

---

## Public Data (Included in Package)

The following datasets are included in `/data/raw/public/`:

### 1. Air Raid Alarm Data
- **File:** `official_data_en.csv`
- **Source:** Ukrainian Air Raid Sirens Dataset
- **URL:** https://github.com/Vadimkin/ukrainian-air-raid-sirens-dataset/tree/main/datasets
- **Notes:** Copyright (c) 2022 Vadym Klymenko

### 2. Population by region (by estimate) as of February 1, 2022. Average annual populations  in January 2022.
- **File:** `kn_0122_ue.xls`
- **Source:** State Statistics Service of Ukraine 
- **URL:** https://ukrstat.gov.ua/operativ/operativ2022/ds/kn/kn_0122_ue.xls

### 3. Geographic Boundaries (shapefiles for map visualization)
- **File:** `geoBoundaries-UKR-ADM1-all/` (folder)
- **Source:** geoBoundaries
- **URL:** https://www.geoboundaries.org/countryDownloads.html
- **Notes:** Download the ADM1 folder

### 4. Consumer Price Index
- **File:** `CPI_m.xlsx`
- **Source:** National Bank of Ukraine
- **URL:** https://bank.gov.ua/en/statistic/macro-indicators#1 
- **Notes:** Updated every month

### 5. Labor force survey (SSSU)
- **File:** `SSSU_lfs.xlsx`
- **Source:** State Statistics Service of Ukraine 
- **URL:** https://stat.gov.ua/en/explorer
- **Notes:** To download the data for the SSSU Labor Force Survey, select Dataset: Labor force survey (annual).

### 6. Number of persons employed of business entities by region in 2010-2024
- **File:** `kzpsg_reg_2010_2020_ue.xlsx`
- **Source:** State Statistics Service of Ukraine 
- **URL:** https://www.ukrstat.gov.ua/operativ/operativ2019/fin/pssg/kzpsg_reg_2010_2020_ue.xlsx 

### 7. Number of persons employed of business entities by type of economic activity in regions in 2014-2023											
- **File:** `kzpsg_ved_15-20.xlsx`
- **Source:** State Statistics Service of Ukraine 
- **URL:** https://ukrstat.gov.ua/operativ/operativ2021/fin/pdsg/kzpsg_ved_15-20.xlsx

### 8. Number of persons employed of enterprises by type of economic activity with a breakdown by large, medium, small and microenterprises in 2010‒2024
- **File:** `Kzp_kved_10_21.xlsx`
- **Source:** State Statistics Service of Ukraine
- **URL:** https://www.ukrstat.gov.ua/operativ/operativ2022/fin/fin_new/Kzp_kved_10_21.xlsx

### 9. Labor statistics enterprises survey (SSSU)
- **File:** `SSSU_lses.xlsx`
- **Source:** State Statistics Service of Ukraine
- **URL:** https://stat.gov.ua/en/explorer
- **Notes:** To download the data, select Dataset: Labor statistics enterprises survey.

### 10. Turnover of enterprises with a breakdown by large, medium, small and microenterprises by regions in 2010─2024											
- **File:** `orpp_roz_reg_10_20_ue.xlsx`
- **Source:** State Statistics Service of Ukraine
- **URL:** https://www.ukrstat.gov.ua/operativ/operativ2021/fin/fin_new/orpp_roz_reg_10_20_ue.xlsx 

### 11. Employee benefits expense of enterprises with a breakdown by large, medium, small and microenterprises by regions in 2010-2024											
- **File:** `vpp_roz_reg_10_20_ue.xlsx`
- **Source:** State Statistics Service of Ukraine
- **URL:** https://www.ukrstat.gov.ua/operativ/operativ2021/fin/fin_new/vpp_roz_reg_10_20_ue.xlsx

### 12. IER monthly enterprise survey, data on capacity/production volumes compared to before the invasion
- **File:** `ier-production_to_before2022.xlsx`
- **Source:** Institute for Economic Research and Policy Consulting (IER)
- **URL:** http://www.ier.com.ua/en/proekt_dilova_dumka/survey
- **Notes:** We manually collected data from monthly publications with survey results. 

### 13. IER monthly enterprise survey, data on difficulty of finding workers
- **File:** `ier_hard_to_find_workers.xlsx`
- **Source:** Institute for Economic Research and Policy Consulting (IER)
- **URL:** http://www.ier.com.ua/en/proekt_dilova_dumka/survey
- **Notes:** We manually collected data from monthly publications with survey results. 

### 14. Unemployment rate (NBU)
- **File:** `Tables_Charts_2025-Q3_en.xlsx`
- **Source:** National Bank of Ukraine
- **URL:** https://stat.gov.ua/en/explorer
- **Notes:** Sheet 3.15 in https://bank.gov.ua/admin_uploads/article/Tables_Charts_2025-Q3_en.xlsx?v=16


---

## Restricted Data (Not Included)

Researchers must obtain the following datasets independently and place them in `/data/raw/restricted/`:

### 1. Work.ua Platform Data
- **Files:** 
  - `Work.ua_data_new.xlsx` (main, up to 2024)
  - `Work.ua_data_2025.xlsx`(main, 2025 addendum)
  - `workua.2020-2021data.xlsx`(2021 first correction)
  - `workua.2020-2021data_new.xlsx`(2021 second correction)
  - `wages-missing.xlsx` (wages correction)
  - `work-ua-changed-regions2023.xlsx` (2023 correction)
- **Source:** Work.ua 
- **URL:** https://work.ua/stat/count/ (data publicly visible but not bulk-downloadable)
- **Notes:** The main data are publicly available from the Work.ua website. The platform does not provide download options, but data can be manually retrieved or collected via automated web scraping. Due to institutional restrictions, we cannot share the specific scraping script used, but the main data are fully replicable using publicly accessible information. The only exception are data corrections that we requested and obtained directly from the platform, for cases in which data were missing or there were reporting issues. Researchers interested in obtaining such data may contact Work.ua directly. 

### 2. Labor Market Assessment 2024–2025 (survey of employers)
- **File:** `SES_firms.xlsx` 
- **Source:** State Employment Service of Ukraine and Helvetas Swiss Intercooperation
- **URL:** https://www.helvetas.org/en/eastern-europe/ukraine and https://stat.gov.ua/en
- **Notes:** Additional information on the survey can be found here: https://www.helvetas.org/Publications-PDFs/Eastern-Europe-Caucasus/Ukraine/Social%20Housing%20Reform%20in%20Ukraine/Labor%20Market%20Assessment%202024-2025%20Business%20Demand,%20Challenges,%20Strategy%20For%20Impact-ENG.pdf. We got confidential access to the survey data. Researchers interested in obtaining such data may contact these institutions. The aggregated data can be accessed at https://www.dcz.gov.ua/stat/statsurvey (in Ukrainian).

### 3. Registered Unemployed in Ukraine 2024–2025 (survey of registered unemployed)
- **File:** `SES_unemployed.xlsx` 
- **Source:** State Employment Service of Ukraine and Helvetas Swiss Intercooperation
- **URL:** https://www.helvetas.org/en/eastern-europe/ukraine and https://stat.gov.ua/en
- **Notes:** Additional information on the survey can be found here: https://www.helvetas.org/Publications-PDFs/Eastern-Europe-Caucasus/Ukraine/Social%20Housing%20Reform%20in%20Ukraine/Registered%20unemployment%20in%20Ukraine%20needs,%20features,%20assessments-ENG.pdf. We got confidential access to the survey data. Researchers interested in obtaining such data may contact these institutions. The aggregated data can be accessed at https://www.dcz.gov.ua/stat/statsurvey (in Ukrainian).

### 4. Budget data (social security contributions and military allowance)
- **File:** `ssc-data.xlsx`
- **Source:** State Treasury Service of Ukraine
- **Notes:** We have had access to budget data from the State Treasury Service of Ukraine through the National Bank of Ukraine but access is generally restricted. Researchers interested in obtaining such data may contact these institutions.


---

*Last updated: January 2026*

