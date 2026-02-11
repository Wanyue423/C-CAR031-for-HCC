Project Name (DSP_Code)

1. System Requirements

Before running this project, please ensure that your system meets the following requirements:
**Operating System**:
  macOS Sequoia **15.3.1**
**R Version**:
  R **4.3.1** (2023-06-16)
> It is strongly recommended to use the specified R version to avoid potential package compatibility issues.

2. Installation Guide

 2.1 macOS Installation
Download and install macOS Sequoia from the official Apple website:
[https://www.apple.com.cn/macos/macos-sequoia/](https://www.apple.com.cn/macos/macos-sequoia/)

 2.2 R Installation
It is recommended to download R from the Tsinghua University CRAN mirror:
[https://mirrors.tuna.tsinghua.edu.cn/CRAN/](https://mirrors.tuna.tsinghua.edu.cn/CRAN/)
Please install **R version 4.3.1** for macOS.

3. Project Structure

The main directories and files in this project are organized as follows:
├── DSP.Rproj/
├── Rdata/
├── GSEA/
├── R_Script/
├── TCellSI-main/
└── README.md

 3.1 Rdata
**Description**:
  The `Rdata` directory contains **essential R data files** required for running the analyses.
**Note**:
  These files are mandatory. Do not delete or rename them.

 3.2 GSEA
**Description**:
  The `GSEA` directory contains R data files required for the **04_Volcano_GSEA analysis**.
**Purpose**:
  Used for Gene Set Enrichment Analysis (GSEA) and volcano plot–related analyses.

 3.3 R_Script
**Description**:
  The `R_Script` directory contains all **R scripts used for data analysis and visualization**.
**Includes**:
  * Data preprocessing
  * Statistical analysis

 3.4 TCellSI-main
**Description**:
  The `TCellSI-main` directory contains the **implementation and resources of the TCellSI algorithm**.
**Purpose**:
  Used to run the TCellSI algorithm and its associated analysis workflows.

4. Usage Instructions (Overview)

1. Verify that your system and R version meet the **System Requirements**
2. Install all required R packages (if an installation script is provided, run it first)
3. Ensure that all required data files are present in the `Rdata` and `GSEA` directories
4. Run the analysis scripts in the `R_Script` directory in the appropriate order
5. For TCellSI-related analysis, navigate to the `TCellSI-main` directory and follow its internal instructions

5. Notes

* Do not modify the directory structure
* If errors occur, please check:
  * R version compatibility
  * Integrity of data files
  * Successful installation of required R packages