# plaeApp_manuscript

Codebase for the PLAE web app manuscript

# Usage
Requires the SQLite database for plae (BIG)

  - 1
  ```
  wget http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/MOARTABLES__anthology_limmaFALSE___4000-counts-universe-study_accession-scANVIprojection-15-5-0.1-50-20__pointRelease01.sqlite.gz
  ```
  - 2
  ```wget http://hpc.nih.gov/~mcgaugheyd/scEiaD/2022_03_22/meta_filter.fst
  ```
  - 3
  ```
  Rscript ~/this/repo/src/ensdb.R
  Rscript ~/this/repo/src/diff_data_pull.R # first edit paths
  Rscript ~/this/repo/src/pubmed_query.R
  ```
  - 4
  Knit manuscript_00.Rmd
