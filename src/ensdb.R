# fix  occasional NA
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
ensembl_id <- keys(EnsDb.Hsapiens.v86, keytype = "GENEID")
ens_sym_table <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembl_id, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
write_tsv(ens_sym_table, file = 'data/hs.ensdb.v86.symbol.geneid.tsv.gz')




library(org.Mm.eg.db)
library(EnsDb.Mmusculus.v79)
ensembl_id <- keys(EnsDb.Mmusculus.v79, keytype = "GENEID")
ens_sym_table <- ensembldb::select(EnsDb.Mmusculus.v79, keys= ensembl_id, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
write_tsv(ens_sym_table, file = 'data/mm.ensdb.v79.symbol.geneid.tsv.gz')
