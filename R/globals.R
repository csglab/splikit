if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".I", ".N", ".GRP", ".",

    "index", "start_cor_id", "chr", "start", "strand",
    "end_cor_id", "row_names_mtx", "sample_id",
    "start_cor_group_id", "start_cor_group_count",
    "end_cor_group_id", "end_cor_group_count",
    "row_names_mtx_new", "raw_row_names_mtx",
    "group_type", "cur_group_id", "n",
    "end", "i", "j", "group_id", "x_tot", "x_1", "x_2",


    "ID", "type", "gene_id", "attribute", "gene_name", "seqid",

    # plot_exclusive_junctions
    "transcript_id", "transcript_name", "exon_number",
    "j_start", "j_end", "exclusive", "targeted",
    "n_tx", "has_excl", "x", "y", "y_label", "xend", "attr", "seqname",
    "i.start", "i.end", "row_names_mtx", "is_annot",
    "observed_in_eventdata"
  ))
}
