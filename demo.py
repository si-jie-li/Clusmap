import clusmap as cm
# bulk
rna = cm.import_data("bulk-data.tsv") # can be other formats
rna = cm.preprocess(rna)
state = cm.gen_mod(rna, outdir=".")
hm = cm.bulk_hm(rna, state, outdir=".",norm_method="z_score",hm_args={'cmap':'coolwarm'})
# sc
pb = cm.compute_pseudo_bulk("sc-data.h5ad",celltype_key = "celltype")
cm.pseudo_bulk_hm(rna, hm, pseudo_bulk_df=pb, log_base=2, 
                  outdir=".", norm_method="z_score", pseudo_hm_args={'cmap':'coolwarm'})

cm.mod_GO("./HM_ModGene.csv", organism="Human", 
          GO_category=["BP", "MF", "CC"], outdir=".")  # default all modules, can specify through mod=[1,3,5]

cm.cluster_sample_stats(rna, "./HM_ModGene.csv", hm=hm, stats=["mean", "std", "cv"],outdir = ".")

# split, merge and assign, e.g.
state.merge(1, 2)
cm.bulk_hm(rna, state,outdir = None)