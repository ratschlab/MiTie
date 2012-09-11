

bams{1} = '/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias1_merged_err_1.new.sorted.paired.bam';
bams{2} = '/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias1_merged_err_1.new.sorted.paired.bam,/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/reads_with_errors/bias1_merged_err_1.new.sorted.paired.bam';

genes = load_regions_bin('/fml/ag-raetsch/nobackup/projects/mip/human_sim/mip_mmr_new_align_gtf_regions50000_enum5_1e4_sample1/graph.bin', 10, bams{:});
