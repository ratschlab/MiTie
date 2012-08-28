function eval_dir(res_dir, mode)

genes = get_eval_genes(mode);

fnsave = sprintf('%s/result_%s.mat', res_dir, mode); 
pgenes = collect_results(res_dir);
pgenes = closed_to_half_open(pgenes);
[SN, SP, F, matched, correct, rank, ab_rank] = eval_genes(genes, pgenes);
fprintf('eval res_realign_GPms_ms_mult%i: SN:%.2f SP:%.2f F:%.2f\n', j, SN*100, SP*100, F*100)
save(fnsave, 'SN', 'SP', 'F', 'matched', 'correct', 'rank', 'ab_rank');

