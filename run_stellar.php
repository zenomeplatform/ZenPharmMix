<?php
//main.nf на php с доработками
$ref_file = $argv["1"];
$bam_path = $argv["2"];
$bam_name = $argv["3"];
$report_path = $argv["4"];
$work_dir = "/analyze/soft/tools/StellarPGx/MIX";
$params_build = 'hg19';
$d_base = "{$work_dir}/database";
$res_base = "{$work_dir}/resources";
$caller_base = "{$work_dir}/scripts";

$debug38 = '';
$debug37 = '';
if ($params_build == 'hg38') $debug38 = '--minimum_extract_score_over_homref=0 ';
else $debug37 = '--minimum_extract_score_over_homref=0 ';

$genes_need = array ("1" => "cyp2d6", "cyp2c19", "cyp2c9", "cyp2b6", "cyp1a2", "cyp3a4", "cyp3a5", "cypor");

exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_reports:v1 rm -rf {$bam_path}/stellar_tmp");
exec("mkdir {$bam_path}/stellar_tmp");
exec("mkdir -p {$bam_path}/{$report_path}");

$main_report = '';

for ($k = 1; isset($genes_need["$k"]); $k++) {
	$params_gene = $genes_need["$k"];
	
echo ("!!! Running {$params_gene}\n");
	//Константы от версии генома
	$coord_file = file("{$work_dir}/fragments_database.txt");
	for ($t = 1; isset($coord_file["$t"]); $t++) {
		$ll = explode("\t", trim($coord_file["$t"]));
		if (($ll["0"] == $params_build)&&($ll["1"] == $params_gene)) break;
	}
	$chrom = $ll["2"];
	$region_a1 = "{$chrom}:{$ll["3"]}";
	$region_a2 = $ll["3"];
	$region_b1 = "{$chrom}:{$ll["4"]}";
	$region_b2 = $ll["4"];
	$db = "{$d_base}/{$params_gene}/{$params_build}";
	if ($params_build == "hg19") $db = "{$d_base}/{$params_gene}/b37";
	$res_dir = "{$res_base}/{$params_gene}/cyp_{$params_build}";
	//call_snvs1
echo ("!!! Call_SNVS1\n");
	exec("mkdir -p {$bam_path}/stellar_tmp/{$bam_name}_var_1/{$params_gene}");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_reports:v1 graphtyper genotype {$ref_file} --sam={$bam_path}/{$bam_name}.MIX.bam --region={$region_a1} --output={$bam_path}/stellar_tmp/{$bam_name}_var_1/{$params_gene} --prior_vcf={$res_dir}/common_plus_core_var.vcf.gz -a {$debug38} --bamshrink_max_fraglen=250 --no_filter_on_coverage --no_filter_on_strand_bias");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 bcftools concat {$bam_path}/stellar_tmp/{$bam_name}_var_1/{$params_gene}/{$chrom}/*.vcf.gz -o {$bam_path}/stellar_tmp/{$bam_name}_var_1/{$params_gene}/{$region_a2}.vcf"); 
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 bgzip -f {$bam_path}/stellar_tmp/{$bam_name}_var_1/{$params_gene}/{$region_a2}.vcf");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 tabix -f {$bam_path}/stellar_tmp/{$bam_name}_var_1/{$params_gene}/{$region_a2}.vcf.gz");
	//call_snvs2
echo ("!!! Call_SNVS2\n");
	exec("mkdir -p {$bam_path}/stellar_tmp/{$bam_name}_var_2/{$params_gene}");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_reports:v1 graphtyper genotype {$ref_file} --sam={$bam_path}/{$bam_name}.MIX.bam --region={$region_a1} --output={$bam_path}/stellar_tmp/{$bam_name}_var_2/{$params_gene} --prior_vcf=${res_dir}/common_plus_core_var.vcf.gz -a {$debug38}{$debug37} --bamshrink_max_fraglen=1000 --bamshrink_min_matching=40 --genotype_aln_min_support=2 --genotype_dis_min_support=2 --no_filter_on_proper_pairs --no_filter_on_read_bias --no_filter_on_strand_bias --minimum_extract_variant_support=1");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 bcftools concat {$bam_path}/stellar_tmp/{$bam_name}_var_2/{$params_gene}/{$chrom}/*.vcf.gz -o {$bam_path}/stellar_tmp/{$bam_name}_var_2/{$params_gene}/{$region_a2}.vcf");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 bgzip -f {$bam_path}/stellar_tmp/{$bam_name}_var_2/{$params_gene}/{$region_a2}.vcf");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 tabix -f {$bam_path}/stellar_tmp/{$bam_name}_var_2/{$params_gene}/{$region_a2}.vcf.gz");
	//call_sv_del
echo ("!!! Call_SV_Del\n");
	exec("mkdir -p {$bam_path}/stellar_tmp/{$bam_name}_sv_del/{$params_gene}");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_reports:v1 graphtyper genotype_sv {$ref_file} --sam={$bam_path}/{$bam_name}.MIX.bam --region={$region_a1} --output={$bam_path}/stellar_tmp/{$bam_name}_sv_del/{$params_gene} {$res_dir}/sv_test.vcf.gz");
	//call_sv_dup
echo ("!!! Call_SV_Dup\n");
	exec("mkdir -p {$bam_path}/stellar_tmp/{$bam_name}_sv_dup/{$params_gene}");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_reports:v1 graphtyper genotype_sv {$ref_file} --sam={$bam_path}/{$bam_name}.MIX.bam --region={$region_a1} --output={$bam_path}/stellar_tmp/{$bam_name}_sv_dup/{$params_gene} {$res_dir}/sv_test3.vcf.gz");
	//get_depth
echo ("!!! Depth\n");
	$test3_file = "test3.bed";
	exec("samtools bedcov --reference {$ref_file} {$res_dir}/{$test3_file} {$bam_path}/{$bam_name}.WGS.bam > {$bam_path}/stellar_tmp/{$bam_name}_{$params_gene}_ctrl.depth");
	//format_snvs
echo ("!!! Format_SNVS\n");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 bcftools isec -p {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene} -Oz {$bam_path}/stellar_tmp/{$bam_name}_var_1/{$params_gene}/{$region_a2}.vcf.gz {$bam_path}/stellar_tmp/{$bam_name}_var_2/{$params_gene}/{$region_a2}.vcf.gz");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 bcftools concat -a -D -r {$region_b1} {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/0000.vcf.gz {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/0001.vcf.gz {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/0002.vcf.gz -Oz -o {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/{$bam_name}_{$region_b2}.vcf.gz");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 tabix {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/{$bam_name}_{$region_b2}.vcf.gz");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 bash {$work_dir}/bcftools_pipe.sh {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/{$bam_name}_{$region_b2}.vcf.gz {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/{$bam_name}_all_norm.vcf.gz");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 tabix {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/{$bam_name}_all_norm.vcf.gz");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 cp {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/{$bam_name}_all_norm.vcf.gz {$bam_path}/stellar_tmp/{$bam_name}_{$params_gene}.vcf.gz");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 tabix {$bam_path}/stellar_tmp/{$bam_name}_{$params_gene}.vcf.gz");
	//get_core_var
echo ("!!! Get_Core_Var\n");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 bcftools isec {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/{$bam_name}_all_norm.vcf.gz {$res_dir}/allele_def_var.vcf.gz -p {$bam_path}/stellar_tmp/{$bam_name}_int/{$params_gene} -Oz");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 bash {$work_dir}/bcftools_pipe.sh {$bam_path}/stellar_tmp/{$bam_name}_int/{$params_gene}/0002.vcf.gz {$bam_path}/stellar_tmp/{$bam_name}_int/{$params_gene}/{$bam_name}_core.vcf.gz");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 tabix {$bam_path}/stellar_tmp/{$bam_name}_int/{$params_gene}/{$bam_name}_core.vcf.gz");
echo ("!!! Analyse_1\n");
	//analyse_1
	$bquery = "-f'%ID\t%ALT\t[%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n'";
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 bcftools query {$bquery} {$bam_path}/stellar_tmp/{$bam_name}_sv_del/{$params_gene}/{$chrom}/{$region_a2}.vcf.gz > {$bam_path}/stellar_tmp/{$bam_name}_sv_del/{$params_gene}/{$bam_name}_gene_del_summary.txt");
echo ("!!! Analyse_2\n");
	//analyse_2
	$bquery = "-f'%POS~%REF>%ALT\t[%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n' -i'GT=\"alt\"'";
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 bcftools query {$bquery} {$bam_path}/stellar_tmp/{$bam_name}_sv_dup/{$params_gene}/{$chrom}/{$region_a2}.vcf.gz > {$bam_path}/stellar_tmp/{$bam_name}_sv_dup/{$params_gene}/{$bam_name}_gene_dup_summary.txt");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 bcftools query {$bquery} {$bam_path}/stellar_tmp/{$bam_name}_int/{$params_gene}/{$bam_name}_core.vcf.gz >> {$bam_path}/stellar_tmp/{$bam_name}_sv_dup/{$params_gene}/{$bam_name}_gene_dup_summary.txt");
echo ("!!! Analyse_3\n");
	//analyse_3
	$bquery1 = "-f'[%POS~%REF>%ALT~%GT\n]'";
	$bquery2 = "-f '%POS~%REF>%ALT\n'";
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 bcftools query {$bquery1} {$bam_path}/stellar_tmp/{$bam_name}_int/{$params_gene}/{$bam_name}_core.vcf.gz -o {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/{$bam_name}_core_snvs.dip");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 bcftools query {$bquery2} {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/{$bam_name}_all_norm.vcf.gz -o {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/{$bam_name}_full.dip");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 bcftools query {$bquery1} {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/{$bam_name}_all_norm.vcf.gz -o {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/{$bam_name}_gt.dip");
	exec("docker run -v '/analyze/':'/analyze/' quantum/zenome_merging:v1 chmod -R 777 {$bam_path}/stellar_tmp/");
	//call_stars
echo ("!!! Call stars\n");
	$caller_dir = "{$caller_base}/{$params_gene}/{$params_build}/bin";
	if ($params_build == "hg19") $caller_dir = "{$caller_base}/{$params_gene}/b37/bin";
	exec("python3 {$caller_dir}/stellarpgx.py {$db}/diplo_db_debugged2.dbs {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/{$bam_name}_core_snvs.dip {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/{$bam_name}_full.dip {$bam_path}/stellar_tmp/{$bam_name}_var/{$params_gene}/{$bam_name}_gt.dip {$db}/genotypes4.dbs {$bam_path}/stellar_tmp/{$bam_name}_sv_del/{$params_gene}/{$bam_name}_gene_del_summary.txt {$bam_path}/stellar_tmp/{$bam_name}_sv_dup/{$params_gene}/{$bam_name}_gene_dup_summary.txt {$bam_path}/stellar_tmp/{$bam_name}_{$params_gene}_ctrl.depth {$db}/haps_var_new.dbs {$db}/a_scores.dbs > {$bam_path}/stellar_tmp/{$bam_name}_{$params_gene}.alleles");  
	exec("cp {$bam_path}/stellar_tmp/{$bam_name}_{$params_gene}.alleles {$bam_path}/{$report_path}/{$bam_name}_{$params_gene}.alleles");
	$main_report .= file_get_contents("{$bam_path}/{$report_path}/{$bam_name}_{$params_gene}.alleles");
	$main_report .= "\n";
}

$test_id_cut = substr($bam_name,0,8);
$fww = fopen("{$bam_path}/{$report_path}/{$test_id_cut}.stellar.txt", "w");
fwrite($fww, $main_report);
fclose($fww);
?>
