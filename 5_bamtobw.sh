#!/bin/bash
#SBATCH --account=PCON0022
#SBATCH --time=20:30:00
#SBATCH --nodes=1 
#SBATCH --mem=32GB

ml python/3.7-2019.10
source activate MACS 

cd /fs/ess/PCON0022/zhenyuwu/Dawei_dualCHIPseq/BMDM_BRD4_GSE151015/align/rep1
bamCoverage --bam SRR11828501.sra.fastq.srt.bam -o RPKM_BMDM_BRD4_rep1.bw --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize 2652783500

cd /fs/ess/PCON0022/zhenyuwu/Dawei_dualCHIPseq/BMDM_BRD4_GSE151015/align/rep2
bamCoverage --bam SRR11828502.sra.fastq.srt.bam -o RPKM_BMDM_BRD4_rep2.bw --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize 2652783500

cd /fs/ess/PCON0022/zhenyuwu/Dawei_dualCHIPseq/BMDM_Spt16_GSE117333/align/rep1
bamCoverage --bam SRR7537468.sra.fastq.srt.bam -o RPKM_BMDM_Spt16_rep1.bw --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize 2652783500

cd /fs/ess/PCON0022/zhenyuwu/Dawei_dualCHIPseq/BMDM_Spt16_GSE117333/align/ctl1
bamCoverage --bam SRR7537462.sra.fastq.srt.bam -o RPKM_BMDM_Spt16_control.bw --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize 2652783500


cd /fs/ess/PCON0022/zhenyuwu/Dawei_dualCHIPseq/C2C12_BRD4_GSE99101/align/rep1
bamCoverage --bam SRR5578772.sra.fastq.srt.bam -o RPKM_C2C12_BRD4_rep1.bw --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize 2652783500

cd /fs/ess/PCON0022/zhenyuwu/Dawei_dualCHIPseq/C2C12_BRD4_GSE99101/align/ctl1
bamCoverage --bam SRR5578786.sra.fastq.srt.bam -o RPKM_C2C12_BRD4_ctrl.bw --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize 2652783500


cd /fs/ess/PCON0022/zhenyuwu/Dawei_dualCHIPseq/C2C12_Spt16_GSE116169/align/rep1
bamCoverage --bam SRR9326257.sra_1.fastq.srt.bam -o RPKM_C2C12_Spt16_rep1.bw --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize 2652783500


cd /fs/ess/PCON0022/zhenyuwu/Dawei_dualCHIPseq/C2C12_Spt16_GSE116169/align/ctl1
bamCoverage --bam SRR9326249.sra_1.fastq.srt.bam -o RPKM_C2C12_Spt16_ctrl.bw --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize 2652783500




