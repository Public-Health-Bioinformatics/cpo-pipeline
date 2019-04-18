# -*- coding: utf-8 -*-
"""cpo-pipeline.assembly.quality_control

This module provides functions for processing QC results
during the QC & Assembly phase of the cpo-pipeline.
"""

def busco_qc_check(busco_results, qc_thresholds):
    """
    Args:
        busco_results (dict): Busco results
        qc_thresholds (dict): Threshold values for determining QC pass/fail
    Returns:
        boolean: Assembly passes our BUSCO quality criteria
    """
    complete_single = busco_results['complete_single']
    complete_duplicate = busco_results['complete_duplicate']
    total = busco_results['total']
    busco_complete_single_cutoff = qc_thresholds['busco_complete_single_cutoff']
    busco_complete_duplicate_cutoff = qc_thresholds['busco_complete_duplicate_cutoff']
    if (complete_single / total) >= busco_complete_single_cutoff and \
       (complete_duplicate / total) <= busco_complete_duplicate_cutoff:
        return True
    else:
        return False
    
def quast_qc_check(quast_results, estimated_genome_size):
    """
    QUAST PASS CRITERIA:
    1. total length within +/- 10% of expected genome size
    Args:
        quast_results (dict): Quast results
        qc_thresholds (dict): Threshold values for determining QC pass/fail
    Returns:
        boolean: Assembly passes our QUAST quality criteria
    """
    total_length = quast_results['total_length']
    return bool(total_length <= (estimated_genome_size * 1.1) and \
                total_length >= (estimated_genome_size * 0.9))

def fastqc_qc_check(fastqc_results):
    """
    Args:
        fastqc_results (dict): FastQC results
    Returns:
        boolean: Sequence data passes our FastQC quality criteria
    """
    if fastqc_results["basic_statistics"] == "PASS" and \
       fastqc_results["per_base_sequence_quality"] == "PASS" and \
       fastqc_results["sequence_length_distribution"] == "PASS":
        return True
    else:
        return False
