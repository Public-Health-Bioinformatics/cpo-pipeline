# -*- coding: utf-8 -*-
"""cpo-pipeline.assembly.qc.check

This module provides functions for processing QC results
during the QC & Assembly phase of the cpo-pipeline.
"""

def busco_qc_check(busco_results, qc_thresholds):
    """
    BUSCO PASS CRITERIA:
    1. complete singles > 90% of total genes
    2. complte duplicates < 90% of total genes
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
    
def quast_qc_check(quast_results, qc_thresholds):
    """
    QUAST PASS CRITERIA:
    1. total length vs reference length +-10%
    2. percent gc versus reference percent gc +- 5%
    3. genome fraction percent > 90
    Args:
        quast_results (dict): Quast results
        qc_thresholds (dict): Threshold values for determining QC pass/fail
    Returns:
        boolean: Assembly passes our QUAST quality criteria
    """
    total_length = quast_results['total_length']
    reference_length = quast_results['reference_length']
    assembly_percent_gc = quast_results['percent_GC']
    reference_percent_gc = quast_results['reference_percent_GC']
    genome_fraction_percent = quast_results['genome_fraction_percent']
    percent_gc_cutoff = qc_thresholds['quast_percent_gc_cutoff']
    assembly_length_cutoff = qc_thresholds['quast_assembly_length_cutoff']
    genome_fraction_percent_cutoff = qc_thresholds['quast_genome_fraction_percent_cutoff']
    if total_length <= (reference_length * (1 + assembly_length_cutoff)) and \
       total_length >= (reference_length * (1 - assembly_length_cutoff)) and \
       assembly_percent_gc <= (reference_percent_gc * (1 + percent_gc_cutoff)) and \
       assembly_percent_gc >= (reference_percent_gc * (1 - percent_gc_cutoff)) and \
       genome_fraction_percent >= genome_fraction_percent_cutoff:
        return True
    else:
        return False

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
