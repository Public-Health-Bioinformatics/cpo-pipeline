# -*- coding: utf-8 -*-
"""cpo-pipeline.typing.parsers.result_parsers

This module provides functions for parsing result files generated by tools
during the Typing phase of the cpo-pipeline.
"""

import csv
import pandas
import numpy
import pprint

def parse_mlst_result(path_to_mlst_result):
    """
    Args:
        path_to_mlst_result (str): Path to the mlst report file.

    Returns:
        list of dict: Parsed mlst report.
        For example:
        [
            {
                'contig_file': 'SAMPLE-ID.fa',
                'scheme_id': 'ecoli',
                'sequence_type': '405',
                'multi_locus_alleles': {
	            'adk': '35',
	            'fumc': '37',
	            'gyrB': '29',
	            'icd': '25',
	            'mdh': '4',
	            'purA': '5',
	            'recA': '73'
                }
            }
        ]
    """
    mlst_results = []
    #read in the mlst result
    with open(path_to_mlst_result) as mlst_result_file:
        reader = csv.reader(mlst_result_file, delimiter='\t')
        for row in reader:
            mlst_result = {}
            mlst_result['contig_file'] = row[0]
            mlst_result['scheme_id'] = row[1]
            mlst_result['sequence_type'] = row[2]
            mlst_result['multi_locus_alleles'] = {}
            for field in row[3:]:
                (locus, allele) = tuple(field.replace(')', '').split('('))
                mlst_result['multi_locus_alleles'][locus] = allele
            mlst_results.append(mlst_result)
    return mlst_results

def parse_mob_recon_contig_report(path_to_mob_recon_contig_report):
    """
    Args:
        path_to_mob_recon_contig_report (str): Path to the mob_recon contig_report file.

    Returns:
        list of dict: Parsed mob_recon contig_report.
        For example:
        [
            {
                'file_id': 'contigs.fa',
                'cluster_id': '683',
                'contig_id': 'contigs.fa|contig00054_len=2672_cov=424.9_corr=0_origname=NODE_54_length_2672_cov_424.949312_pilon_sw=shovill-spades/1.0.1_date=20181024',
                'contig_length': 2672,
                'circularity_status': 'Incomplete',
                'rep_type': 'IncL/M',
                'rep_type_accession': '000148__NC_021488_00028',
                'relaxase_type': 'MOBP',
                'relaxase_type_accession': 'NC_004464_00056',
                'mash_nearest_neighbor': 'JX988621',
                'mash_neighbor_distance': 0.00560872,
                'repetitive_dna_id': 'KX646543',
                'match_type': 'ISL3',
                'score': 10178.0,
                'contig_match_start': 22,
                'contig_match_end': 5535
            },
            ...
        ]
        Note: several of these fields may be empty for some records.
    """
    def parse_value_maybe(value):
        if value == "":
            return None
        else:
            return value

    # Ideally we would let the csv.DictReader pick these up automatically
    # from the header, but there's a bug in mob_recon that adds extra whitespace
    # in front of 'mash_nearest_neighbor' in the header, so we list them explicitly.
    mob_recon_contig_report_fieldnames = [
        'file_id',
        'cluster_id',
        'contig_id',
        'contig_length',
        'circularity_status',
        'rep_type',
        'rep_type_accession',
        'relaxase_type',
        'relaxase_type_accession',
        'mash_nearest_neighbor',
        'mash_neighbor_distance',
        'repetitive_dna_id',
        'match_type',
        'score',
        'contig_match_start',
        'contig_match_end'
    ]
    mob_recon_contig_report_results = []
    with open(path_to_mob_recon_contig_report) as mob_recon_contig_report_file:
        reader = csv.DictReader(mob_recon_contig_report_file, fieldnames=mob_recon_contig_report_fieldnames, delimiter='\t')
        next(reader) # skip header
        integer_fields = ['contig_length', 'contig_match_start', 'contig_match_end']
        float_fields = ['mash_neighbor_distance', 'score']
        for row in reader:
            for key in row.keys():
                row[key] = parse_value_maybe(row[key])
            for key in integer_fields:
                if row[key]:
                    row[key] = int(row[key])
            for key in float_fields:
                if row[key]:
                    row[key] = float(row[key])
            mob_recon_contig_report_results.append(row)

    return mob_recon_contig_report_results

def parse_mob_recon_mobtyper_aggregate_report(path_to_mob_recon_mobtyper_aggregate_report):
    """
    Args:
        path_to_mob_recon_mobtyper_aggregate_report (str): Path to the mob_recon mobtyper_aggregate_report file.

    Returns:
        list of dict: Parsed mob_recon mobtyper_aggregate_report.
        For example:
        [
            {
                'file_id': 'plasmid_683.fa',
                'num_contigs': 35,
                'total_length': 532060,
                'percent_gc': 52.48656166597752,
                'rep_types': [
                    'IncL/M',
                    'rep_cluster_1254'
                ],
                'rep_type_accessions': [
                    '000148__NC_021488_00028',
                    '000562__KT346360_00001'
                ],
                'relaxase_types': [
                    'MOBP'
                ],
                'relaxase_type_accessions': [
                    'NC_004464_00056'
                ],
                'mpf_type': 'MPF_1',
                'mpf_type_accessions': [
                    'NC_004464_00073',
                    'NC_005246_00069',
                    'NC_019154_00069',
                    'NC_004464_00068',
                    'NC_004464_00066',
                    'NC_019344_00078',
                    'NC_019063_00094',
                    'NC_005246_00057',
                    'NC_004464_00105'
                ],
                'orit_types': [
                    '-'
                ],
                'orit_accessions': [
                    '-'
                ],
                'predicted_mobility': 'Conjugative',
                'mash_nearest_neighbor': 'JX988621',
                'mash_neighbor_distance': 0.00560872,
                'mash_neighbor_cluster': '683'
            },
            ...
        ]
    """
    # A few of the field names in the mobtyper_aggregate_report
    # include special characters '(' and ')' or aren't very descriptive
    # 'gc' -> 'percent_gc' so we list them explicitly here instead of reading
    # from the report header.
    mob_recon_mobtyper_aggregate_report_fieldnames = [
        'file_id',
        'num_contigs',
        'total_length',
        'percent_gc',
        'rep_types',
        'rep_type_accessions',
        'relaxase_types',
        'relaxase_type_accessions',
        'mpf_type',
        'mpf_type_accessions',
        'orit_types',
        'orit_accessions',
        'predicted_mobility',
        'mash_nearest_neighbor',
        'mash_neighbor_distance',
        'mash_neighbor_cluster'
    ]
    mob_recon_mobtyper_aggregate_report_results = []
    with open(path_to_mob_recon_mobtyper_aggregate_report) as mob_recon_mobtyper_aggregate_report_file:
        reader = csv.DictReader(mob_recon_mobtyper_aggregate_report_file,
                                fieldnames=mob_recon_mobtyper_aggregate_report_fieldnames,
                                delimiter='\t')
        next(reader) # skip header
        integer_fields = ['num_contigs', 'total_length']
        float_fields = ['percent_gc', 'mash_neighbor_distance']
        array_fields = [
            'rep_types',
            'rep_type_accessions',
            'relaxase_types',
            'relaxase_type_accessions',
            'mpf_type_accessions',
            'orit_types',
            'orit_accessions'
        ]
        for row in reader:
            for key in integer_fields:
                row[key] = int(row[key])
            for key in float_fields:
                row[key] = float(row[key])
            for key in array_fields:
                row[key] = row[key].split(',')
            mob_recon_mobtyper_aggregate_report_results.append(row)

    return mob_recon_mobtyper_aggregate_report_results


def parse_mobsuite_result(path_to_mobsuite_result):
    """
    Args:
        path_to_mobsuite_result (str): Path to the kraken report file.

    Returns:
        dict: Parsed mobsuite report.
        For example:
        { 'ecoli': { 'contig_file': '/path/to/contig.fa',
                     'species_id': 'ecoli',
                     'sequence_type': '405',
                     'species':, 'Escherichia  ;Shigella ',
                     'row': '/path/to/contig.fa\tecoli\t405\tadk(35)\tfumC(37)\tgyrB(29)\ticd(25)\tmdh(4)\tpurA(5)\trecA(73)'
                   }
        }
    """
    mobsuite = {}
    mResult = pandas.read_csv(path_to_mobsuite_result, delimiter='\t', header=0)
    mResult = mResult.replace(numpy.nan, '', regex=True)

    for i in range(len(mResult.index)):
        mr = {}
        mr['file_id'] = str(mResult.iloc[i,0])
        mr['cluster_id'] = str(mResult.iloc[i,1])
        if (mr['cluster_id'] == "chromosome"):
            break
        mr['contig_id'] = str(mResult.iloc[i,2])
        mr['contig_num'] = mr['contig_id'][(mr['contig_id'].find("contig")+6):mr['contig_id'].find("_len=")]
        mr['contig_length'] = int(mResult.iloc[i,3])
        mr['circularity_status'] = str(mResult.iloc[i,4])
        mr['rep_type'] = str(mResult.iloc[i,5])
        mr['rep_type_accession'] = str(mResult.iloc[i,6])
        mr['relaxase_type'] = str(mResult.iloc[i,7])
        mr['relaxase_type_accession'] = str(mResult.iloc[i,8])
        mr['mash_nearest_neighbor'] = str(mResult.iloc[i,9])
        mr['mash_neighbor_distance'] = float(mResult.iloc[i,10])
        mr['repetitive_dna_id'] = str(mResult.iloc[i,11])
        mr['match_type'] = str(mResult.iloc[i,12])
        if (mr['match_type'] == ""):
            mr['score'] = -1
            mr['contig_match_start'] = -1
            mr['contig_match_end'] = -1
        else:
            mr['score'] = int(mResult.iloc[i,13])
            mr['contig_match_start'] = int(mResult.iloc[i,14])
            mr['contig_match_end'] = int(mResult.iloc[i,15])
        mr['row'] = "\t".join(str(x) for x in mResult.ix[i].tolist())
        mobsuite[mr['contig_id']]=(mr)
    return mobsuite

def parse_mobsuite_plasmids(path_to_mobsuite_result):
    """
    Args:
        path_to_mobsuite_result (str): Path to the mobsuite result file.

    Returns:
        dict: Parsed mobsuite result.
        For example:

        {'plasmid_133.fasta': {'predicted_mobility': 'Mobilizable',
                               'file_id': 'plasmid_133.fasta',
                               'gc': 56,
                               'mash_nearest_neighbor': 'CP016867',
                               'mash_neighbor_cluster': 133,
                               'mash_neighbor_distance': 0.000264084,
                               'mpf_type': '-',
                               'mpf_type_accession': '-',
                               'num_contigs': 16,
                               'orit_accession': 'CP018941',
                               'orit_type': 'MOBP',
                               'relaxase_type': 'MOBP',
                               'relaxase_type_accession': 'NC_013090_00004',
                               'rep_typeAccession': '000019__CP000651_00005',
                               'rep_types': 'ColRNAI_rep_cluster_1987',
                               'row': [ 'plasmid_133.fasta\t16\t53600\t56.1791044776\t',
                                        'ColRNAI_rep_cluster_1987\t',
                                        '000019__CP000651_00005\tMOBP\tNC_013090_00004\t',
                                        '-\t-\tMOBP\tCP018941\tMobilizable\tCP016867\t',
                                        '0.000264084\t133'
                                      ],
                               'total_length': 53600
                               },
         'plasmid_182.fasta': {...
                              },
         ...
        }
    """
    mobsuite = {}
    mResults = pandas.read_csv(path_to_mobsuite_result, delimiter='\t', header=0)
    mResults = mResults.replace(numpy.nan, '', regex=True)

    for i in range(len(mResults.index)):
        mr = {}
        mr['file_id'] = str(mResults.iloc[i,0])
        mr['num_contigs'] = int(mResults.iloc[i,1])
        mr['total_length'] = int(mResults.iloc[i,2])
        mr['gc'] = int(mResults.iloc[i,3])
        mr['rep_types'] = str(mResults.iloc[i,4])
        mr['rep_typeAccession'] = str(mResults.iloc[i,5])
        mr['relaxase_type'] = str(mResults.iloc[i,6])
        mr['relaxase_type_accession'] = str(mResults.iloc[i,7])
        mr['mpf_type'] = str(mResults.iloc[i,8])
        mr['mpf_type_accession'] = str(mResults.iloc[i,9])
        mr['orit_type'] = str(mResults.iloc[i,10])
        mr['orit_accession'] = str(mResults.iloc[i,11])
        mr['predicted_mobility'] = str(mResults.iloc[i,12])
        mr['mash_nearest_neighbor'] = str(mResults.iloc[i,13])
        mr['mash_neighbor_distance'] = float(mResults.iloc[i,14])
        mr['mash_neighbor_cluster'] = int(mResults.iloc[i,15])
        mr['row'] = "\t".join(str(x) for x in mResults.ix[i].tolist())
        mobsuite[mr['file_id']] = mr
    return mobsuite

def parse_abricate_result(path_to_abricate_result):
    """
    Args:
        path_to_resfinder_result (str): Path to the abricate report file.
    
    Returns:
        list of dict: Parsed abricate report.
        For example:
        [
            {
                'file': 'contigs.fa',
                'sequence': 'contig00044',
                'start': 3183,
                'end': 3995,
                'gene': 'NDM-1',
                'coverage': '1-813/813',
                'coverage_map': '===============',
                'gaps': '0/0',
                'percent_coverage': 100.00,
                'percent_identity': 100.00,
                'database': 'bccdc',
                'accession': 'CAZ39946.1',
                'product': '  subclass B1 metallo-beta-lactamase NDM-1 '
            },
            ...
        ]
    """
    abricate_report_fieldnames = [
        'file',
        'sequence',
        'start',
        'end',
        'gene',
        'coverage',
        'coverage_map',
        'gaps',
        'percent_coverage',
        'percent_identity',
        'database',
        'accession',
        'product'
    ]
    abricate_report_results = []
    with open(path_to_abricate_result) as abricate_report_file:
        reader = csv.DictReader(abricate_report_file, fieldnames=abricate_report_fieldnames, delimiter='\t')
        next(reader) # skip header
        integer_fields = ['start', 'end']
        float_fields = ['percent_coverage', 'percent_identity']
        for row in reader:
            for key in integer_fields:
                row[key] = int(row[key])
            for key in float_fields:
                row[key] = float(row[key])
            abricate_report_results.append(row)

    return abricate_report_results
    
def parse_resfinder_result(path_to_resfinder_results, plasmid_contigs, likely_plasmid_contigs):
    """
    Args:
        path_to_resfinder_result (str): Path to the resfinder report file.
        plasmid_contings (str):
        likely_plasmid_contigs ():

    Returns:
        dict: Parsed resfinder report.
        For example:

        {'NDM-1': {'accession': 'CAZ39946.1',
                   'coverage': '1-813/813',
                   'coverage_map': '===============',
                   'database': 'bccdc',
                   'end': 3995,
                   'file': '/projects/cpo/analyses/2018-10-01_validation/Cfr-ST22/cpo_pipeline_module_1_assembly/2018-10-24_7d537af_output/contigs/BC11-Cfr001.fa',
                   'gaps': '0/0',
                   'gene': 'NDM-1',
                   'percent_coverage': 100.0,
                   'percent_identity': 100.0,
                   'product': '  subclass B1 metallo-beta-lactamase NDM-1 ',
                   'row': [ '/projects/cpo/analyses/2018-10-01_validation/Cfr-ST22/cpo_pipeline_module_1_assembly/2018-10-24_7d537af_output/contigs/BC11-Cfr001.fa\t',
                            'contig00044\t3183\t3995\tNDM-1\t1-813/813\t===============\t',
                            '0/0\t100.0\t100.0\tbccdc\tCAZ39946.1\t  subclass B1 ',
                            'metallo-beta-lactamase NDM-1 ',
                          ],
                   'sequence': 'contig00044',
                   'short_gene': 'NDM-1',
                   'source': 'likely plasmid',
                   'start': 3183
                  }
        }
    """
    rFinder = {}
    resFinder = pandas.read_csv(path_to_resfinder_results, delimiter='\t', header=0)
    resFinder = resFinder.replace(numpy.nan, '', regex=True)

    for i in range(len(resFinder.index)):
        rf = {}
        rf['file'] = str(resFinder.iloc[i,0])
        rf['sequence'] = str(resFinder.iloc[i,1])
        rf['start'] = int(resFinder.iloc[i,2])
        rf['end'] = int(resFinder.iloc[i,3])
        rf['gene'] = str(resFinder.iloc[i,4])
        rf['short_gene'] = rf['gene']
        rf['coverage'] = str(resFinder.iloc[i,5])
        rf['coverage_map'] = str(resFinder.iloc[i,6])
        rf['gaps'] = str(resFinder.iloc[i,7])
        rf['percent_coverage'] = float(resFinder.iloc[i,8])
        rf['percent_identity'] = float(resFinder.iloc[i,9])
        rf['database'] = str(resFinder.iloc[i,10])
        rf['accession'] = str(resFinder.iloc[i,11])
        rf['product'] = str(resFinder.iloc[i,12])
        rf['row'] = "\t".join(str(x) for x in resFinder.ix[i].tolist())
        if (rf['sequence'][6:] in plasmid_contigs):
            rf['source'] = "plasmid"
        elif (rf['sequence'][6:] in likely_plasmid_contigs):
            rf['source'] = "likely plasmid"
        else:
            rf['source'] = "likely chromosome"
        rFinder[rf['gene']]=rf
    return rFinder

def parse_rgi_result(path_to_rgi_results, plasmid_contigs, likely_plasmid_contigs):
    """
    Args:
        path_to_rgi_result (str): Path to the rgi report file.
        plasmid_contings (str):
        likely_plasmid_contigs ():

    Returns:
        dict: Parsed rgi report.
        For example:

    {2886: {'AMR_Gene_Family': 'Penicillin-binding protein mutations conferring resistance to beta-lactam antibiotics',
            'ARO': 3004446,
            'Best_Hit_ARO': 'Haemophilus influenzae PBP3 conferring resistance to beta-lactam antibiotics',
            'Best_Hit_Bitscore': 592.038,
            'Best_Identities': 52.75,
            'CARD_Protein_Sequence': 'MVKFNSSRKSGKSKKTI...',
            'Contig': 'contig00006_112 ',
            'Contig_Num': '00006',
            'Cut_Off': 'Strict',
            'Drug_Class': 'carbapenem; cephalosporin; monobactam; penam; cephamycin',
            'ID': 'gnl|BL_ORD_ID|849|hsp_num:1',
            'Model_ID': 2886,
            'Model_type': 'protein variant model',
            'ORF_ID': 'contig00006_112 # 126897 # 128663 # -1 # ID=6_112;partial=00;start_type=ATG;rbs_motif=AGGA/GGAG/GAGG;rbs_spacer=11-12bp;gc_cont=0.546',
            'Orientation': '-',
            'Other_SNPs': '',
            'Pass_Bitscore': 500,
            'Percentage_Length_of_Reference_Sequence': 0.0,
            'Predicted_DNA': 'ATGAAAGCAGCGGCAAAAACGCACAA...',
            'Predicted_Protein': 'MKAAAKTHKPKRQEE...',
            'Resistance_Mechanism': 'antibiotic target alteration',
            'SNPs_in_Best_Hit_ARO': 'S357N, D350N',
            'Start': 126897,
            'Stop': 128663,
            'row': 'contig00006_112 # 126897 # 128663 # -1 # '
               'ID=6_112;partial=00;start_type=ATG;rbs_motif=AGGA/GGAG/GAGG;rbs_spacer=11-12bp;gc_cont=0.546\t'
               'contig00006_112 \t126897\t128663\t-\tStrict\t500\t592.038\t'
               'Haemophilus influenzae PBP3 conferring resistance to '
               'beta-lactam antibiotics\t52.75\t3004446\tprotein variant '
               'model\tS357N, D350N\t\tcarbapenem; cephalosporin; monobactam; '
               'penam; cephamycin\tantibiotic target alteration\t'
               'Penicillin-binding protein mutations conferring resistance to '
               'beta-lactam antibiotics\t'
               'ATGAAAGCAGCGGCAAAAACGCACAAACCAAAACGCCAGGAAGAACAAGCCAACTTTATCAGTTGGCGTTTTGCGTTACTGTGCGGCTGTATTTTGCTGGCACTGGGTTTTCTGCTGGGTCGCGTTGCCTGGCTGCAAATCATCGCGCCGGACATGCTGGTGCGTCAGGGTGATATGCGCTCTCTACGCGTCCAGGAAGTGTCTACATCGCGCGGAATGATTACCGACCGCTCTGGTCGTCCGCTGGCGGTGAGTGTTCCGGTTAAAGCTATATGGGCCGACCCGAAAGAAGTACATGATGCCGGCGGAGTGAGCGTTGGCGAACGCTGGAGAGCGCTGTCAACCGCGCTGAATCTCCCGCTCGATCAACTGGCTTCCCGCATTAACGCGAATCCGAAAGGGCGCTTTATCTATCTGGCGCGTCAGGTAAACCCTGACATGGCTGATTACATCAAAAAACTGAAGCTGCCAGGTATCCATCTGCGCGAAGAATCCCGCCGTTACTACCCTTCCGGGGAAGTGACTGCTCACCTCATCGGATTTACGAACGTCGACAGCCAGGGCATTGAAGGCGTTGAGAAGAGCTTCGACAAGTGGCTCACCGGACAACCGGGCGAGCGTATTGTACGTAAAGACCGGTATGGCCGCGTCATTGAAGATATCTCCTCTACCGACAGCCAGGCGGCGCATAACCTCGCGTTGAGCATTGATGAGCGTTTACAGGCACTGGTCTACCGTGAACTGAATAACGCCGTGGCGTTCAACAAGGCGGAGTCAGGCAGTGCGGTACTGGTGGATGTGAACACCGGTGAAGTGCTGGCAATGGCCAACAGTCCGTCCTACAACCCGAACAATCTCACCGGTACGCCAAAAGACGCGATGCGTAACCGCACCATTACCGACGTGTTTGAACCGGGTTCTACCGTTAAACCGATGGTGGTGATGACCGCGCTGCAGCGCGGTGTGGTACGCGAAAATACGGTCCTCAACACTATCCCTTACCGAATTAATGGTCACGAAATCAAAGACGTGGCACGTTATAGCGAATTAACCCTCACCGGGGTTTTGCAGAAGTCGAGTAACGTCGGTGTTTCCAAGCTGGCGTTAGCGATGCCGTCCTCAGCGTTAGTAGATACTTACTCACGTTTTGGGCTGGGAAAAGCGACCAATTTGGGGTTGGTCGGAGAACGCAGTGGCTTATATCCTCAAAAACAACGGTGGTCTGACATAGAGAGGGCCACCTTCTCTTTCGGCTACGGGCTAATGGTAACGCCGTTACAGTTAGCGCGAGTCTATGCAACGATCGGCAGTTACGGCGTTTATCGTCCACTGTCGATTACTAAAGTTGACCCTCCAGTTCCGGGCGAGCGTATCTTCCCGGAAGCGACCGTACGTACCGTGGTGCACATGATGGAAAGCGTGGCGCTGCCCGGCGGCGGCGGCGTGAAGGCGGCGATTAAAGGTTATCGTATCGCCATTAAAACCGGTACAGCGAAAAAAGTAGGGCCGGATGGCCGCTACATCAACAAATACATTGCTTATACCGCAGGTGTTGCGCCTGCGAGTCAGCCGCGCTTCGCGCTGGTTGTTGTTATCAACGATCCGCAGGCGGGTAAATACTACGGTGGCGCCGTTTCCGCGCCGGTCTTTGGTGCCATCATGGGCGGCGTACTGCGTACCATGAACATCGAGCCGGATGCGCTGACAACGGGCGATAAAAATGAATTTGTGAATAATCAAGGCGAGGCAACAGGTGGCAGATCGTAA\t'
               'MKAAAKTHKPKRQEEQANFISWRFALLCGCILLALGFLLGRVAWLQIIAPDMLVRQGDMRSLRVQEVSTSRGMITDRSGRPLAVSVPVKAIWADPKEVHDAGGVSVGERWRALSTALNLPLDQLASRINANPKGRFIYLARQVNPDMADYIKKLKLPGIHLREESRRYYPSGEVTAHLIGFTNVDSQGIEGVEKSFDKWLTGQPGERIVRKDRYGRVIEDISSTDSQAAHNLALSIDERLQALVYRELNNAVAFNKAESGSAVLVDVNTGEVLAMANSPSYNPNNLTGTPKDAMRNRTITDVFEPGSTVKPMVVMTALQRGVVRENTVLNTIPYRINGHEIKDVARYSELTLTGVLQKSSNVGVSKLALAMPSSALVDTYSRFGLGKATNLGLVGERSGLYPQKQRWSDIERATFSFGYGLMVTPLQLARVYATIGSYGVYRPLSITKVDPPVPGERIFPEATVRTVVHMMESVALPGGGGVKAAIKGYRIAIKTGTAKKVGPDGRYINKYIAYTAGVAPASQPRFALVVVINDPQAGKYYGGAVSAPVFGAIMGGVLRTMNIEPDALTTGDKNEFVNNQGEATGGRS\t'
               'MVKFNSSRKSGKSKKTIRKLTAPETVKQNKPQKVFEKCFMRGRYMLSTVLILLGLCALVARAAYVQSINADTLSNEADKRSLRKDEVLSVRGSILDRNGQLLSVSVPMSAIVADPKTMLKENSLADKERIAALAEELGMTENDLVKKIEKNSKSGYLYLARQVELSKANYIRRLKIKGIILETEHRRFYPRVEEAAHVVGYTDIDGNGIEGIEKSFNSLLVGKDGSRTVRKDKRGNIVAHISDEKKYDAQDVTLSIDEKLQSMVYREIKKAVSENNAESGTAVLVDVRTGEVLAMATAPSYNPNNRVGVKSELMRNRAITDTFEPGSTVKPFVVLTALQRGVVKRDEIIDTTSFKLSGKEIVDVAPRAQQTLDEILMNSSNRGVSRLALRMPPSALMETYQNAGLSKPTDLGLIGEQVGILNANRKRWADIERATVAYGYGITATPLQIARAYATLGSFGVYRPLSITKVDPPVIGKRVFSEKITKDIVGILEKVAIKNKRAMVEGYRVGVKTGTARKIENGHYVNKYVAFTAGIAPISDPRYALVVLINDPKAGEYYGGAVSAPVFSNIMGYALRANAIPQDAEAAENTTTKSAKRIVYIGEHKNQKVN\t'
               '0.0\tgnl|BL_ORD_ID|849|hsp_num:1\t2886',
            'source': 'likely chromosome'
        }
    2661: {...
          }
    ...
    }
    """
    rgiR = {}
    RGI = pandas.read_csv(path_to_rgi_results, delimiter='\t', header=0)
    RGI = RGI.replace(numpy.nan, '', regex=True)

    for i in range(len(RGI.index)):
        r = {}
        r['orf_id'] = str(RGI.iloc[i,0])
        r['contig'] = str(RGI.iloc[i,1])
        r['contig_num'] = r['contig'][6:r['contig'].find("_")]
        r['start'] = int(RGI.iloc[i,2])
        r['stop'] = int(RGI.iloc[i,3])
        r['orientation'] = str(RGI.iloc[i,4])
        r['cut_off'] = str(RGI.iloc[i,5])
        r['pass_bitscore'] = int(RGI.iloc[i,6])
        r['best_hit_bitscore'] = float(RGI.iloc[i,7])
        r['best_hit_aro'] = str(RGI.iloc[i,8])
        r['best_identities'] = float(RGI.iloc[i,9])
        r['aro'] = int(RGI.iloc[i,10])
        r['model_type'] = str(RGI.iloc[i,11])
        r['snps_in_best_hit_aro'] = str(RGI.iloc[i,12])
        r['other_snps'] = str(RGI.iloc[i,13])
        r['drug_class'] = str(RGI.iloc[i,14])
        r['resistance_mechanism'] = str(RGI.iloc[i,15])
        r['amr_gene_family'] = str(RGI.iloc[i,16])
        r['predicted_dna'] = str(RGI.iloc[i,17])
        r['predicted_protein'] = str(RGI.iloc[i,18])
        r['card_protein_sequence'] = str(RGI.iloc[i,19])
        r['percentage_length_of_reference_sequence'] = float(RGI.iloc[i,20])
        r['id'] = str(RGI.iloc[i,21])
        r['model_id'] = int(RGI.iloc[i,22])
        r['row'] = "\t".join(str(x) for x in RGI.ix[i].tolist())
        if (r['contig_num'] in plasmid_contigs):
            r['source'] = "plasmid"
        elif (r['contig_num'] in likely_plasmid_contigs):
            r['source'] = "likely plasmid"
        else:
            r['source'] = "likely chromosome"
        rgiR[r['model_id']] = r
    return rgiR

def parse_rgi_result_txt(path_to_rgi_result):
    """
    Args:
        path_to_rgi_report (str): Path to the tabular rgi report file (SAMPLE-ID.rgi.txt).

    Returns:
        list of dict: Parsed rgi report.
        For example:
        [
            {
                'orf_id': 'contig00007_44 # 46711 # 49833 # -1 # ID=7_44;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;gc_cont=0.547',
                'contig': 'contig00007_44',
                'start': 46711,
                'stop': 49833,
                'orientation': '-',
                'cut_off': 'Strict',
                'pass_bitscore': 1800,
                'best_hit_bitscore': 1894.78,
                'best_hit_aro': 'mdtB',
                'best_identities': 92.7,
                'aro': '3000793',
                'model_type': 'protein homolog model',
                'snps_in_best_hit_aro': [
                    'S357N',
                    'D350N'
                ],
                'other_snps': None,
                'drug_class': 'aminocoumarin antibiotic',
                'resistance mechanism': 'antibiotic efflux',
                'amr_gene_family': 'resistance-nodulation-cell division (RND) antibiotic efflux pump',
                'predicted_dna': 'ATGCAGGTGTTACCTCCTGACAACACAGGCGGACCATCGC...',
                'predicted_protein': 'MQVLPPDNTGGPSRLFILRPVATTLLMVAILLAGII...',
                'card_protein_sequence': 'MQVLPPSSTGGPSRLFIMRPVATTLLMVAILL...',
                'percentage_length_of_reference_sequence': 100.00,
                'id': 'gnl|BL_ORD_ID|776|hsp_num:0',
                'model_id': '820'
            },
            ...
        ]
    """
    rgi_report_fieldnames = [
        'orf_id',
        'contig',
        'start',
        'stop',
        'orientation',
        'cut_off',
        'pass_bitscore',
        'best_hit_bitscore',
        'best_hit_aro',
        'best_identities',
        'aro',
        'model_type',
        'snps_in_best_hit_aro',
        'other_snps',
        'drug_class',
        'resistance mechanism',
        'amr_gene_family',
        'predicted_dna',
        'predicted_protein',
        'card_protein_sequence',
        'percentage_length_of_reference_sequence',
        'id',
        'model_id'
    ]
    rgi_report_results = []
    
    def parse_value_maybe(value):
        if value == "n/a":
            return None
        else:
            return value
        
    with open(path_to_rgi_result) as rgi_report_file:
        reader = csv.DictReader(rgi_report_file, fieldnames=rgi_report_fieldnames, delimiter='\t')
        next(reader) # skip header
        integer_fields = [
            'start',
            'stop',
        ]
        float_fields = [
            'pass_bitscore',
            'best_hit_bitscore',
            'best_identities',
            'percentage_length_of_reference_sequence'
        ]
        array_fields = [
            'snps_in_best_hit_aro',
            'other_snps'
        ]
        for row in reader:
            for key in integer_fields:
                row[key] = int(row[key])
            for key in float_fields:
                row[key] = float(row[key])
            for key in array_fields:
                # 'n/a' => None
                # 'S80I' => ['S80I']
                # 'S357N, D350N' => ['S357N', 'D350N']
                row[key] = row[key].split(', ') if parse_value_maybe(row[key]) else None
            rgi_report_results.append(row)

    return rgi_report_results
