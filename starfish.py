#!/usr/bin/env python3

import argparse
import itertools
import shutil
import os
import string
import subprocess as sp
import pysam as ps
from os.path import join, basename, exists

can_draw_venn = True

try:
    import pyvenn.venn as pv
    import matplotlib.pyplot as plt
except ImportError:
    can_draw_venn = False

def call(cmd, quite=True):
    if quite:
        with open(os.devnull, "w") as f:
            sp.call(cmd, stdout=f, stderr=f)
    else:
        sp.call(cmd)

def index_vcf(vcf):
    call(['tabix', vcf])

def remove_vcf(vcf):
    os.remove(vcf)
    if exists(vcf + '.tbi'):
        os.remove(vcf + '.tbi')

def count_records(vcf_name):
    vcf = ps.VariantFile(vcf_name)
    return sum(1 for rec in vcf)

def run_rtg(rtg, ref_sdf, lhs_vcf, rhs_vcf, out_dir,
            bed_regions=None, all_records=False, squash_ploidy=False, sample=None, threads=None):
    cmd = [rtg, 'vcfeval', '-t', ref_sdf, '-b', lhs_vcf, '-c', rhs_vcf, '-o', out_dir]
    if bed_regions is not None:
        cmd += ['--bed-regions', bed_regions]
    if all_records:
        cmd.append('--all-records')
    if sample is not None:
        cmd += ['--sample', sample]
    if threads is not None:
        cmd += ['--threads', str(threads)]
    call(cmd)

def intersect(lhs_label, lhs_vcf, rhs_label, rhs_vcf, args):
    tmp_dir = join(args.output, 'temp')
    run_rtg(args.rtg, args.sdf, lhs_vcf, rhs_vcf, tmp_dir,
            bed_regions=args.regions,
            all_records=args.all_records,
            squash_ploidy=args.ignore_genotypes,
            sample=args.sample,
            threads=args.threads)
    lhs_and_rhs = join(args.output, lhs_label + '_and_' + rhs_label + '.vcf.gz')
    rhs_and_lhs = join(args.output, rhs_label + '_and_' + lhs_label + '.vcf.gz')
    lhs_not_rhs = join(args.output, lhs_label + '_not_' + rhs_label + '.vcf.gz')
    rhs_not_lhs = join(args.output, rhs_label + '_not_' + lhs_label + '.vcf.gz')
    tp_baseline = join(tmp_dir, 'tp-baseline.vcf.gz')
    tp = join(tmp_dir, 'tp.vcf.gz')
    fn = join(tmp_dir, 'fn.vcf.gz')
    fp = join(tmp_dir, 'fp.vcf.gz')
    shutil.move(tp_baseline, lhs_and_rhs)
    shutil.move(tp_baseline + '.tbi', lhs_and_rhs + '.tbi')
    shutil.move(tp, rhs_and_lhs)
    shutil.move(tp + '.tbi', rhs_and_lhs + '.tbi')
    shutil.move(fn, lhs_not_rhs)
    shutil.move(fn + '.tbi', lhs_not_rhs + '.tbi')
    shutil.move(fp, rhs_not_lhs)
    shutil.move(fp + '.tbi', rhs_not_lhs + '.tbi')
    shutil.rmtree(tmp_dir)
    return lhs_and_rhs, rhs_and_lhs, lhs_not_rhs, rhs_not_lhs

def naive_common(vcfs, out):
    assert len(vcfs) > 1
    cmd = ['bcftools', 'isec', '-Oz', '-o', out, '-n=' + str(len(vcfs)), '-w1'] + vcfs
    call(cmd)

def naive_complement(head_vcf, tail_vcfs, out):
    assert len(tail_vcfs) > 0
    cmd = ['bcftools', 'isec', '-C', '-n=1', '-w1', '-Oz', '-o', out, head_vcf] + tail_vcfs
    call(cmd)

def naive_intersect(head_vcf, tail_vcfs, tail_mask, out):
    positive_tails = [tail_vcfs[i] for i, mask in enumerate(tail_mask) if mask]
    assert len(positive_tails) > 0
    negative_tails = [tail_vcfs[i] for i, mask in enumerate(tail_mask) if not mask]
    if len(negative_tails) > 0:
        common_positive_temp_vcf = out + '.temp'
        naive_common([head_vcf] + positive_tails, common_positive_temp_vcf)
        index_vcf(common_positive_temp_vcf)
        naive_complement(common_positive_temp_vcf, negative_tails, out)
        remove_vcf(common_positive_temp_vcf)
    else:
        naive_common([head_vcf] + positive_tails, out)

def plot_ven(labels, names):
    if len(names) == 2:
        return pv.venn2(labels, names=names)
    elif len(names) == 3:
        return pv.venn3(labels, names=names)
    elif len(names) == 4:
        return pv.venn4(labels, names=names)
    elif len(names) == 5:
        return pv.venn5(labels, names=names)
    elif len(names) == 6:
        return pv.venn6(labels, names=names)
    else:
        raise Error("Unsupported Venn")

def main(args):
    if len(args.variants) < 2:
        print("There must be at least two VCFs to intersect")
        return
    
    if args.names is not None:
        if len(args.variants) != len(args.names):
            print("There must be one name per VCF")
            return
    
    if not exists(args.output):
        os.makedirs(args.output)
    
    vcfs = args.variants
    labels = string.ascii_uppercase[:len(vcfs)]
    
    if len(vcfs) > len(labels):
        print("Maximum number of VCFs to intersect is", len(labels))
        return
    
    labelled_vcfs = [(labels[i], vcf) for i, vcf in enumerate(vcfs)]
    
    # Intersect all VCFs
    intersections = {}
    for (lhs_label, lhs_vcf), (rhs_label, rhs_vcf) in itertools.combinations(labelled_vcfs, 2):
        if lhs_label not in intersections:
            intersections[lhs_label] = {}
        intersections[lhs_label][rhs_label] = intersect(lhs_label, lhs_vcf, rhs_label, rhs_vcf, args)
    
    for i in range(len(vcfs)):
        isecs = [intersections[labels[i]][rhs_label][0] for rhs_label in labels[i + 1:]]
        tail_len = len(vcfs) - (i + 1)
        for tail_mask in itertools.product([0, 1], repeat=tail_len):
            if any(tail_mask):
                hot_labels = labels[i] + ''.join([labels[j + i + 1] for j in range(tail_len) if tail_mask[j]])
                out_vcf = join(args.output, hot_labels + '.vcf.gz')
                
                if i > 0:
                    tail_mask = list(tail_mask)
                    # Need to remove matches with previous VCFs
                    for head in labels[:i]:
                        isecs.append(intersections[head][labels[i]][1])
                        tail_mask.append(0)
                
                naive_intersect(vcfs[i], isecs, tail_mask, out_vcf)
                index_vcf(out_vcf)
    
    # Find unique records
    for i, label in enumerate(labels):
        isecs = []
        if label in intersections:
            isecs += [isec[2] for _, isec in intersections[label].items()]
        isecs += [intersections[lhs_label][label][3] for lhs_label in labels[:i]]
        out_vcf = join(args.output, label + '.vcf.gz')
        if len(isecs) > 1:
            naive_common(isecs, out_vcf)
            index_vcf(out_vcf)
        else:
            shutil.move(isecs[0], out_vcf)
            shutil.move(isecs[0] + '.tbi', out_vcf + '.tbi')
    
    with open(join(args.output, 'README.txt'), 'w') as readme:
        for label, vcf in labelled_vcfs:
            readme.write(label + ': ' + vcf + '\n')
    
    # Cleanup
    for lhs_label, rhs_isecs in intersections.items():
        for rhs_label, isecs in rhs_isecs.items():
            for isec in isecs:
                if exists(isec):
                    remove_vcf(isec)
    
    # Make plots
    if can_draw_venn and args.names is not None:
        if len(vcfs) < 7:
            venn_points = [[] for _ in range(len(vcfs))]
            venn_labels = {}
            for mask in itertools.product([0, 1], repeat=len(vcfs)):
                if any(mask):
                    venn_key = ''.join([str(bit) for bit in mask])
                    hot_labels = ''.join([labels[i] for i, bit in enumerate(mask) if bit])
                    vcf = join(args.output, hot_labels + '.vcf.gz')
                    nrecs = count_records(vcf)
                    venn_labels[venn_key] = str(nrecs)
            fig, ax = plot_ven(venn_labels, args.names)
            ax.legend_.remove()
            plt.tight_layout()
            if args.vennout is None:
                plt.show()
            else:
                plt.savefig(args.vennout, format='pdf', )
        else:
            print("Venn plots only supported for up to 6 VCFs")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-V', '--variants',
                        nargs='+',
                        type=str,
                        required=True,
                        help='VCF files to intersect')
    parser.add_argument('--rtg', 
                        type=str,
                        required=True,
                        help='RTG Tools binary')
    parser.add_argument('--sdf',
                        type=str,
                        required=True,
                        help='RTG Tools SDF reference index')
    parser.add_argument('-O', '--output',
                        type=str,
                        required=True,
                        help='Output directory')
    parser.add_argument('--regions',
                        type=str,
                        required=False,
                        help='regions in BED format to perform intersection')
    parser.add_argument('--ignore_genotypes',
                        default=False,
                        action='store_true',
                        help='Ignore genotypes during intersection; use variants only')
    parser.add_argument('--all_records',
                        default=False,
                        action='store_true',
                        help='Intersect all records')
    parser.add_argument('--sample',
                        type=str,
                        required=False,
                        help='Sample to compare, only required if multiple samples in VCFs')
    parser.add_argument('--names',
                        nargs='+',
                        type=str,
                        required=False,
                        help='Display names of the VCF files')
    parser.add_argument('--vennout',
                        type=str,
                        required=False,
                        help='Save Venn diagram in PDF format')
    parser.add_argument('--threads',
                        type=int,
                        required=False,
                        help='Maximum number of threads to use (default is all cores)')
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
