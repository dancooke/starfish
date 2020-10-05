#!/usr/bin/env python3

import argparse
import itertools
import shutil
import os
import string
import sys
import subprocess as sp
import pysam as ps
from os.path import join, basename, exists, dirname

can_draw_venn = True

try:
    import matplotlib
    matplotlib.use('Agg')
    import pyvenn.venn as pv
    import matplotlib.pyplot as plt
except ImportError:
    can_draw_venn = False

bcftools = 'bcftools'
tabix = 'tabix'
rtg = 'rtg'

def call(cmd, print_cmd=False, quite=True):
    if print_cmd:
        print(' '.join(cmd))
    if quite:
        with open(os.devnull, "w") as f:
            sp.call(cmd, stdout=f, stderr=f)
    else:
        sp.call(cmd)

def vcf_index_exists(vcf_fname):
    return exists(vcf_fname + '.tbi')

def index_vcf(vcf_fname, overwrite=True):
    if overwrite:
        call([tabix, '-f', vcf_fname])
    else:
        call([tabix, vcf_fname])

def remove_vcf_index(vcf_fname):
    os.remove(vcf_fname + '.tbi')

def remove_vcf(vcf_fname):
    os.remove(vcf_fname)
    if vcf_index_exists(vcf_fname):
        remove_vcf_index(vcf_fname)

def count_records(vcf_fname):
    vcf = ps.VariantFile(vcf_fname, 'r')
    return sum(1 for rec in vcf)

def run_rtg(ref_sdf, lhs_vcf, rhs_vcf, out_dir,
            bed_regions=None, all_records=False,
            ref_overlap=False, squash_ploidy=False,
            ignore_genotypes=False, sample=None,
            score_field='QUAL',
            ploidy=None,
            threads=None,
            debug=False):
    cmd = [rtg, 'vcfeval', '-t', ref_sdf, '-b', lhs_vcf, '-c', rhs_vcf, '-o', out_dir]
    if bed_regions is not None:
        cmd += ['--bed-regions', bed_regions]
    if all_records:
        cmd.append('--all-records')
    if ref_overlap:
        cmd.append('--ref-overlap')
    if squash_ploidy:
        cmd.append('--squash-ploidy')
    if ignore_genotypes:
        cmd += ['--sample', 'ALT']
    elif sample is not None:
        cmd += ['--sample', sample]
    if threads is not None:
        cmd += ['--threads', str(threads)]
    if ploidy is not None:
        cmd += ['--Xdefault-ploidy', str(ploidy)]
    cmd += ['--vcf-score-field', score_field]
    call(cmd, print_cmd=debug, quite=not debug)

def make_empty_vcf(vcf_fname, template_vcf_fname, index=True):
    bcftools_cmd = [bcftools, 'view', '-h', '-Oz', '-o', vcf_fname, template_vcf_fname]
    call(bcftools_cmd)
    if index:
        index_vcf(vcf_fname)

def rtg_intersect(lhs_label, lhs_vcf, rhs_label, rhs_vcf, args, debug=False):
    tmp_dir = join(args.output, 'temp')
    run_rtg(args.sdf, lhs_vcf, rhs_vcf, tmp_dir,
            bed_regions=args.regions,
            all_records=args.all_records,
            ref_overlap=args.ref_overlap,
            squash_ploidy=args.squash_ploidy,
            sample=args.sample,
            ploidy=args.ploidy,
            threads=args.threads,
            debug=debug)
    lhs_and_rhs = join(args.output, lhs_label + '_and_' + rhs_label + '.vcf.gz')
    rhs_and_lhs = join(args.output, rhs_label + '_and_' + lhs_label + '.vcf.gz')
    lhs_not_rhs = join(args.output, lhs_label + '_not_' + rhs_label + '.vcf.gz')
    rhs_not_lhs = join(args.output, rhs_label + '_not_' + lhs_label + '.vcf.gz')
    tp_baseline = join(tmp_dir, 'tp-baseline.vcf.gz')
    tp = join(tmp_dir, 'tp.vcf.gz')
    fn = join(tmp_dir, 'fn.vcf.gz')
    fp = join(tmp_dir, 'fp.vcf.gz')
    if exists(tp_baseline):
        shutil.move(tp_baseline, lhs_and_rhs)
        shutil.move(tp_baseline + '.tbi', lhs_and_rhs + '.tbi')
        shutil.move(tp, rhs_and_lhs)
        shutil.move(tp + '.tbi', rhs_and_lhs + '.tbi')
        shutil.move(fn, lhs_not_rhs)
        shutil.move(fn + '.tbi', lhs_not_rhs + '.tbi')
        shutil.move(fp, rhs_not_lhs)
        shutil.move(fp + '.tbi', rhs_not_lhs + '.tbi')
    else:
        # then there were no sequences in common (or one/both of the input are empty)
        make_empty_vcf(lhs_and_rhs, lhs_vcf)
        make_empty_vcf(rhs_and_lhs, rhs_vcf)
        shutil.copyfile(lhs_vcf, lhs_not_rhs)
        shutil.copyfile(lhs_vcf + '.tbi', lhs_not_rhs + '.tbi')
        shutil.copyfile(rhs_vcf, rhs_not_lhs)
        shutil.copyfile(rhs_vcf + '.tbi', rhs_not_lhs + '.tbi')
    shutil.rmtree(tmp_dir)
    return lhs_and_rhs, rhs_and_lhs, lhs_not_rhs, rhs_not_lhs

def naive_common(vcfs, out, debug=False):
    assert len(vcfs) > 1
    cmd = [bcftools, 'isec', '-Oz', '-o', out, '-n=' + str(len(vcfs)), '-w1'] + vcfs
    call(cmd, print_cmd=debug, quite=not debug)

def naive_complement(head_vcf, tail_vcfs, out, debug=False):
    assert len(tail_vcfs) > 0
    cmd = [bcftools, 'isec', '-C', '-w1', '-Oz', '-o', out, head_vcf] + tail_vcfs
    call(cmd, print_cmd=debug, quite=not debug)

def naive_intersect(head_vcf, tail_vcfs, tail_mask, out, debug=False):
    positive_tails = [tail_vcfs[i] for i, mask in enumerate(tail_mask) if mask]
    assert len(positive_tails) > 0
    negative_tails = [tail_vcfs[i] for i, mask in enumerate(tail_mask) if not mask]
    if len(negative_tails) > 0:
        common_positive_temp_vcf = out + '.temp'
        naive_common([head_vcf] + positive_tails, common_positive_temp_vcf, debug=debug)
        index_vcf(common_positive_temp_vcf)
        naive_complement(common_positive_temp_vcf, negative_tails, out, debug=debug)
        remove_vcf(common_positive_temp_vcf)
    else:
        naive_common([head_vcf] + positive_tails, out, debug=debug)

def concat(vcfs, out, remove_duplicates=True, debug=False):
    assert len(vcfs) > 1
    cmd = [bcftools, 'concat', '-a', '-Oz', '-o', out]
    if remove_duplicates:
        cmd.append('-D')
    cmd += vcfs
    call(cmd, print_cmd=debug, quite=not debug)

def merge(vcfs, out, force_samples=True, debug=False):
    assert len(vcfs) > 1
    cmd = [bcftools, 'merge', '-Oz', '-o', out]
    if force_samples:
        cmd.append('--force-samples')
    cmd += vcfs
    call(cmd, print_cmd=debug, quite=not debug)

def read_samples(vcf_filename):
    vcf = ps.VariantFile(vcf_filename)
    return list(vcf.header.samples)

def can_concat(vcfs):
    samples = [read_samples(vcf) for vcf in vcfs]
    return samples[1:] == samples[:-1]

def combine(vcfs, out, debug=False):
    if can_concat(vcfs):
        concat(vcfs, out, debug=debug)
    else:
        merge(vcfs, out, debug=debug)

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
    
    global rtg
    rtg = args.rtg
    global bcftools
    bcftools = args.bcftools
    
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
        intersections[lhs_label][rhs_label] = rtg_intersect(lhs_label, lhs_vcf, rhs_label, rhs_vcf, args, debug=args.verbose)
    
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
                naive_intersect(vcfs[i], isecs, tail_mask, out_vcf, debug=args.verbose)
                index_vcf(out_vcf)
        
    # Find unique records
    for i, label in enumerate(labels):
        isecs = []
        if label in intersections:
            isecs += [isec[2] for _, isec in intersections[label].items()]
        isecs += [intersections[lhs_label][label][3] for lhs_label in labels[:i]]
        out_vcf = join(args.output, label + '.vcf.gz')
        if len(isecs) > 1:
            naive_common(isecs, out_vcf, debug=args.verbose)
            index_vcf(out_vcf)
        else:
            shutil.move(isecs[0], out_vcf)
            shutil.move(isecs[0] + '.tbi', out_vcf + '.tbi')
    
    # Write README
    with open(join(args.output, 'README.txt'), 'w') as readme:
        readme.write("Command: " + ' '.join(sys.argv) + '\n\n')
        for label, vcf in labelled_vcfs:
            readme.write(label + ': ' + vcf + '\n')
    
    # Cleanup
    for lhs_label, rhs_isecs in intersections.items():
        for rhs_label, isecs in rhs_isecs.items():
            for isec in isecs:
                if exists(isec):
                    remove_vcf(isec)
    
    if len(vcfs) > 2:
        partial_supported_vcfs = [[] for _ in range(len(vcfs) + 1)]
        for mask in itertools.product([0, 1], repeat=len(vcfs)):
            if any(mask):
                hot_labels = ''.join([labels[i] for i, bit in enumerate(mask) if bit])
                vcf_fname = join(args.output, hot_labels + '.vcf.gz')
                partial_supported_vcfs[len(hot_labels)].append(vcf_fname)
        for i, ivcfs in enumerate(partial_supported_vcfs[2:-1], 2):
            ivcf = join(args.output, str(i) + '.vcf.gz')
            combine(ivcfs, ivcf, debug=args.verbose)
            index_vcf(ivcf)
            iplus_vcfs = ivcfs + [vcf for vcfs in partial_supported_vcfs[i:] for vcf in vcfs]
            iplus_vcf = join(args.output, str(i) + '+.vcf.gz')
            combine(iplus_vcfs, iplus_vcf, debug=args.verbose)
            index_vcf(iplus_vcf)
            iminus_vcfs = ivcfs + [vcf for vcfs in partial_supported_vcfs[:i] for vcf in vcfs]
            iminus_vcf = join(args.output, str(i) + '-.vcf.gz')
            combine(iminus_vcfs, iminus_vcf, debug=args.verbose)
            index_vcf(iminus_vcf)
    
    # Make plots
    if can_draw_venn and args.vennout is not None:
        names = args.names if args.names is not None else labels
        if len(vcfs) < 7:
            venn_points = [[] for _ in range(len(vcfs))]
            venn_labels = {}
            for mask in itertools.product([0, 1], repeat=len(vcfs)):
                if any(mask):
                    venn_key = ''.join([str(bit) for bit in mask])
                    hot_labels = ''.join([labels[i] for i, bit in enumerate(mask) if bit])
                    vcf_fname = join(args.output, hot_labels + '.vcf.gz')
                    nrecs = count_records(vcf_fname)
                    venn_labels[venn_key] = str(nrecs)
            fig, ax = plot_ven(venn_labels, names)
            ax.legend_.remove()
            plt.tight_layout()
            if args.vennout is None:
                plt.show()
            else:
                if len(dirname(args.vennout)) == 0:
                    args.vennout = join(args.output, args.vennout)
                plt.savefig(args.vennout, format=args.vennout_format)
        else:
            print("Venn plots only supported for up to 6 VCFs")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-V', '--variants',
                        nargs='+',
                        type=str,
                        required=True,
                        help='VCF files to intersect')
    parser.add_argument('-t', '--sdf',
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
    parser.add_argument('--squash-ploidy',
                        default=False,
                        action='store_true',
                        help='Perform haplploid matching - ignore genotype mismatches')
    parser.add_argument('--sample',
                        type=str,
                        required=False,
                        help='Sample to compare (if multiple samples in VCFs) or ALT to ignore genotypes')
    parser.add_argument('--all-records',
                        default=False,
                        action='store_true',
                        help='Intersect all records')
    parser.add_argument('--ploidy', 
                        type=int,
                        required=False,
                        help='Set the default ploidy for comparison')
    parser.add_argument('--names',
                        nargs='+',
                        type=str,
                        required=False,
                        help='Display names of the VCF files')
    parser.add_argument('--vennout',
                        type=str,
                        required=False,
                        help='Save Venn diagram to file')
    parser.add_argument('--vennout-format',
                        type=str,
                        default='pdf',
                        help='File format to save Venn diagram')
    parser.add_argument('--threads',
                        type=int,
                        required=False,
                        help='Maximum number of threads to use (default is all cores)')
    parser.add_argument('--ref-overlap',
                        default=False,
                        action='store_true',
                        help='Call RTG vcfeval with "ref-overlap" option')
    parser.add_argument('--rtg', 
                        type=str,
                        default='rtg',
                        help='RTG Tools binary')
    parser.add_argument('--bcftools', 
                        type=str,
                        default='bcftools',
                        help='bcftools binary')
    parser.add_argument('--verbose',
                        default=False,
                        action='store_true',
                        help='Print executed commands and command output')
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
