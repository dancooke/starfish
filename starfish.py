#!/usr/bin/env python3

import argparse
import itertools
import shutil
import os
import string
import sys
import subprocess as sp
import pysam as ps
from pathlib import Path

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
        print(' '.join(str(c) for c in cmd))
    if quite:
        with open(os.devnull, "w") as f:
            sp.call(cmd, stdout=f, stderr=f)
    else:
        sp.call(cmd)

def vcf_index_exists(vcf_filename):
    vcf_index_filaneme = vcf_filename.with_suffix(vcf_filename.suffix + '.tbi')
    return vcf_index_filaneme.exists()

def index_vcf(vcf_filename, overwrite=True):
    if overwrite:
        call([tabix, '-f', vcf_filename])
    else:
        call([tabix, vcf_filename])

def remove_vcf_index(vcf_filename):
    vcf_index_filaneme = vcf_filename.with_suffix(vcf_filename.suffix + '.tbi')
    vcf_index_filaneme.unlink()

def remove_vcf(vcf_filename, index=True):
    vcf_filename.unlink()
    if index and vcf_index_exists(vcf_filename):
        remove_vcf_index(vcf_filename)

def decompose_multiallelic(source_vcf, dest_vcf, debug=False):
    cmd = [bcftools, 'norm', '-m', '-', '--force', '-Oz', '-o', dest_vcf, source_vcf]
    call(cmd, print_cmd=debug, quite=not debug)

def count_records(vcf_filename):
    vcf = ps.VariantFile(vcf_filename, 'r')
    return sum(1 for rec in vcf)

def run_vcfeval(ref_sdf, lhs_vcf, rhs_vcf, out_dir,
                bed_regions=None, 
                all_records=False,
                ref_overlap=False, 
                squash_ploidy=False,
                sample=None,
                score_field='GQ',
                decompose=False,
                ploidy=None,
                output_mode=None,
                flag_alternates=False,
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
    if sample is not None:
        cmd += ['--sample', sample]
    if decompose:
        cmd.append('--decompose')
    if threads is not None:
        cmd += ['--threads', str(threads)]
    if ploidy is not None:
        cmd += ['--Xdefault-ploidy', str(ploidy)]
    if output_mode is not None:
        cmd += ['--output-mode', output_mode]
    if flag_alternates:
        cmd.append("--XXcom.rtg.vcf.eval.flag-alternates=true")
    cmd += ['--vcf-score-field', score_field]
    call(cmd, print_cmd=debug, quite=not debug)

def split_baseline_annotated(baseline_vcf_filename, tp_vcf_filename, fn_vcf_filename):
    baseline_vcf = ps.VariantFile(baseline_vcf_filename)
    tp_vcf = ps.VariantFile(tp_vcf_filename, 'wz', header=baseline_vcf.header)
    fn_vcf = ps.VariantFile(fn_vcf_filename, 'wz', header=baseline_vcf.header)
    for rec in baseline_vcf:
        if rec.info["BASE"] == "TP":
            tp_vcf.write(rec)
        elif rec.info["BASE"] == "FN":
            if "BASE_ALTERNATE" in rec.info:
                tp_vcf.write(rec)
            else:
                fn_vcf.write(rec)
    tp_vcf.close()
    index_vcf(tp_vcf_filename)
    fn_vcf.close()
    index_vcf(fn_vcf_filename)

def split_calls_annotated(calls_vcf_filename, tp_vcf_filename, fp_vcf_filename):
    calls_vcf = ps.VariantFile(calls_vcf_filename)
    tp_vcf = ps.VariantFile(tp_vcf_filename, 'wz', header=calls_vcf.header)
    fp_vcf = ps.VariantFile(fp_vcf_filename, 'wz', header=calls_vcf.header)
    for rec in calls_vcf:
        if rec.info["CALL"] == "TP":
            tp_vcf.write(rec)
        elif rec.info["CALL"] == "FP":
            if "CALL_ALTERNATE" in rec.info:
                tp_vcf.write(rec)
            else:
                fp_vcf.write(rec)
    tp_vcf.close()
    index_vcf(tp_vcf_filename)
    fp_vcf.close()
    index_vcf(fp_vcf_filename)

def rename_vcf(source, dest):
    source.rename(dest)
    source_index = source.with_suffix(source.suffix + '.tbi')
    if source_index.exists():
        dest_index = dest.with_suffix(dest.suffix + '.tbi')
        source_index.rename(dest_index)

def copy_vcf(source, dest):
    shutil.copyfile(source, dest)
    source_index = source.with_suffix(source.suffix + '.tbi')
    if source_index.exists():
        dest_index = dest.with_suffix(dest.suffix + '.tbi')
        shutil.copyfile(source_index, dest_index)

def make_empty_vcf(vcf_filename, template_vcf_filename, index=True):
    bcftools_cmd = [bcftools, 'view', '-h', '-Oz', '-o', vcf_filename, template_vcf_filename]
    call(bcftools_cmd)
    if index:
        index_vcf(vcf_filename)

def rtg_intersect(lhs_label, lhs_vcf, rhs_label, rhs_vcf, args, debug=False):
    tmp_dir = args.output / 'temp'
    annotate_and_flag_alternative = args.squash_ploidy and (args.ploidy is not None and args.ploidy > 2)
    run_vcfeval(args.sdf, lhs_vcf, rhs_vcf, tmp_dir,
                bed_regions=args.regions,
                all_records=args.all_records,
                ref_overlap=args.ref_overlap,
                squash_ploidy=args.squash_ploidy,
                sample=args.sample,
                decompose=args.decompose,
                ploidy=args.ploidy,
                output_mode="annotate" if annotate_and_flag_alternative else None,
                flag_alternates=annotate_and_flag_alternative,
                threads=args.threads,
                debug=debug)
    lhs_and_rhs = args.output / (lhs_label + '_and_' + rhs_label + '.vcf.gz')
    rhs_and_lhs = args.output / (rhs_label + '_and_' + lhs_label + '.vcf.gz')
    lhs_not_rhs = args.output / (lhs_label + '_not_' + rhs_label + '.vcf.gz')
    rhs_not_lhs = args.output / (rhs_label + '_not_' + lhs_label + '.vcf.gz')
    tp_baseline = tmp_dir / 'tp-baseline.vcf.gz'
    tp = tmp_dir / 'tp.vcf.gz'
    fn = tmp_dir / 'fn.vcf.gz'
    fp = tmp_dir / 'fp.vcf.gz'
    if annotate_and_flag_alternative:
        split_baseline_annotated(tmp_dir / 'baseline.vcf.gz', tp_baseline, fn)
        split_calls_annotated(tmp_dir / 'calls.vcf.gz', tp, fp)
    if tp_baseline.exists():
        rename_vcf(tp_baseline, lhs_and_rhs)
        rename_vcf(tp, rhs_and_lhs)
        rename_vcf(fn, lhs_not_rhs)
        rename_vcf(fp, rhs_not_lhs)
    else:
        # then there were no sequences in common (or one/both of the input are empty)
        make_empty_vcf(lhs_and_rhs, lhs_vcf)
        make_empty_vcf(rhs_and_lhs, rhs_vcf)
        copy_vcf(lhs_vcf, lhs_not_rhs)
        copy_vcf(rhs_vcf, rhs_not_lhs)
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
        common_positive_temp_vcf = out.with_suffix(out.suffix + '.temp')
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
    
    if not args.output.exists():
        args.output.mkdir(parents=True)
    
    vcfs = args.variants
    labels = string.ascii_uppercase[:len(vcfs)]
    
    if len(vcfs) > len(labels):
        print("Maximum number of VCFs to intersect is", len(labels))
        return
    
    tmp_vcfs = []
    if args.decompose:
        for vcf in vcfs:
            tmp_vcf = args.output / vcf.stem.replace('.vcf', '.decompose.tmp.vcf.gz')
            decompose_multiallelic(vcf, tmp_vcf)
            index_vcf(tmp_vcf)
            tmp_vcfs.append(tmp_vcf)
        vcfs = tmp_vcfs
    
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
                out_vcf = args.output / (hot_labels + '.vcf.gz')
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
        out_vcf = args.output / (label + '.vcf.gz')
        if len(isecs) > 1:
            naive_common(isecs, out_vcf, debug=args.verbose)
            index_vcf(out_vcf)
        else:
            rename_vcf(isecs[0], out_vcf)
    
    # Write README
    with open(args.output / 'README.txt', 'w') as readme:
        readme.write("Command: " + ' '.join(sys.argv) + '\n\n')
        for label, vcf in labelled_vcfs:
            readme.write(label + ': ' + str(vcf) + '\n')
    
    # Cleanup
    for lhs_label, rhs_isecs in intersections.items():
        for rhs_label, isecs in rhs_isecs.items():
            for isec in isecs:
                if isec.exists():
                    remove_vcf(isec)
    for vcf in tmp_vcfs:
        remove_vcf(vcf)
    
    if len(vcfs) > 2:
        partial_supported_vcfs = [[] for _ in range(len(vcfs) + 1)]
        for mask in itertools.product([0, 1], repeat=len(vcfs)):
            if any(mask):
                hot_labels = ''.join([labels[i] for i, bit in enumerate(mask) if bit])
                vcf_fname = args.output / (hot_labels + '.vcf.gz')
                partial_supported_vcfs[len(hot_labels)].append(vcf_fname)
        for i, ivcfs in enumerate(partial_supported_vcfs[2:-1], 2):
            ivcf = args.output / (str(i) + '.vcf.gz')
            combine(ivcfs, ivcf, debug=args.verbose)
            index_vcf(ivcf)
            iplus_vcfs = ivcfs + [vcf for vcfs in partial_supported_vcfs[i:] for vcf in vcfs]
            iplus_vcf = args.output / (str(i) + '+.vcf.gz')
            combine(iplus_vcfs, iplus_vcf, debug=args.verbose)
            index_vcf(iplus_vcf)
            iminus_vcfs = ivcfs + [vcf for vcfs in partial_supported_vcfs[:i] for vcf in vcfs]
            iminus_vcf = args.output / (str(i) + '-.vcf.gz')
            combine(iminus_vcfs, iminus_vcf, debug=args.verbose)
            index_vcf(iminus_vcf)
    
    # Make plots
    if can_draw_venn and args.venn:
        names = args.names if args.names is not None else labels
        if len(vcfs) < 7:
            venn_points = [[] for _ in range(len(vcfs))]
            venn_labels = {}
            for mask in itertools.product([0, 1], repeat=len(vcfs)):
                if any(mask):
                    venn_key = ''.join([str(bit) for bit in mask])
                    hot_labels = ''.join([labels[i] for i, bit in enumerate(mask) if bit])
                    vcf_fname = args.output / (hot_labels + '.vcf.gz')
                    nrecs = count_records(vcf_fname)
                    venn_labels[venn_key] = str(nrecs)
            fig, ax = plot_ven(venn_labels, names)
            ax.legend_.remove()
            plt.tight_layout()
            plt.savefig(args.output / 'venn.pdf', format=args.vennout_format)
        else:
            print("Venn plots only supported for up to 6 VCFs")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-V', '--variants',
                        nargs='+',
                        type=Path,
                        required=True,
                        help='VCF files to intersect')
    parser.add_argument('-t', '--sdf',
                        type=str,
                        required=True,
                        help='RTG Tools SDF reference index')
    parser.add_argument('-O', '--output',
                        type=Path,
                        required=True,
                        help='Output directory')
    parser.add_argument('--regions',
                        type=Path,
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
    parser.add_argument('--venn',
                        default=False,
                        action='store_true',
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
    parser.add_argument('--decompose',
                        default=False,
                        action='store_true',
                        help='Decompose multi-allelic and complex alleles')
    parser.add_argument('--rtg', 
                        type=Path,
                        default='rtg',
                        help='RTG Tools binary')
    parser.add_argument('--bcftools', 
                        type=Path,
                        default='bcftools',
                        help='bcftools binary')
    parser.add_argument('--verbose',
                        default=False,
                        action='store_true',
                        help='Print executed commands and command output')
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
