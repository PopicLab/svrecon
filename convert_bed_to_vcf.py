import argparse

import pysam
import math
from collections import defaultdict


def create_vcf_header():
    header = pysam.VariantHeader()
    for i in range(22):
        header.contigs.add("chr%d" % (i + 1))
    header.contigs.add("chrX")
    header.contigs.add("chrY")
    ##contig=<ID=chr1,length=248956422>
    header.info.add('SVTYPE', number=1, type='String', description='Type of structural variant')
    header.info.add('END', number=1, type='Integer', description='End position of the variant')
    header.info.add('SVLEN', number=1, type='Integer', description='Difference in length between REF and ALT alleles')
    header.info.add('CPXTYPE', number=1, type='String', description='CPX SV')
    header.add_line('##ALT=<ID=BND,Description="Breakend">')
    header.formats.add('GT', number=1, type='String', description='Genotype')
    header.add_sample('SAMPLE')
    return header


def convert_bed_record_to_variant(header, chrom, start, end, variant_type, qual, record_id):
    """
    Internal helper to convert a BED record to a samtools variant. This method retains zero-indexing, contrary
    to VCF standard, for better compatibility with Python.

    :param header:
    :param chrom:
    :param start:
    :param end:
    :param variant_type:
    :param qual:
    :param record_id:
    :return:
    """
    svtype = variant_type
    if variant_type == "AcbaC-A":
        svtype = "dupINVdup-A"
    # if variant_type == "nrTRA-A":
    #     svtype = "nrTRA-B"
    # if variant_type == "nrTRA-B":
    #     svtype = "nrTRA-A"
    # if variant_type == "INVdel-A":
    #     svtype = "INVdel-B"
    # elif variant_type == "INVdel-B":
    #     svtype = "INVdel-A"

    svtype = svtype.replace("-", "_")

    # pos = int(start)
    # end_pos = int(end)
    s, e = sorted([int(start), int(end)])
    # print(chrom, start, end, variant_type, record_id)
    svlen = e - s
    variant = header.new_record()
    variant.contig = chrom
    variant.start = s
    variant.stop = e
    variant.id = svtype
    variant.ref = 'N'
    variant.alts = ('<cpx>',)
    variant.qual = round(-10 * math.log10(1 - float(qual))) if float(qual) < 0.999999 else 60
    variant.filter.add('PASS')
    variant.info['SVTYPE'] = svtype
    variant.info['SVLEN'] = svlen
    variant.samples['SAMPLE']['GT'] = (None, None)  # ./.
    # print(variant)
    return variant


def convert_bed_to_vcf(input_file, output_dir):
    record_count = 0
    header = create_vcf_header()
    seen_rows = set()

    duplicated_lines = 0

    vcf_records = defaultdict(list)

    try:
        with open(input_file, 'r') as bed_file:
            for line_num, line in enumerate(bed_file, 1):
                line = line.strip()
                if not line or line.startswith('#'): continue
                fields = line.split('\t')
                conf_id, chrom, start, end, variant_type, qual = fields[0], fields[1], fields[2], fields[3], fields[4], fields[5]

                if (conf_id, chrom, min(start, end), max(start, end), variant_type) in seen_rows:
                    duplicated_lines += 1
                    continue
                seen_rows.add((conf_id, chrom, min(start, end), max(start, end), variant_type))

                variant = convert_bed_record_to_variant(header, chrom, start, end, variant_type, qual, record_count)

                vcf_records[conf_id].append(variant)

                record_count += 1
    except Exception as e:
        print(f"Error: {e}")
        raise
    print(f"Finished reading {input_file}. \nOmitted {duplicated_lines} duplicated lines out of {line_num} lines.")

    for conf_id, records in vcf_records.items():
        output_file = f"{output_dir}/{conf_id}.vcf"
        print(f"Writing {output_file}")
        with pysam.VariantFile(output_file, 'w', header=header) as vcf_out:
            for record in records:
                vcf_out.write(record)


    print(f"Finished conversion")



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", help="Input BED file")
    parser.add_argument("--output_dir", help="Output directory")
    args = parser.parse_args()
    convert_bed_to_vcf(args.input_file, args.output_dir)
