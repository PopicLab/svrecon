from typing import List, TYPE_CHECKING

import pysam

if TYPE_CHECKING:
    from svrecon.scoring import SVValidationResult


def write_annotated_vcf(results: List['SVValidationResult'], vcf_path: str) -> None:
    if not results:
        return
    header = results[0].records[0].header.copy()
    header.add_line('##INFO=<ID=SVRECON,Number=1,Type=String,'
                     'Description="svrecon evaluation outcome: hit, miss, or inconclusive">')
    if any(sv_validation.grammars for sv_validation in results):
        header.add_line('##INFO=<ID=SVRECON_GRAMMAR,Number=1,Type=String,'
                        'Description="svrecon outcome of this record\'s GRAMMAR_ID; '
                        'SVRECON is the outcome of the SV">')

    with pysam.VariantFile(vcf_path, 'w', header=header) as vcf_out:
        for sv_validation in results:
            for record in sv_validation.records:
                record.translate(header)
                record.info['SVRECON'] = sv_validation.outcome.value
                if sv_validation.grammars:
                    record.info['SVRECON_GRAMMAR'] = sv_validation.grammars[record.info['GRAMMAR_ID']].outcome.value
                vcf_out.write(record)
