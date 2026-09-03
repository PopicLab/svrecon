"""One test per CigarQueryCheck, over CIGARs of COLUMNS aligned columns built from a column
string. Every check reads only one of its two arguments -- the CIGAR or the query -- so the
unused one is passed as None."""
import random
from itertools import groupby

import pytest

from svrecon.config import DEFAULTS
from svrecon.constants import QueryValidationStatus
from svrecon.reconstruct import Query
from svrecon.scorers.base import (CigarQueryAlignmentSimilarityCheck, CigarQueryMappableCheck,
                                  CigarQueryMaximumErrorWindowCheck, CigarQueryNoLargeErrorsCheck)
from svrecon.scorers.utils import Cigar

NUM_COLUMNS = 1000  # columns in all test CIGARs
NUM_CLIP = 50
ERROR_THRESHOLD = DEFAULTS['match_error_threshold']
WINDOW, WINDOW_ERROR = DEFAULTS['max_window_size'], DEFAULTS['max_window_error']
MIN_MAPPABLE, MAX_UNMAPPABLE = DEFAULTS['min_mappable_fraction'], DEFAULTS['max_unmappable_size']
MAX_CONTIGUOUS_ERROR = 50
NM_OPS = 'XID'              # all Cigar.NM counts
COLUMN_ERROR_OPS = 'XIDN'  # error columns to the window check, which counts target gaps too

# Helper functions for cigar scoring
def string_to_cigar(columns: str) -> Cigar:
    return Cigar([(Cigar._OP_CHARS.index(op), len(list(run))) for op, run in groupby(columns)])


def pad_matches(run: str = '') -> str:
    """``run`` centred in COLUMNS aligned columns, padded out with matches."""
    pad = NUM_COLUMNS - len(run)
    return 'M' * (pad // 2) + run + 'M' * (pad - pad // 2)


def get_dummy_query(sequence: str) -> Query:
    """Only ``sequence`` is ever read off it."""
    return Query(chrom='chrT', svtype='sv', svid='sv', grammar='dummy', sequence=sequence,
                 ref_start=0, ref_end=0, recon_segments=[], ref_segments=[], ref_sequence='',
                 buffer=0)


def pad_bases(seq: str = '') -> str:
    """returned padded seq with 'A'"""
    pad = NUM_COLUMNS - len(seq)
    return 'A' * (pad // 2) + seq + 'A' * (pad - pad // 2)


def spread_unmappable(mappable_fraction: float, seed: int = 0) -> str:
    """NUM_COLUMNS bases with ``mappable_fraction`` of them A, the N's scattered at random. Too
    few to bunch into a run anywhere near MAX_UNMAPPABLE, so only the proportion trips a check."""
    num_unmappable = NUM_COLUMNS - int(NUM_COLUMNS * mappable_fraction)
    positions = set(random.Random(seed).sample(range(NUM_COLUMNS), num_unmappable))
    return ''.join('N' if column in positions else 'A' for column in range(NUM_COLUMNS))


SIMILARITY_CHECK = CigarQueryAlignmentSimilarityCheck(match_error_threshold=ERROR_THRESHOLD)
NO_LARGE_ERRORS_CHECK = CigarQueryNoLargeErrorsCheck(max_size=MAX_CONTIGUOUS_ERROR)
ERROR_WINDOW_CHECK = CigarQueryMaximumErrorWindowCheck(max_window_size=WINDOW,
                                                       max_window_error=WINDOW_ERROR)
MAPPABLE_CHECK = CigarQueryMappableCheck(min_mappable_fraction=MIN_MAPPABLE,
                                         max_unmappable_size=MAX_UNMAPPABLE)

AT_THRESHOLD = int(NUM_COLUMNS * ERROR_THRESHOLD)  # error columns the check still passes on (it is <=)
UNDER = MAX_CONTIGUOUS_ERROR - 1


@pytest.mark.parametrize('errors,expected', [(AT_THRESHOLD, QueryValidationStatus.PASS),
                                             (AT_THRESHOLD + 1, QueryValidationStatus.FAIL)],
                         ids=['at-threshold', 'over-threshold'])
@pytest.mark.parametrize('seed', range(10))
def test_alignment_similarity(seed, errors, expected):
    """Tests alignment similarities at or above error thresholds"""
    errors = ''.join(random.Random(seed).choices(NM_OPS, k=errors))
    columns = 'S' * NUM_CLIP + pad_matches(errors) + 'S' * NUM_CLIP
    assert SIMILARITY_CHECK.validate(string_to_cigar(columns), None).status is expected


@pytest.mark.parametrize('columns,expected', [
    (pad_matches('I' * MAX_CONTIGUOUS_ERROR), QueryValidationStatus.FAIL),
    (pad_matches('D' * MAX_CONTIGUOUS_ERROR), QueryValidationStatus.FAIL),
    (pad_matches('N' * MAX_CONTIGUOUS_ERROR), QueryValidationStatus.FAIL),
    ('S' * MAX_CONTIGUOUS_ERROR + pad_matches(), QueryValidationStatus.FAIL),
    (pad_matches() + 'S' * MAX_CONTIGUOUS_ERROR, QueryValidationStatus.FAIL),
    (pad_matches('X' * MAX_CONTIGUOUS_ERROR), QueryValidationStatus.PASS),
    ('S' * UNDER + pad_matches('I' * UNDER + 'M' * UNDER + 'D' * UNDER + 'M' * UNDER + 'N' * UNDER)
     + 'S' * UNDER, QueryValidationStatus.PASS),
], ids=['insertion', 'deletion', 'ref-skip', 'left-clip', 'right-clip', 'mismatches-not-a-run',
        'all-under'])
def test_no_large_errors(columns, expected):
    """Tests large indels/clips directly at or under error thresholds"""
    assert NO_LARGE_ERRORS_CHECK.validate(string_to_cigar(columns), None).status is expected


@pytest.mark.parametrize('errors,expected', [(WINDOW_ERROR, QueryValidationStatus.FAIL),
                                             (WINDOW_ERROR - 1, QueryValidationStatus.PASS)],
                         ids=['at-limit', 'under-limit'])
@pytest.mark.parametrize('seed', range(10))
def test_maximum_error_window(seed, errors, expected):
    """tests columns with errors scattered inside one random window."""
    rng = random.Random(seed)
    columns = list(pad_matches())
    start = rng.randrange(NUM_COLUMNS - WINDOW + 1)
    for column in rng.sample(range(start, start + WINDOW), errors):
        columns[column] = rng.choice(COLUMN_ERROR_OPS)
    assert ERROR_WINDOW_CHECK.validate(string_to_cigar(''.join(columns)), None).status is expected


@pytest.mark.parametrize('run_length,expected', [
    (0, QueryValidationStatus.PASS),
    (MAX_UNMAPPABLE - 1, QueryValidationStatus.PASS),
    (MAX_UNMAPPABLE, QueryValidationStatus.INCONCLUSIVE),
], ids=['no-ns', 'under-size', 'at-size'])
def test_mappable_run_size(run_length, expected):
    """sequences of varying lengths of N blocks"""
    sequence = pad_bases('N' * run_length)
    assert MAPPABLE_CHECK.validate(None, get_dummy_query(sequence)).status is expected


@pytest.mark.parametrize('mappable_fraction,expected', [
    (MIN_MAPPABLE, QueryValidationStatus.PASS),
    (MIN_MAPPABLE - 1 / NUM_COLUMNS, QueryValidationStatus.INCONCLUSIVE),
], ids=['at-threshold', 'one-base-under'])
def test_mappable_fraction(mappable_fraction, expected):
    """sequences of varying proportions sequence consisting of N"""
    sequence = spread_unmappable(mappable_fraction)
    assert MAPPABLE_CHECK.validate(None, get_dummy_query(sequence)).status is expected

@pytest.mark.parametrize('fraction,size,sequence', [
    (None, MAX_UNMAPPABLE, spread_unmappable(MIN_MAPPABLE - 1 / NUM_COLUMNS)),
    (MIN_MAPPABLE, None, pad_bases('N' * MAX_UNMAPPABLE)),
    (None, None, pad_bases('N' * (NUM_COLUMNS // 2))),
], ids=['fraction-off', 'run-off', 'both-off'])
def test_mappable_conditions_are_individually_optional(fraction, size, sequence):
    """Every sequence here trips whichever condition is switched off, so a None threshold passes."""
    check = CigarQueryMappableCheck(min_mappable_fraction=fraction, max_unmappable_size=size)
    assert check.validate(None, get_dummy_query(sequence)).status is QueryValidationStatus.PASS
