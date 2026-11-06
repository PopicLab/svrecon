"""Grouping of one call's alternative grammars by GRAMMAR_ID, and the event roll-up over them."""
import pytest

from svrecon.constants import Outcome, QueryValidationStatus, ValidationSource
from svrecon.scorers.base import QueryValidation
from svrecon.scoring import QueryValidationResult, SVValidationResult
from svrecon.utils import group_records_by_grammar_id


class Record:
    def __init__(self, grammar_id=None):
        self.info = {'SVID': 'sv0'}
        if grammar_id is not None:
            self.info['GRAMMAR_ID'] = grammar_id


def grammar_result(grammar_id, status, source=ValidationSource.READS):
    validation = QueryValidation(source=source, passed=status is QueryValidationStatus.PASS, status=status)
    result = SVValidationResult('sv0', 'delINVdup', [QueryValidationResult(query=None, validations=[validation])], [])
    result.grammar_id = grammar_id
    return result


def test_untagged_records_are_one_grammar():
    assert group_records_by_grammar_id([Record(), Record()]) is None


def test_records_group_by_grammar_id_in_order():
    records = [Record(1), Record(0), Record(1), Record(0)]
    grouped = group_records_by_grammar_id(records)
    assert list(grouped) == [0, 1]
    assert grouped[0] == [records[1], records[3]] and grouped[1] == [records[0], records[2]]


@pytest.mark.parametrize('statuses,outcome,grammar_id', [
    ((QueryValidationStatus.FAIL, QueryValidationStatus.PASS), Outcome.HIT, 1),
    ((QueryValidationStatus.PASS, QueryValidationStatus.PASS), Outcome.HIT, 0),
    ((QueryValidationStatus.FAIL, QueryValidationStatus.INCONCLUSIVE), Outcome.INCONCLUSIVE, 1),
    ((QueryValidationStatus.FAIL, QueryValidationStatus.FAIL), Outcome.MISS, 0),
])
def test_sv_takes_best_grammar(statuses, outcome, grammar_id):
    grammars = {i: grammar_result(i, status) for i, status in enumerate(statuses)}
    event = SVValidationResult.from_grammars('sv0', [], grammars)
    assert event.outcome == outcome
    assert event.grammar_id == grammar_id
    assert event.query_validation_results is grammars[grammar_id].query_validation_results


def test_read_hit_preferred_over_assembly_hit():
    grammars = {0: grammar_result(0, QueryValidationStatus.PASS, ValidationSource.ASSEMBLY),
                1: grammar_result(1, QueryValidationStatus.PASS, ValidationSource.READS)}
    event = SVValidationResult.from_grammars('sv0', [], grammars)
    assert event.grammar_id == 1
    assert event.validation_source == ValidationSource.READS
