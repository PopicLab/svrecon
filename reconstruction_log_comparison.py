from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple

# --- Configuration ---
# Tolerance for location comparison (in base pairs)
# If the start and end of two intervals are within this tolerance, they are considered a match.
LOCATION_TOLERANCE = 10


@dataclass(frozen=True)
class LogEntry:
    """Represents a structured entry from the log file."""
    sv_id: str
    sv_type: str
    chrom: str
    start: int
    end: int
    hit_status: str

    @property
    def location_key(self) -> str:
        """Returns a string key for grouping by type and chromosome."""
        return f"{self.sv_type}_{self.chrom}"


def parse_log_line(line: str) -> LogEntry | None:
    """Parses a single log line into a LogEntry object."""
    # Pattern to capture the fields: ID, Type, Chrom:Start-End, HitStatus, [Score]
    # Example: INFO:__main__:sv8277    DEL    chr7:151517579-151517778   hit    [0.009950248756218905]
    pattern = re.compile(
        r"INFO:__main__:(\S+)\s+"  # sv_id (group 1)
        r"(\S+)\s+"  # sv_type (group 2)
        r"(\w+):(\d+)-(\d+)\s+"  # chrom:start-end (groups 3, 4, 5)
        r"(\S+)\s+"  # hit_status (group 6)
        r"\[([\d.,\s]+)\]"  # [score] (group 7)
    )

    match = pattern.search(line)
    if match:
        try:
            sv_id, sv_type, chrom, start_str, end_str, hit_status, score_str = match.groups()
            return LogEntry(
                sv_id=sv_id,
                sv_type=sv_type,
                chrom=chrom,
                start=int(start_str),
                end=int(end_str),
                hit_status=hit_status,
            )
        except (ValueError, IndexError) as e:
            # Handle lines that match the pattern but have bad data types
            print(f"Skipping malformed data in line: {line.strip()} | Error: {e}", file=sys.stderr)
            return None
    return None


def load_log_file(filepath: str) -> Dict[str, List[LogEntry]]:
    """
    Loads a log file and groups entries by (SV Type, Chromosome).

    Returns: A dictionary where keys are 'TYPE_CHROM' and values are lists of LogEntry objects.
    """
    entries_by_key = {}

    try:
        with open(filepath, 'r') as f:
            for line in f:
                entry = parse_log_line(line)
                if entry:
                    key = entry.location_key
                    if key not in entries_by_key:
                        entries_by_key[key] = []
                    entries_by_key[key].append(entry)
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}", file=sys.stderr)
        sys.exit(1)

    return entries_by_key

def find_unique_locations(
        file1_entries: Dict[str, List[LogEntry]],
        file2_entries: Dict[str, List[LogEntry]],
        tolerance: int,
        skip: List,
) -> Tuple[Set[str], Set[str]]:
    """
    Compares two sets of log entries based on type and chromosome,
    using the given tolerance for start/end coordinates.

    Returns: A tuple (unique_in_file1, unique_in_file2) where each is a set of formatted location strings,
             now including the hit status.
    """
    unique_in_file1 = set()
    unique_in_file2 = set()

    # Combine all keys (Type_Chrom) from both files
    all_keys = set(file1_entries.keys()) | set(file2_entries.keys())

    for key in all_keys:
        # Get entries for the current key, defaulting to empty list if not present
        list1 = file1_entries.get(key, [])
        list2 = file2_entries.get(key, [])

        # Track which indices in each list have been matched
        matched1 = [False] * len(list1)
        matched2 = [False] * len(list2)

        # --- Matching Loop ---
        for i, entry1 in enumerate(list1):
            for j, entry2 in enumerate(list2):
                # Skip if entry2 has already been matched to a preceding entry1
                if matched2[j]:
                    continue

                # Check if start and end are within tolerance
                start_diff = abs(entry1.start - entry2.start)
                end_diff = abs(entry1.end - entry2.end)

                if start_diff <= tolerance and end_diff <= tolerance:
                    # Found a match
                    matched1[i] = True
                    matched2[j] = True
                    # Optimization: Since the entry1 is matched, break the inner loop
                    break

        # --- Collect Unmatched Entries ---
        # MODIFICATION: Include the hit_status in the unique location string
        for i, entry in enumerate(list1):
            if not matched1[i]:
                # Format: TYPE CHROM:START-END (STATUS)
                if entry.sv_type not in skip:
                    unique_in_file1.add(f"{entry.sv_type} {entry.chrom}:{entry.start}-{entry.end} ({entry.hit_status})")

        for j, entry in enumerate(list2):
            if not matched2[j]:
                # Format: TYPE CHROM:START-END (STATUS)
                if entry.sv_type not in skip:
                    unique_in_file2.add(f"{entry.sv_type} {entry.chrom}:{entry.start}-{entry.end} ({entry.hit_status})")

    return unique_in_file1, unique_in_file2

# --- Main Execution ---
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <file1_path> <file2_path>", file=sys.stderr)
        sys.exit(1)

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]

    print(f"--- Loading and Parsing Files (Tolerance: {LOCATION_TOLERANCE} bp) ---")
    file1_data = load_log_file(file1_path)
    file2_data = load_log_file(file2_path)

    print(f"File 1 unique keys (Type_Chrom): {len(file1_data)}")
    print(f"File 2 unique keys (Type_Chrom): {len(file2_data)}")

    print("\n--- Comparing Locations ---")
    unique_1, unique_2 = find_unique_locations(file1_data, file2_data, LOCATION_TOLERANCE, skip=['DEL', 'INV', 'DUP'])

    # --- Output Results ---
    # The output section remains clean and simply prints the new, comprehensive strings.

    print(f"\n=======================================================")
    print(f"Unique locations in {file1_path} (not in {file2_path}): {len(unique_1)}")
    print("=======================================================")
    if unique_1:
        for location in sorted(unique_1):
            print(location)
    else:
        print("None found.")

    print(f"\n=======================================================")
    print(f"Unique locations in {file2_path} (not in {file1_path}): {len(unique_2)}")
    print("=======================================================")
    if unique_2:
        for location in sorted(unique_2):
            print(location)
    else:
        print("None found.")