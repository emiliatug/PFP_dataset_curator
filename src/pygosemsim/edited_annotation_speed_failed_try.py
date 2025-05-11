# import gzip
# from pathlib import Path
import gzip
from collections import defaultdict
from pathlib import Path

# resource_dir = Path(__file__).resolve().parent / "_resources"
resource_dir = Path(__file__).resolve().parent / "_resources"


def parse_gaf_line(line, qualified_only=True):
    """Parse a single GAF line and return relevant fields."""
    row = line.split("\t")
    qualifiers = row[3].split("|")
    if qualified_only and "NOT" in qualifiers:
        return None
    return {
        "db_object_id": row[1],
        "db_object_symbol": row[2],
        "db_object_name": row[9],
        "db_object_type": row[11],
        "go_id": row[4],
        "qualifiers": qualifiers,
        "evidence_code": row[6],
    }


def from_gaf_lines(lines, qualified_only=True):
    """Read gene association entries from an iterable of lines."""
    annots = defaultdict(
        lambda: {
            "db_object_symbol": None,
            "db_object_name": None,
            "db_object_type": None,
            "annotations": defaultdict(dict),
        }
    )

    format_ver = None
    for line in lines:
        if line.startswith("!"):
            if "gaf-version" in line:
                format_ver = line.split(":")[1].strip()
            continue
        parsed_line = parse_gaf_line(line, qualified_only)
        if not parsed_line:
            continue
        uid = parsed_line["db_object_id"]
        annots[uid]["db_object_symbol"] = parsed_line["db_object_symbol"]
        annots[uid]["db_object_name"] = parsed_line["db_object_name"]
        annots[uid]["db_object_type"] = parsed_line["db_object_type"]
        annots[uid]["annotations"][parsed_line["go_id"]] = {
            "go_id": parsed_line["go_id"],
            "qualifiers": parsed_line["qualifiers"],
            "evidence_code": parsed_line["evidence_code"],
        }

    print(f"gaf-version: {format_ver}")
    return annots


def from_gaf(pathlike, **kwargs):
    """Read gene association entries from a GAF file."""
    with gzip.open(pathlike, "rt") as f:
        annot = from_gaf_lines(f, **kwargs)
    return annot


# +
# def from_gaf_lines(lines, qualified_only=True):
#     """Read gene association entries
#     Reference:
#         http://www.geneontology.org/page/go-annotation-file-gaf-format-21
#     """
#     lines_iter = iter(lines)

#     # Header
#     fv_line = next(lines_iter)
#     format_ver = fv_line.split(":")[1].strip()
#     print(f"gaf-version: {format_ver}")
#     annots = {}
#     # Records
#     for line in lines_iter:
#         if line.startswith("!"):
#             continue
#         row = line.split("\t")
#         uid = row[1]  # DB Object ID (= UniProt ID)
#         if uid not in annots:
#             annots[uid] = {
#                 "db_object_id": uid,
#                 "db_object_symbol": row[2],
#                 "db_object_name": row[9],
#                 "db_object_type": row[11],
#                 "annotation": {}
#             }
#         qualifiers = row[3].split("|")
#         if qualified_only:
#             if "NOT" in qualifiers:
#                 continue
#         # Add GO annotation
#         go_id = row[4]
#         annots[uid]["annotation"][go_id] = {
#             "go_id": go_id,
#             "qualifier": qualifiers,
#             "evidence_code": row[6],
#         }
#     return annots

# +
# def from_gaf(pathlike, **kwargs):
#     with gzip.open(pathlike, "rt") as f:
#         annot = from_gaf_lines(f, **kwargs)
#     return annot

# +
# def from_resource(name, **kwargs):
#     filename = f"{name}.gaf.gz"
#     return from_gaf(resource_dir / filename, **kwargs)
