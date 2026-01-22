from collections import Counter

import gffutils


def summarize_db(db: gffutils.FeatureDB) -> dict:
    """
    Create a deterministic summary of a gffutils FeatureDB suitable for comparison.
    """
    summary = {}

    print("Basic metadata")
    summary["version"] = db.version
    summary["dialect"] = tuple(sorted(db.dialect.items()))

    print("Feature counts by type")
    type_counts = Counter(f.featuretype for f in db.all_features())
    summary["featuretype_counts"] = dict(sorted(type_counts.items()))

    print("Feature counts by (featuretype, strand)")
    type_strand_counts = Counter((f.featuretype, f.strand) for f in db.all_features())
    summary["featuretype_strand_counts"] = dict(sorted(type_strand_counts.items()))

    # Attributes seen per feature type
    # attrs_by_type = {}
    # for f in db.all_features():
    #     attrs_by_type.setdefault(f.featuretype, set()).update(f.attributes.keys())
    # summary["attributes_by_featuretype"] = {
    #     k: tuple(sorted(v)) for k, v in sorted(attrs_by_type.items())
    # }
    #
    # # Parent → child featuretype relationships
    # parent_child = Counter()
    # for f in db.all_features():
    #     for parent_id in f.attributes.get("Parent", []):
    #         parent_child[(parent_id, f.featuretype)] += 1
    # summary["parent_child_counts"] = dict(sorted(parent_child.items()))

    return summary


def compare_dbs(db1_path: str, db2_path: str) -> None:
    db1 = gffutils.FeatureDB(db1_path)
    db2 = gffutils.FeatureDB(db2_path)

    s1 = summarize_db(db1)
    s2 = summarize_db(db2)

    if s1 == s2:
        print("Databases are equivalent under summary comparison.")
        return

    print("Databases differ. Differences:")
    keys = sorted(set(s1) | set(s2))
    for k in keys:
        if s1.get(k) != s2.get(k):
            print(f"\n== {k} ==")
            print("DB1:", s1.get(k))
            print("DB2:", s2.get(k))


compare_dbs(
    db1_path="/tmp/gencode.v38.annotation.db",
    db2_path="/Users/dmccabe/Data/ref/gencode.v38.primary_assembly.annotation.db",
)

"""
== featuretype_counts ==
DB1: {'CDS': 825028, 'Selenocysteine': 119, 'UTR': 349402, 'exon': 1499012, 'gene': 60649, 'start_codon': 92997, 'stop_codon': 86205, 'transcript': 237012}
DB2: {'CDS': 825232, 'Selenocysteine': 119, 'UTR': 349477, 'exon': 1499267, 'gene': 60708, 'start_codon': 93028, 'stop_codon': 86231, 'transcript': 237079}
== featuretype_strand_counts ==
DB1: {('CDS', '+'): 421443, ('CDS', '-'): 403585, ('Selenocysteine', '+'): 49, ('Selenocysteine', '-'): 70, ('UTR', '+'): 177274, ('UTR', '-'): 172128, ('exon', '+'): 769867, ('exon', '-'): 729145, ('gene', '+'): 30769, ('gene', '-'): 29880, ('start_codon', '+'): 47326, ('start_codon', '-'): 45671, ('stop_codon', '+'): 44198, ('stop_codon', '-'): 42007, ('transcript', '+'): 121919, ('transcript', '-'): 115093}
DB2: {('CDS', '+'): 421476, ('CDS', '-'): 403756, ('Selenocysteine', '+'): 49, ('Selenocysteine', '-'): 70, ('UTR', '+'): 177288, ('UTR', '-'): 172189, ('exon', '+'): 769921, ('exon', '-'): 729346, ('gene', '+'): 30801, ('gene', '-'): 29907, ('start_codon', '+'): 47335, ('start_codon', '-'): 45693, ('stop_codon', '+'): 44205, ('stop_codon', '-'): 42026, ('transcript', '+'): 121952, ('transcript', '-'): 115127}
"""
