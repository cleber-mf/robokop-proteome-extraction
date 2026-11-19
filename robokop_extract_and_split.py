#!/usr/bin/env python3
"""Extract ROBOKOP protein-centred assertions and make leakage-safe splits.

Two scopes use the same pipeline:
  --all-proteins            discover every protein-like node in ROBOKOP KG
  --protein-ids FILE        use a future user-supplied set of protein CURIEs

The state database makes the long all-protein job restartable.  The final
deliverables retain the original structure:
  interactions_raw.csv, examples.csv, train.csv, heldout.csv, audit.json.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import sqlite3
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Iterable, Iterator


ENDPOINT = "https://robokop-automat.apps.renci.org/robokopkg/cypher"
METADATA_ENDPOINT = "https://robokop-automat.apps.renci.org/robokopkg/metadata"
DEFAULT_PROTEIN_NODE_CATEGORIES = {"Gene", "Protein", "GeneOrProtein"}

# Identity and hierarchy links are useful for normalization, but not normally
# training targets.  The namespace-free terminal value is compared.
STRUCTURAL_PREDICATES = {
    "subclass_of", "superclass_of", "same_as", "exact_match", "close_match",
    "broad_match", "narrow_match", "equivalent_to", "maps_to",
}
LITERATURE_ONLY_PREDICATES = {"occurs_together_in_literature_with"}
SYMMETRIC_PREDICATES = {
    "related_to", "interacts_with", "physically_interacts_with",
    "molecularly_interacts_with", "associated_with", "correlated_with",
    "coexpressed_with", "coexists_with", "overlaps", "similar_to", "same_as",
    "exact_match", "close_match", "occurs_together_in_literature_with",
}

RAW_FIELDS = [
    "protein_id", "protein_name", "protein_categories", "partner_id", "partner_name",
    "partner_categories", "predicate", "subject_id", "object_id", "edge_id",
    "edge_properties", "example_id", "split",
]
EXAMPLE_FIELDS = [
    "example_id", "subject_id", "predicate", "object_id", "anchor_proteins",
    "partner_ids", "partner_names", "partner_categories", "evidence_edge_count",
    "edge_ids", "split",
]


def predicate_leaf(value: str) -> str:
    return value.rsplit(":", 1)[-1]


def json_cell(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def categories(value: Any) -> list[str]:
    """Normalize Neo4j labels or CSV-serialized labels to category names."""
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except json.JSONDecodeError:
            value = [value]
    if not isinstance(value, list):
        return []
    return sorted({str(label).rsplit(":", 1)[-1] for label in value})


def cypher_string(value: str) -> str:
    return json.dumps(value)


def post_cypher(query: str, timeout: int) -> list[dict[str, Any]]:
    request = urllib.request.Request(
        ENDPOINT,
        data=json.dumps({"query": query}).encode("utf-8"),
        headers={"Accept": "application/json", "Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            body = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as error:
        detail = error.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"ROBOKOP HTTP {error.code}: {detail[:1000]}") from error
    except urllib.error.URLError as error:
        raise RuntimeError(f"Unable to reach ROBOKOP: {error}") from error

    if body.get("errors"):
        raise RuntimeError(f"ROBOKOP returned Cypher errors: {json.dumps(body['errors'])}")
    results = body.get("results", [])
    if not results:
        raise RuntimeError(f"Unexpected ROBOKOP response: {json.dumps(body)[:1000]}")
    columns = results[0]["columns"]
    return [dict(zip(columns, item["row"], strict=True)) for item in results[0].get("data", [])]


def fetch_metadata(timeout: int) -> dict[str, Any]:
    request = urllib.request.Request(METADATA_ENDPOINT, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = json.loads(response.read().decode("utf-8"))
        return payload if isinstance(payload, dict) else {"response": payload}
    except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError) as error:
        return {"metadata_retrieval_error": str(error)}


def make_edge_query(ids: list[str]) -> str:
    literals = ", ".join(cypher_string(value) for value in ids)
    return f"""
UNWIND [{literals}] AS requested_protein_id
MATCH (protein {{id: requested_protein_id}})-[edge]-(partner)
RETURN
  protein.id AS protein_id,
  protein.name AS protein_name,
  labels(protein) AS protein_categories,
  partner.id AS partner_id,
  partner.name AS partner_name,
  labels(partner) AS partner_categories,
  type(edge) AS predicate,
  startNode(edge).id AS subject_id,
  endNode(edge).id AS object_id,
  elementId(edge) AS edge_id,
  properties(edge) AS edge_properties
""".strip()


def make_protein_page_query(protein_categories: set[str], offset: int, limit: int) -> str:
    labels = ", ".join(cypher_string(f"biolink:{category}") for category in sorted(protein_categories))
    return f"""
MATCH (protein)
WHERE any(label IN labels(protein) WHERE label IN [{labels}])
  AND protein.id IS NOT NULL
RETURN protein.id AS protein_id
ORDER BY protein.id
SKIP {offset}
LIMIT {limit}
""".strip()


def chunks(values: list[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def read_protein_ids(path: Path) -> list[str]:
    ids = [
        line.strip() for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    return list(dict.fromkeys(ids))


def selected(row: dict[str, Any], args: argparse.Namespace) -> bool:
    predicate = predicate_leaf(str(row["predicate"]))
    if not args.include_ontology and predicate in STRUCTURAL_PREDICATES:
        return False
    if not args.include_literature_cooccurrence and predicate in LITERATURE_ONLY_PREDICATES:
        return False
    return not args.partner_categories or bool(
        set(categories(row["partner_categories"])).intersection(args.partner_categories)
    )


def example_id(row: dict[str, Any]) -> str:
    predicate = str(row["predicate"])
    subject, object_ = str(row["subject_id"]), str(row["object_id"])
    if predicate_leaf(predicate) in SYMMETRIC_PREDICATES:
        subject, object_ = sorted((subject, object_))
    return "|".join((subject, predicate, object_))


def stable_order_key(seed: int, value: str) -> str:
    return hashlib.sha256(f"{seed}|{value}".encode("utf-8")).hexdigest()


class StateStore:
    """SQLite-backed extraction state so full-KG jobs can resume safely."""

    def __init__(self, path: Path, signature: str, reset: bool) -> None:
        if reset:
            for candidate in (path, Path(f"{path}-wal"), Path(f"{path}-shm")):
                if candidate.exists():
                    candidate.unlink()
        path.parent.mkdir(parents=True, exist_ok=True)
        self.connection = sqlite3.connect(path)
        self.connection.row_factory = sqlite3.Row
        self.connection.execute("PRAGMA journal_mode=WAL")
        self.connection.execute("PRAGMA synchronous=NORMAL")
        self.connection.execute("CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
        saved = self.get_metadata("signature")
        if saved is not None and saved != signature:
            raise ValueError(
                "Existing state belongs to a different extraction configuration. "
                "Use a different --state-db or intentionally pass --reset-state."
            )
        self.set_metadata("signature", signature)
        self.connection.executescript("""
            CREATE TABLE IF NOT EXISTS processed_proteins (protein_id TEXT PRIMARY KEY);
            CREATE TABLE IF NOT EXISTS raw (
                row_id INTEGER PRIMARY KEY,
                protein_id TEXT, protein_name TEXT, protein_categories TEXT,
                partner_id TEXT, partner_name TEXT, partner_categories TEXT,
                predicate TEXT, subject_id TEXT, object_id TEXT, edge_id TEXT,
                edge_properties TEXT, example_id TEXT, selected INTEGER NOT NULL
            );
            CREATE INDEX IF NOT EXISTS raw_selected_example ON raw(selected, example_id);
            CREATE TABLE IF NOT EXISTS examples (
                example_id TEXT PRIMARY KEY,
                subject_id TEXT, predicate TEXT, object_id TEXT,
                anchor_proteins TEXT, partner_ids TEXT, partner_names TEXT,
                partner_categories TEXT, evidence_edge_count INTEGER, edge_ids TEXT,
                sort_key TEXT, split TEXT
            );
            CREATE TABLE IF NOT EXISTS example_anchors (
                example_id TEXT, protein_id TEXT,
                PRIMARY KEY (example_id, protein_id)
            );
        """)
        self.connection.commit()

    def close(self) -> None:
        self.connection.close()

    def get_metadata(self, key: str) -> str | None:
        row = self.connection.execute("SELECT value FROM metadata WHERE key = ?", (key,)).fetchone()
        return None if row is None else str(row[0])

    def set_metadata(self, key: str, value: str) -> None:
        self.connection.execute(
            "INSERT INTO metadata(key, value) VALUES (?, ?) "
            "ON CONFLICT(key) DO UPDATE SET value = excluded.value",
            (key, value),
        )

    def unprocessed(self, ids: list[str]) -> list[str]:
        if not ids:
            return []
        placeholders = ",".join("?" for _ in ids)
        seen = {
            row[0] for row in self.connection.execute(
                f"SELECT protein_id FROM processed_proteins WHERE protein_id IN ({placeholders})", ids
            )
        }
        return [value for value in ids if value not in seen]

    def insert_batch(self, rows: list[dict[str, Any]], processed_ids: list[str], args: argparse.Namespace) -> None:
        with self.connection:
            values = []
            for row in rows:
                values.append((
                    str(row["protein_id"]), row["protein_name"], json_cell(row["protein_categories"]),
                    str(row["partner_id"]), row["partner_name"], json_cell(row["partner_categories"]),
                    str(row["predicate"]), str(row["subject_id"]), str(row["object_id"]),
                    str(row["edge_id"]), json_cell(row["edge_properties"]), example_id(row),
                    int(selected(row, args)),
                ))
            self.connection.executemany(
                """INSERT INTO raw(
                    protein_id, protein_name, protein_categories, partner_id, partner_name,
                    partner_categories, predicate, subject_id, object_id, edge_id,
                    edge_properties, example_id, selected
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                values,
            )
            self.connection.executemany(
                "INSERT OR IGNORE INTO processed_proteins(protein_id) VALUES (?)",
                ((value,) for value in processed_ids),
            )


def iter_target_batches(store: StateStore, args: argparse.Namespace) -> Iterator[list[str]]:
    if args.protein_ids is not None:
        ids = read_protein_ids(args.protein_ids)
        if not ids:
            raise ValueError("--protein-ids contains no CURIEs.")
        for batch in chunks(ids, args.batch_size):
            pending = store.unprocessed(batch)
            if pending:
                yield pending
        return

    offset = 0
    page_number = 0
    while True:
        page_number += 1
        query = make_protein_page_query(args.protein_node_categories, offset, args.protein_page_size)
        page = post_cypher(query, args.timeout)
        ids = [str(row["protein_id"]) for row in page]
        if not ids:
            return
        for batch in chunks(ids, args.batch_size):
            pending = store.unprocessed(batch)
            if pending:
                yield pending
        offset += len(ids)
        print(f"Discovered protein page {page_number}: {len(ids)} nodes (offset {offset})", file=sys.stderr)
        if len(ids) < args.protein_page_size:
            return


def make_example(group: list[dict[str, Any]], seed: int) -> dict[str, Any]:
    first = group[0]
    partner_categories = sorted({item for row in group for item in categories(row["partner_categories"])})
    identifier = first["example_id"]
    return {
        "example_id": identifier,
        "subject_id": first["subject_id"],
        "predicate": first["predicate"],
        "object_id": first["object_id"],
        "anchor_proteins": ";".join(sorted({row["protein_id"] for row in group})),
        "partner_ids": ";".join(sorted({row["partner_id"] for row in group})),
        "partner_names": ";".join(sorted({str(row["partner_name"]) for row in group if row["partner_name"]})),
        "partner_categories": ";".join(partner_categories),
        # In --all-proteins mode a protein--protein edge is encountered from
        # both protein anchors.  Count underlying KG edges, not those repeated
        # anchor encounters.
        "evidence_edge_count": len({row["edge_id"] for row in group}),
        "edge_ids": ";".join(sorted({row["edge_id"] for row in group})),
        "sort_key": stable_order_key(seed, identifier),
        "split": "train",
    }


def rebuild_examples(store: StateStore, seed: int) -> None:
    connection = store.connection
    with connection:
        connection.execute("DELETE FROM examples")
        connection.execute("DELETE FROM example_anchors")
        cursor = connection.execute(
            """SELECT protein_id, protein_name, protein_categories, partner_id, partner_name,
                      partner_categories, predicate, subject_id, object_id, edge_id,
                      edge_properties, example_id
               FROM raw WHERE selected = 1 ORDER BY example_id"""
        )
        current_id: str | None = None
        group: list[dict[str, Any]] = []

        def persist(rows: list[dict[str, Any]]) -> None:
            if not rows:
                return
            item = make_example(rows, seed)
            connection.execute(
                """INSERT INTO examples VALUES
                   (:example_id, :subject_id, :predicate, :object_id, :anchor_proteins,
                    :partner_ids, :partner_names, :partner_categories, :evidence_edge_count,
                    :edge_ids, :sort_key, :split)""",
                item,
            )
            connection.executemany(
                "INSERT INTO example_anchors(example_id, protein_id) VALUES (?, ?)",
                ((item["example_id"], protein) for protein in item["anchor_proteins"].split(";")),
            )

        for db_row in cursor:
            row = dict(db_row)
            if current_id is not None and row["example_id"] != current_id:
                persist(group)
                group = []
            current_id = row["example_id"]
            group.append(row)
        persist(group)


def apply_edge_split(store: StateStore, holdout_fraction: float) -> None:
    connection = store.connection
    with connection:
        connection.execute("UPDATE examples SET split = 'train'")
        strata = connection.execute(
            "SELECT predicate, partner_categories, COUNT(*) AS count FROM examples GROUP BY predicate, partner_categories"
        ).fetchall()
        for stratum in strata:
            count = int(stratum["count"])
            if count < 2:
                continue
            n_test = min(max(1, round(count * holdout_fraction)), count - 1)
            chosen = connection.execute(
                """SELECT example_id FROM examples
                   WHERE predicate = ? AND partner_categories = ?
                   ORDER BY sort_key LIMIT ?""",
                (stratum["predicate"], stratum["partner_categories"], n_test),
            ).fetchall()
            connection.executemany(
                "UPDATE examples SET split = 'heldout' WHERE example_id = ?",
                ((row["example_id"],) for row in chosen),
            )


def apply_protein_split(store: StateStore, holdout_fraction: float, seed: int) -> None:
    connection = store.connection
    anchors = [row[0] for row in connection.execute("SELECT DISTINCT protein_id FROM example_anchors")]
    if len(anchors) < 2:
        raise ValueError("--split-mode protein needs at least two proteins with selected assertions.")
    anchors.sort(key=lambda value: stable_order_key(seed, value))
    n_test = min(max(1, round(len(anchors) * holdout_fraction)), len(anchors) - 1)
    heldout_anchors = anchors[:n_test]
    with connection:
        connection.execute("UPDATE examples SET split = 'train'")
        for batch in chunks(heldout_anchors, 500):
            placeholders = ",".join("?" for _ in batch)
            connection.execute(
                f"""UPDATE examples SET split = 'heldout'
                    WHERE example_id IN (
                      SELECT example_id FROM example_anchors WHERE protein_id IN ({placeholders})
                    )""",
                batch,
            )


def write_csv(path: Path, fields: list[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def export_deliverables(store: StateStore, output_dir: Path, audit: dict[str, Any]) -> None:
    connection = store.connection
    raw_cursor = connection.execute(
        """SELECT r.protein_id, r.protein_name, r.protein_categories, r.partner_id,
                  r.partner_name, r.partner_categories, r.predicate, r.subject_id,
                  r.object_id, r.edge_id, r.edge_properties, r.example_id,
                  COALESCE(e.split, 'excluded') AS split
           FROM raw r LEFT JOIN examples e ON r.example_id = e.example_id
           ORDER BY r.row_id"""
    )
    write_csv(output_dir / "interactions_raw.csv", RAW_FIELDS, (dict(row) for row in raw_cursor))
    examples = """SELECT example_id, subject_id, predicate, object_id, anchor_proteins,
                         partner_ids, partner_names, partner_categories, evidence_edge_count,
                         edge_ids, split FROM examples ORDER BY example_id"""
    write_csv(output_dir / "examples.csv", EXAMPLE_FIELDS, (dict(row) for row in connection.execute(examples)))
    for split, filename in (("train", "train.csv"), ("heldout", "heldout.csv")):
        cursor = connection.execute(examples.replace("ORDER BY example_id", "WHERE split = ? ORDER BY example_id"), (split,))
        write_csv(output_dir / filename, EXAMPLE_FIELDS, (dict(row) for row in cursor))
    (output_dir / "audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    scope = parser.add_mutually_exclusive_group(required=True)
    scope.add_argument("--all-proteins", action="store_true",
                       help="Discover every ROBOKOP Gene, Protein, or GeneOrProtein node.")
    scope.add_argument("--protein-ids", type=Path,
                       help="Text file of normalized protein CURIEs, one per line.")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--state-db", type=Path,
                        help="Restartable SQLite state (default: .robokop_extract_state.sqlite3 in --out-dir).")
    parser.add_argument("--reset-state", action="store_true",
                        help="Discard an existing state database for this extraction.")
    parser.add_argument("--holdout-fraction", type=float, default=0.20)
    parser.add_argument("--seed", type=int, default=20261112)
    parser.add_argument("--batch-size", type=int, default=25,
                        help="Protein CURIEs per edge query; lower it for high-degree nodes.")
    parser.add_argument("--protein-page-size", type=int, default=500,
                        help="Node-discovery page size for --all-proteins (maximum 900).")
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--split-mode", choices=("edge", "protein"), default="edge")
    parser.add_argument("--protein-node-categories", default="Gene,Protein,GeneOrProtein",
                        help="Biolink categories treated as proteins in --all-proteins mode.")
    parser.add_argument("--partner-categories", default="",
                        help="Comma-separated categories to retain; empty retains every partner type.")
    parser.add_argument("--include-ontology", action="store_true")
    parser.add_argument("--include-literature-cooccurrence", action="store_true")
    args = parser.parse_args()
    if not 0 < args.holdout_fraction < 1:
        parser.error("--holdout-fraction must be between 0 and 1.")
    if args.batch_size < 1 or not 1 <= args.protein_page_size <= 900:
        parser.error("--batch-size must be positive and --protein-page-size must be 1..900.")
    args.partner_categories = {
        item.strip().rsplit(":", 1)[-1] for item in args.partner_categories.split(",") if item.strip()
    }
    args.protein_node_categories = {
        item.strip().rsplit(":", 1)[-1] for item in args.protein_node_categories.split(",") if item.strip()
    } or DEFAULT_PROTEIN_NODE_CATEGORIES
    return args


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.state_db is None:
        args.state_db = args.out_dir / ".robokop_extract_state.sqlite3"
    input_ids = read_protein_ids(args.protein_ids) if args.protein_ids is not None else None
    signature = json_cell({
        "scope": "all-proteins" if args.all_proteins else "protein-ids",
        "input_protein_ids": input_ids,
        "protein_node_categories": sorted(args.protein_node_categories),
        "partner_categories": sorted(args.partner_categories),
        "include_ontology": args.include_ontology,
        "include_literature_cooccurrence": args.include_literature_cooccurrence,
    })
    store = StateStore(args.state_db, signature, args.reset_state)
    try:
        store.set_metadata("kg_metadata", json_cell(fetch_metadata(args.timeout)))
        completed = 0
        for batch in iter_target_batches(store, args):
            completed += 1
            print(f"Querying batch {completed} ({len(batch)} proteins)...", file=sys.stderr)
            rows = post_cypher(make_edge_query(batch), args.timeout)
            store.insert_batch(rows, batch, args)

        rebuild_examples(store, args.seed)
        if args.split_mode == "edge":
            apply_edge_split(store, args.holdout_fraction)
        else:
            apply_protein_split(store, args.holdout_fraction, args.seed)

        connection = store.connection
        counts = {"raw_edges_retrieved": connection.execute("SELECT COUNT(*) FROM raw").fetchone()[0]}
        counts["logical_examples_selected"] = connection.execute("SELECT COUNT(*) FROM examples").fetchone()[0]
        counts["train_examples"] = connection.execute("SELECT COUNT(*) FROM examples WHERE split = 'train'").fetchone()[0]
        counts["heldout_examples"] = connection.execute("SELECT COUNT(*) FROM examples WHERE split = 'heldout'").fetchone()[0]
        counts["excluded_raw_edges"] = connection.execute("SELECT COUNT(*) FROM raw WHERE selected = 0").fetchone()[0]
        predicate_counts = dict(connection.execute("SELECT predicate, COUNT(*) FROM examples GROUP BY predicate ORDER BY predicate"))
        category_counts = dict(connection.execute("SELECT partner_categories, COUNT(*) FROM examples GROUP BY partner_categories ORDER BY partner_categories"))
        audit = {
            "endpoint": ENDPOINT,
            "metadata_endpoint": METADATA_ENDPOINT,
            "kg_metadata": json.loads(store.get_metadata("kg_metadata") or "{}"),
            "retrieved_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "scope": "all-proteins" if args.all_proteins else "protein-ids",
            "input_protein_ids": input_ids,
            "options": {
                "holdout_fraction": args.holdout_fraction, "seed": args.seed,
                "split_mode": args.split_mode,
                "protein_node_categories": sorted(args.protein_node_categories),
                "partner_categories": sorted(args.partner_categories),
                "include_ontology": args.include_ontology,
                "include_literature_cooccurrence": args.include_literature_cooccurrence,
            },
            "counts": counts,
            "examples_by_predicate": predicate_counts,
            "examples_by_partner_categories": category_counts,
            "cypher_templates": {
                "protein_discovery": make_protein_page_query(args.protein_node_categories, 0, args.protein_page_size),
                "edge_extraction": make_edge_query(["<protein CURIE batch>"]),
            },
            "state_db": str(args.state_db),
            "leakage_check": "passed: logical example IDs are disjoint across train and heldout",
        }
        export_deliverables(store, args.out_dir, audit)
        print(json.dumps(counts, indent=2))
        return 0
    finally:
        store.close()


if __name__ == "__main__":
    raise SystemExit(main())
