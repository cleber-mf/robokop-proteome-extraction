# ROBOKOP protein-neighbourhood dataset

The accompanying script exports one-hop protein interactions from ROBOKOP KG,
then makes a reproducible train/held-out split.  It has two scopes with the
same output schema and filtering/splitting rules:

- `--all-proteins` discovers every ROBOKOP node labelled `Gene`, `Protein`, or
  `GeneOrProtein` and extracts its one-hop neighbourhood.
- `--protein-ids` retains the original focused mode for a future user-supplied
  list of normalized protein CURIEs.

This is deliberately a **protein-centred knowledge-graph dataset**, not a
strict physical interactome.  It can include gene/protein, chemical, biological
process/activity, pathway, disease, phenotype, anatomy, and other endpoint
types.  In a strict protein-protein interactome, filter to gene/protein partners
and then keep only interaction predicates supported by the sources you accept.

## 1. Select the scope

For the full protein set, no input file is required.  The all-protein mode is
the default workflow described here.  It is a large job: its output and its
restartable SQLite state can require substantial disk space because ROBOKOP KG
contains many millions of nodes and edges.

For a focused future run, make a UTF-8 text file containing one normalized
CURIE per line.  For example:

```text
# Examples only: replace with your proteins
HGNC:9236
UniProtKB:P37231
```

Resolve aliases first with ROBOKOP's Node Normalizer and record the final
CURIEs.  Do not use gene symbols alone: the graph identifier is what makes the
query unambiguous and reproducible.

## 2. Run extraction

Run, using the full paths appropriate for your machine (the script uses only
the Python standard library):

```bash
python3 robokop_extract_and_split.py \
  --all-proteins \
  --out-dir robokop_all_proteins
```

This retains every non-structural partner type ROBOKOP returns: gene/protein,
chemical, biological process/activity, pathway, disease, phenotype, anatomy,
and others.  Add `--partner-categories
GeneOrProtein,ChemicalEntity,BiologicalProcessOrActivity` only when a narrower
target domain is wanted.  To retain all proteins but use a narrower definition
of a protein-like source node, change
`--protein-node-categories`; its default is `Gene,Protein,GeneOrProtein`.

The state file `robokop_all_proteins/.robokop_extract_state.sqlite3` is updated
atomically after each successful batch.  If a long run stops, rerun the same
command and it will skip completed proteins.  Do not change filters or scope
when resuming; use a different output directory, or deliberately add
`--reset-state` to start again.

For a future focused query, use the same command with the scope changed:

```bash
python3 robokop_extract_and_split.py \
  --protein-ids proteins.txt \
  --out-dir robokop_selected_proteins \
  --partner-categories GeneOrProtein,ChemicalEntity,BiologicalProcessOrActivity
```

In both modes, the default excludes ontology/identity mappings such as `subclass_of` and the
literature-co-occurrence predicate; these are not biological interactions.
Use `--include-ontology` or `--include-literature-cooccurrence` only when that
extra signal is a deliberate part of the task.

The ROBOKOP endpoint used is:

```text
POST https://robokop-automat.apps.renci.org/robokopkg/cypher
```

The KG documentation identifies the current release, node/edge counts, and
direct Cypher/TRAPI/download endpoints.  Save `audit.json` beside every model
run because it retains graph metadata, query templates, retrieval date, scope,
filters, and random seed.

## 3. Pick the evaluation you mean

`--split-mode edge` (default) is a stratified held-out-edge test.  It asks:
given a protein already seen in training, can the model recover a missing
protein--predicate--partner assertion?  Evidence records with the same logical
assertion stay together, so an assertion cannot appear in training via one
source and in the held-out set via another.

`--split-mode protein` is appropriate for several proteins and asks a harder
question: can the model generalize to a protein whose entire requested
neighbourhood was not in training?  It cannot be used for one protein.

Neither random split is a prospective discovery benchmark.  For that, rebuild
the graph from source releases available before a cutoff date and test on later
assertions; remove all sources/evidence released after the cutoff from training.

## 4. Use the output correctly

| File | What it contains |
| --- | --- |
| `interactions_raw.csv` | Every returned protein-anchored graph edge, including properties/provenance. `excluded` rows show what filtering removed. A protein–protein edge can occur under both protein anchors; use `examples.csv` for de-duplicated modeling records. |
| `examples.csv` | De-duplicated logical assertions and their split assignment. |
| `train.csv` / `heldout.csv` | The model examples. Their `example_id`s are guaranteed disjoint. |
| `audit.json` | Counts, predicate/category audit, exact Cypher, split options, and leakage check. |

Before modeling, inspect `audit.json` and decide which sources and predicates
are scientifically admissible.  In particular, `related_to` is very broad;
text-mined and literature-co-occurrence edges are noisy; and GO/pathway
annotations test association prediction rather than a physical interaction.

## Direct Cypher audit

To inspect what kinds of assertions a particular protein has before exporting,
replace the CURIE below and submit it to the same endpoint:

```cypher
MATCH (p {id: "HGNC:9236"})-[e]-(x)
RETURN labels(x) AS partner_categories, type(e) AS predicate, count(*) AS edges
ORDER BY edges DESC
```

Then use the observed labels and predicates—not guessed names—to decide the
final partner and evidence filters.  For the full dataset, inspect this audit
on a small focused run first, then launch `--all-proteins` with the chosen
filters.
