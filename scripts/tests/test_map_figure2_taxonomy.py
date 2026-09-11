import importlib.util
from pathlib import Path
import unittest


SPEC = importlib.util.spec_from_file_location("mapper", Path(__file__).resolve().parents[1] / "map_figure2_taxonomy.py")
mapper = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mapper)
PAYLOAD = b"""<TaxaSet><Taxon><TaxId>20</TaxId><ScientificName>Example strain</ScientificName>
<Rank>strain</Rank><ParentTaxId>10</ParentTaxId><AkaTaxIds><TaxId>19</TaxId></AkaTaxIds>
<LineageEx><Taxon><TaxId>2</TaxId><ScientificName>Bacteria</ScientificName><Rank>domain</Rank></Taxon>
<Taxon><TaxId>10</TaxId><ScientificName>Example species</ScientificName><Rank>species</Rank></Taxon>
</LineageEx></Taxon></TaxaSet>"""


class TaxonomyTests(unittest.TestCase):
    def test_rank_and_species_ancestor_are_separate(self):
        record = mapper.parse_taxonomy(PAYLOAD, ["20"])["20"]
        self.assertEqual((record["scientific_name"], record["taxonomic_rank"], record["domain"]),
                         ("Example strain", "strain", "Bacteria"))
        self.assertEqual(record["species_taxid"], "10")

    def test_merged_id_preserves_original_group(self):
        records = mapper.parse_taxonomy(PAYLOAD, ["19", "20"])
        self.assertEqual(records["19"]["mapping_status"], "merged")
        rows = [{"taxon_group": "taxid:19", "barcode_accession_pairs": "3"},
                {"taxon_group": "taxid:20", "barcode_accession_pairs": "7"}]
        result = mapper.annotate(rows, records, {"retrieved_utc": "test", "response_sha256": "test"})
        self.assertEqual(len(result), 2)
        self.assertEqual([{key: row[key] for key in rows[0]} for row in result], rows)

    def test_missing_taxid_fails(self):
        with self.assertRaisesRegex(ValueError, "Unresolved"):
            mapper.parse_taxonomy(PAYLOAD, ["999"])

    def test_nested_lineage_is_not_a_requested_record(self):
        with self.assertRaisesRegex(ValueError, "Unresolved"):
            mapper.parse_taxonomy(PAYLOAD, ["10"])

    def test_duplicate_and_error_responses_fail(self):
        taxon = PAYLOAD.decode().removeprefix("<TaxaSet>").removesuffix("</TaxaSet>")
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            mapper.parse_taxonomy(("<TaxaSet>" + taxon + taxon + "</TaxaSet>").encode(), ["20"])
        with self.assertRaisesRegex(ValueError, "invalid"):
            mapper.parse_taxonomy(b"<TaxaSet><ERROR>invalid</ERROR></TaxaSet>", ["20"])


if __name__ == "__main__":
    unittest.main()