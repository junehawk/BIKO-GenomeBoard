import logging
from typing import Dict, Optional

import requests

from scripts.common.cache import get_cached, set_cached
from scripts.common.config import get
from scripts.common.models import Variant

GNOMAD_API = get("api.gnomad", "https://gnomad.broadinstitute.org/api")

# Try multiple datasets - gnomAD v4 (GRCh38) first, then v2.1 (GRCh37)
VARIANT_QUERY = """
query($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id
    genome {
      ac
      an
      populations { id ac an }
    }
    exome {
      ac
      an
      populations { id ac an }
    }
  }
}
"""

logger = logging.getLogger(__name__)


def _graphql_query(query: str, variables: dict) -> Optional[dict]:
    try:
        resp = requests.post(GNOMAD_API, json={"query": query, "variables": variables}, timeout=30)
        resp.raise_for_status()
        return resp.json()
    except Exception:
        return None


def _calc_af(ac, an):
    """Calculate allele frequency from allele count and number."""
    if an and an > 0 and ac is not None:
        return ac / an
    return None


def _extract_frequencies(variant_data: dict) -> Dict:
    """Extract AF from genome or exome data, preferring genome."""
    # Try genome first, then exome
    for source_key in ["genome", "exome"]:
        source = variant_data.get(source_key)
        if not source or not source.get("an"):
            continue

        gnomad_all = _calc_af(source.get("ac"), source.get("an"))
        gnomad_eas = None
        for pop in source.get("populations", []):
            if pop["id"] == "eas" and pop.get("an"):
                gnomad_eas = _calc_af(pop.get("ac"), pop.get("an"))
                break

        return {"gnomad_all": gnomad_all, "gnomad_eas": gnomad_eas, "api_available": True}

    return {"gnomad_all": None, "gnomad_eas": None, "api_available": False}


def query_gnomad(variant: Variant) -> Dict:
    # Check cache first
    cached = get_cached(variant.chrom, variant.pos, variant.ref, variant.alt, "gnomad")
    if cached:
        return cached

    chrom_num = variant.chrom.replace("chr", "")
    variant_id = f"{chrom_num}-{variant.pos}-{variant.ref}-{variant.alt}"

    # Try gnomad_r4 (GRCh38) first, then gnomad_r2_1 (GRCh37)
    for dataset in get("api.gnomad_datasets", ["gnomad_r4", "gnomad_r2_1"]):
        data = _graphql_query(VARIANT_QUERY, {"variantId": variant_id, "dataset": dataset})

        if data and "errors" in data:
            logger.debug(f"gnomAD {dataset} error for {variant_id}: {data['errors']}")
            continue

        if not data or not data.get("data", {}).get("variant"):
            continue

        result = _extract_frequencies(data["data"]["variant"])
        if result["gnomad_all"] is not None:
            logger.info(f"gnomAD hit: {variant_id} in {dataset}, AF={result['gnomad_all']}")
            # Store in cache only if API returned data
            set_cached(variant.chrom, variant.pos, variant.ref, variant.alt, "gnomad", result)
            return result

    # No data found in any dataset
    return {"gnomad_all": None, "gnomad_eas": None, "api_available": False}


if __name__ == "__main__":
    import json
    import sys

    if len(sys.argv) < 2:
        print(json.dumps({"error": "Usage: python -m scripts.population.query_gnomad 'chr17:7577120 G>A'"}))
        sys.exit(1)
    v = Variant.from_string(sys.argv[1])
    print(json.dumps(query_gnomad(v), indent=2))
