import requests
from typing import Dict, Optional
from scripts.common.models import Variant

GNOMAD_API = "https://gnomad.broadinstitute.org/api"

VARIANT_QUERY = """
query GnomadVariant($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    genome {
      af
      populations { id af }
    }
  }
}
"""

def _graphql_query(query: str, variables: dict) -> Optional[dict]:
    try:
        resp = requests.post(GNOMAD_API, json={"query": query, "variables": variables}, timeout=30)
        resp.raise_for_status()
        return resp.json()
    except Exception:
        return None

def query_gnomad(variant: Variant) -> Dict:
    chrom_num = variant.chrom.replace("chr", "")
    variant_id = f"{chrom_num}-{variant.pos}-{variant.ref}-{variant.alt}"

    data = _graphql_query(VARIANT_QUERY, {"variantId": variant_id, "dataset": "gnomad_r4"})

    if not data or not data.get("data", {}).get("variant"):
        return {"gnomad_all": None, "gnomad_eas": None}

    genome = data["data"]["variant"].get("genome") or {}
    gnomad_all = genome.get("af")
    gnomad_eas = None
    for pop in genome.get("populations", []):
        if pop["id"] == "eas":
            gnomad_eas = pop["af"]
            break

    return {"gnomad_all": gnomad_all, "gnomad_eas": gnomad_eas}


if __name__ == "__main__":
    import sys, json
    if len(sys.argv) < 2:
        print(json.dumps({"error": "Usage: python -m scripts.korean_pop.query_gnomad 'chr17:7577120 G>A'"}))
        sys.exit(1)
    v = Variant.from_string(sys.argv[1])
    print(json.dumps(query_gnomad(v), indent=2))
