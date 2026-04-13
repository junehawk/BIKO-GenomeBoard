"""Clinical Board runner — orchestrates the full diagnostic synthesis process."""

import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional

from scripts.clinical_board.models import AgentOpinion, BoardOpinion
from scripts.clinical_board.ollama_client import OllamaClient
from scripts.clinical_board.case_briefing import build_case_briefing
from scripts.common.config import get

logger = logging.getLogger(__name__)


def _load_agents(client: OllamaClient, model: str):
    """Lazy-load domain agents to avoid import errors if not installed."""
    from scripts.clinical_board.agents.variant_pathologist import VariantPathologist
    from scripts.clinical_board.agents.disease_geneticist import DiseaseGeneticist
    from scripts.clinical_board.agents.pgx_specialist import PGxSpecialist as PgxSpecialist
    from scripts.clinical_board.agents.literature_analyst import LiteratureAnalyst
    return [
        VariantPathologist(client=client, model=model),
        DiseaseGeneticist(client=client, model=model),
        PgxSpecialist(client=client, model=model),
        LiteratureAnalyst(client=client, model=model),
    ]


def _load_chair(client: OllamaClient, model: str):
    from scripts.clinical_board.agents.board_chair import BoardChair
    return BoardChair(client=client, model=model)


def run_clinical_board(
    report_data: dict,
    mode: str,
    ollama_url: str = None,
    agent_model: str = None,
    chair_model: str = None,
) -> Optional[BoardOpinion]:
    """Run the full Clinical Board analysis.

    Args:
        report_data: Output from run_pipeline()
        mode: "cancer" or "rare-disease"
        ollama_url: Ollama server URL (default from config)
        agent_model: Model for domain agents (default: medgemma:27b)
        chair_model: Model for Board Chair (default: gemma4:31b)

    Returns:
        BoardOpinion with synthesized diagnostic opinion, or None on failure.
    """
    start = time.time()
    agent_model = agent_model or get("clinical_board.agent_model", "medgemma:27b")
    chair_model = chair_model or get("clinical_board.chair_model", "gemma4:31b")

    # Initialize client
    client = OllamaClient(base_url=ollama_url)
    if not client.is_available():
        logger.error("[Clinical Board] Ollama server is not available. Skipping.")
        return None

    # Check models
    for model_name in [agent_model, chair_model]:
        if not client.has_model(model_name):
            logger.warning(f"[Clinical Board] Model '{model_name}' not found. "
                          f"Run 'ollama pull {model_name}' to download.")

    # Step 1: Build case briefing
    logger.info("[Clinical Board] Building case briefing...")
    briefing = build_case_briefing(report_data, mode)
    logger.debug(f"[Clinical Board] Briefing length: {len(briefing)} chars")

    # Step 2: Run domain agents in parallel
    logger.info("[Clinical Board] Running domain agents...")
    agents = _load_agents(client, agent_model)
    opinions: list[AgentOpinion] = []

    # Run agents sequentially on single GPU to avoid resource contention
    # (parallel execution causes timeouts when agents compete for GPU)
    for agent in agents:
        try:
            opinion = agent.analyze(briefing)
            opinions.append(opinion)
            logger.info(f"  [Board] {agent.agent_name}: {opinion.confidence} confidence, "
                       f"{len(opinion.findings)} findings")
        except Exception as e:
            logger.error(f"  [Board] {agent.agent_name} failed: {e}")
            opinions.append(AgentOpinion(
                agent_name=agent.agent_name, domain="error",
                findings=[{"finding": f"Agent failed: {e}", "evidence": "", "confidence": "low"}],
                confidence="low",
            ))

    # Step 3: Board Chair synthesis
    logger.info("[Clinical Board] Board Chair synthesizing opinions...")
    chair = _load_chair(client, chair_model)
    board_opinion = chair.synthesize(briefing, opinions)
    board_opinion.agent_opinions = opinions

    elapsed = time.time() - start
    logger.info(f"[Clinical Board] Complete in {elapsed:.1f}s — "
                f"Primary diagnosis: {board_opinion.primary_diagnosis} "
                f"({board_opinion.confidence} confidence)")

    return board_opinion
