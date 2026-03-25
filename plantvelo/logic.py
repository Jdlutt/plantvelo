"""plantvelo.logic
~~~~~~~~~~~~~~~~
Read-classification logic for PlantVelocity.

Three-state splicing model
--------------------------
Standard velocyto collapses two biologically distinct read populations into
a single "unspliced" layer:

  (A) reads whose segment *spans* an exon-intron boundary  → nascent pre-mRNA
  (B) reads that fall *entirely within* an intron, with no boundary crossing
      → intron-retained mRNA  (IR, plant-prevalent ~15-20%)

PlantVelocity separates (B) into a dedicated ``intron_retained`` layer so
that downstream kinetic modelling can treat U→IR→S explicitly.

Classification table (PlantPermissive10X)
-----------------------------------------
Condition                                         | Layer
--------------------------------------------------|------------------
Only exon-compatible models, no introns, no mix   | spliced
All models have exon–intron spanning reads        | unspliced
Only intron models, *validated*, NO edge-spanning | intron_retained  ← NEW
Only intron models (not validated), NO edge span  | intron_retained  ← NEW
Only intron models, IS edge-spanning              | unspliced  (pre-mRNA)
Ambiguity between intron-only & exon-only models  | ambiguous
Mixed models                                      | ambiguous / unspliced

IR flanking threshold
---------------------
``ir_flanking`` (default 5 bp) is the minimum distance a read segment must
lie from either end of an intron to be counted as fully-internal (IR).
Reads closer than ``ir_flanking`` to an intron edge could be mis-mapped
junction reads and are kept as ``unspliced``.

This mirrors the existing ``MIN_FLANK`` concept used in ``feature.intersects``
but is applied independently so users can tune it without touching the
exon-intron edge-spanning logic.
"""
from __future__ import annotations

from typing import Dict, List, Union
import abc
import numpy as np

# ---------------------------------------------------------------------------
# Re-export the velocyto Logic base so that plantvelo.logic is self-contained
# when used as the single import.
# ---------------------------------------------------------------------------
import velocyto as vcy
from velocyto.logic import Logic  # noqa: F401 – re-export for convenience


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _segment_is_fully_internal(segment, feature, ir_flanking: int) -> bool:
    """Return True when *segment* lies strictly inside *feature* by at least
    *ir_flanking* bases on each side.

    This excludes reads that just clip an intron boundary and would otherwise
    be mis-classified as IR.

    Parameters
    ----------
    segment : tuple[int, int]
        (start, end) of the aligned read segment (0-based, half-open or
        closed – we use the same convention as velocyto Feature coordinates).
    feature : vcy.Feature
        The intron feature the segment was matched to.
    ir_flanking : int
        Minimum bases of buffer required on each side of the segment inside
        the intron.
    """
    seg_start, seg_end = segment[0], segment[-1]
    return (seg_start >= feature.start + ir_flanking) and \
           (seg_end <= feature.end - ir_flanking)


# ---------------------------------------------------------------------------
# PlantPermissive10X  –  the primary logic for PlantVelocity
# ---------------------------------------------------------------------------

class PlantPermissive10X(Logic):
    """Permissive 10X logic extended with an *intron_retained* layer.

    Behaviour differences from ``Permissive10X``:
    - Reads that map only to intron features **and** do not span an
      exon-intron boundary **and** whose aligned segment lies fully inside
      the intron (≥ ``ir_flanking`` bases from each edge) are counted in
      the new **intron_retained** layer instead of unspliced.
    - All other classification rules are identical to ``Permissive10X``.

    Parameters
    ----------
    ir_flanking : int
        Minimum internal buffer (bp) from intron edge required to classify a
        read as IR rather than unspliced.  Default: 5.
    """

    def __init__(self, ir_flanking: int = 5) -> None:
        self.name = "PlantPermissive10X"
        self.ir_flanking = ir_flanking

    @property
    def layers(self) -> List[str]:
        return ["spliced", "unspliced", "ambiguous", "intron_retained"]

    @property
    def stranded(self) -> bool:
        return True

    @property
    def perform_validation_markup(self) -> bool:
        return True

    @property
    def accept_discordant(self) -> bool:
        return False

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _has_ir_signal(self, segments_list, ir_flanking: int) -> bool:
        """Return True if *all* segments in this transcript model mapping
        are:
          1. mapped to intron features, AND
          2. lie fully internal (≥ ir_flanking from each edge), AND
          3. do not show any exon-intron spanning.

        Called only when ``has_onlyintron_model`` is already confirmed.
        """
        for sm in segments_list:
            if not sm.maps_to_intron:
                return False
            if not _segment_is_fully_internal(sm.segment, sm.feature, ir_flanking):
                return False
        return True

    # ------------------------------------------------------------------
    # count()  –  mirrors Permissive10X.count() with IR branch inserted
    # ------------------------------------------------------------------

    def count(
        self,
        molitem: vcy.Molitem,
        cell_bcidx: int,
        dict_layers_columns: Dict[str, np.ndarray],
        geneid2ix: Dict[str, int],
    ) -> int:
        spliced = dict_layers_columns["spliced"]
        unspliced = dict_layers_columns["unspliced"]
        ambiguous = dict_layers_columns["ambiguous"]
        intron_retained = dict_layers_columns["intron_retained"]

        # No compatible transcript model
        if len(molitem.mappings_record) == 0:
            return 2

        # Compatible with ≥1 transcript model
        if len(set(i.geneid for i in molitem.mappings_record.keys())) != 1:
            # Multi-gene ambiguity
            return 3

        gene_check = set()
        has_onlyintron_model = 0
        has_only_span_exin_model = 1
        has_onlyintron_and_valid_model = 0
        has_valid_mixed_model = 0
        has_invalid_mixed_model = 0
        has_onlyexo_model = 0
        has_mixed_model = 0
        multi_gene = 0

        # Track whether *all* intron-only models also satisfy the IR signal
        # (fully internal, no boundary spanning).
        all_intrononly_are_ir = True

        for transcript_model, segments_list in molitem.mappings_record.items():
            gene_check.add(transcript_model.geneid)
            if len(gene_check) > 1:
                multi_gene = 1

            has_introns = 0
            has_exons = 0
            has_exseg_with_spliced_flag = 0
            has_validated_intron = 0
            has_exin_intron_span = 0
            has_non3prime = 0

            for segment_match in segments_list:
                if segment_match.maps_to_intron:
                    has_introns = 1
                    if segment_match.feature.is_validated:
                        has_validated_intron = 1
                        if segment_match.feature.end_overlaps_with_part_of(segment_match.segment):
                            downstream_exon = segment_match.feature.get_downstream_exon()
                            if downstream_exon.start_overlaps_with_part_of(segment_match.segment):
                                has_exin_intron_span = 1
                        if segment_match.feature.start_overlaps_with_part_of(segment_match.segment):
                            upstream_exon = segment_match.feature.get_upstream_exon()
                            if upstream_exon.end_overlaps_with_part_of(segment_match.segment):
                                has_exin_intron_span = 1
                elif segment_match.maps_to_exon:
                    has_exons = 1
                    if not segment_match.feature.is_last_3prime:
                        has_non3prime = 1
                    if segment_match.is_spliced:
                        has_exseg_with_spliced_flag = 1

            if has_validated_intron and not has_exons:
                has_onlyintron_and_valid_model = 1
            if has_introns and not has_exons:
                has_onlyintron_model = 1
                # Check IR criterion for this transcript model
                if has_exin_intron_span or not self._has_ir_signal(segments_list, self.ir_flanking):
                    all_intrononly_are_ir = False
            if has_exons and not has_introns:
                has_onlyexo_model = 1
                all_intrononly_are_ir = False  # mixed situation, not pure IR
            if has_exons and has_introns and not has_validated_intron and not has_exin_intron_span:
                has_invalid_mixed_model = 1
                has_mixed_model = 1
            if has_exons and has_introns and has_validated_intron and not has_exin_intron_span:
                has_valid_mixed_model = 1
                has_mixed_model = 1
            if not has_exin_intron_span:
                has_only_span_exin_model = 0

        if multi_gene:
            return 1

        if not len(molitem.mappings_record):
            return 2

        gene_ix = geneid2ix[transcript_model.geneid]

        # ----------------------------------------------------------------
        # Classification tree
        # ----------------------------------------------------------------

        if has_onlyexo_model and not has_onlyintron_model and not has_mixed_model:
            # Normal exonic read → spliced
            spliced[gene_ix, cell_bcidx] += 1
            return 0

        if has_only_span_exin_model:
            # All models have exon–intron spanning reads → nascent pre-mRNA
            unspliced[gene_ix, cell_bcidx] += 1
            return 0

        # ----------------------------------------------------------------
        # Intron-only models: decide IR vs unspliced
        # ----------------------------------------------------------------
        if has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
            if all_intrononly_are_ir:
                # Validated intron, fully internal → intron_retained
                intron_retained[gene_ix, cell_bcidx] += 1
            else:
                # Validated intron but spans edge (or too close) → unspliced
                unspliced[gene_ix, cell_bcidx] += 1
            return 0

        if has_onlyintron_model and not has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
            if all_intrononly_are_ir:
                # Non-validated intron, fully internal → intron_retained
                intron_retained[gene_ix, cell_bcidx] += 1
            else:
                # Non-validated, edge-spanning → unspliced (same as Permissive10X)
                unspliced[gene_ix, cell_bcidx] += 1
            return 0

        # ----------------------------------------------------------------
        # Mixed / ambiguous cases (unchanged from Permissive10X)
        # ----------------------------------------------------------------
        if has_invalid_mixed_model and not has_valid_mixed_model and not has_onlyintron_model and not has_onlyexo_model and not has_only_span_exin_model:
            unspliced[gene_ix, cell_bcidx] += 1
            return 0

        if has_valid_mixed_model and not has_onlyintron_model and not has_onlyexo_model and not has_only_span_exin_model:
            unspliced[gene_ix, cell_bcidx] += 1
            return 0

        if has_onlyintron_model and has_onlyexo_model and not has_mixed_model:
            ambiguous[gene_ix, cell_bcidx] += 1
            return 0

        if has_onlyintron_model and not has_onlyexo_model and has_mixed_model:
            unspliced[gene_ix, cell_bcidx] += 1
            return 0

        if not has_onlyintron_model and has_onlyexo_model and has_mixed_model:
            ambiguous[gene_ix, cell_bcidx] += 1
            return 0

        if has_onlyintron_model and has_onlyexo_model and has_mixed_model:
            ambiguous[gene_ix, cell_bcidx] += 1
            return 0

        return 4


# ---------------------------------------------------------------------------
# PlantValidated10X  –  stricter variant: only validated introns count as IR
# ---------------------------------------------------------------------------

class PlantValidated10X(Logic):
    """Like ``PlantPermissive10X`` but only *validated* introns (confirmed by
    spanning reads in the first BAM pass) are counted as ``intron_retained``.
    Non-validated intron-only reads are discarded (same policy as
    ``ValidatedIntrons10X``).

    Parameters
    ----------
    ir_flanking : int
        Same as ``PlantPermissive10X``.  Default: 5.
    """

    def __init__(self, ir_flanking: int = 5) -> None:
        self.name = "PlantValidated10X"
        self.ir_flanking = ir_flanking

    @property
    def layers(self) -> List[str]:
        return ["spliced", "unspliced", "ambiguous", "intron_retained"]

    @property
    def stranded(self) -> bool:
        return True

    @property
    def perform_validation_markup(self) -> bool:
        return True

    @property
    def accept_discordant(self) -> bool:
        return False

    def count(
        self,
        molitem: vcy.Molitem,
        cell_bcidx: int,
        dict_layers_columns: Dict[str, np.ndarray],
        geneid2ix: Dict[str, int],
    ) -> None:
        spliced = dict_layers_columns["spliced"]
        unspliced = dict_layers_columns["unspliced"]
        ambiguous = dict_layers_columns["ambiguous"]
        intron_retained = dict_layers_columns["intron_retained"]

        if len(molitem.mappings_record) == 0:
            return

        if len(set(i.geneid for i in molitem.mappings_record.keys())) != 1:
            return

        gene_check = set()
        has_onlyintron_model = 0
        has_only_span_exin_model = 1
        has_onlyintron_and_valid_model = 0
        has_valid_mixed_model = 0
        has_invalid_mixed_model = 0
        has_onlyexo_model = 0
        has_mixed_model = 0
        multi_gene = 0
        all_intrononly_are_ir = True

        for transcript_model, segments_list in molitem.mappings_record.items():
            gene_check.add(transcript_model.geneid)
            if len(gene_check) > 1:
                multi_gene = 1

            has_introns = 0
            has_exons = 0
            has_exseg_with_spliced_flag = 0
            has_validated_intron = 0
            has_exin_intron_span = 0
            has_non3prime = 0

            for segment_match in segments_list:
                if segment_match.maps_to_intron:
                    has_introns = 1
                    if segment_match.feature.is_validated:
                        has_validated_intron = 1
                        if segment_match.feature.end_overlaps_with_part_of(segment_match.segment):
                            downstream_exon = segment_match.feature.get_downstream_exon()
                            if downstream_exon.start_overlaps_with_part_of(segment_match.segment):
                                has_exin_intron_span = 1
                        if segment_match.feature.start_overlaps_with_part_of(segment_match.segment):
                            upstream_exon = segment_match.feature.get_upstream_exon()
                            if upstream_exon.end_overlaps_with_part_of(segment_match.segment):
                                has_exin_intron_span = 1
                elif segment_match.maps_to_exon:
                    has_exons = 1
                    if not segment_match.feature.is_last_3prime:
                        has_non3prime = 1
                    if segment_match.is_spliced:
                        has_exseg_with_spliced_flag = 1

            if has_validated_intron and not has_exons:
                has_onlyintron_and_valid_model = 1
            if has_introns and not has_exons:
                has_onlyintron_model = 1
                if has_exin_intron_span or not all(
                    _segment_is_fully_internal(sm.segment, sm.feature, self.ir_flanking)
                    for sm in segments_list if sm.maps_to_intron
                ):
                    all_intrononly_are_ir = False
            if has_exons and not has_introns:
                has_onlyexo_model = 1
                all_intrononly_are_ir = False
            if has_exons and has_introns and not has_validated_intron and not has_exin_intron_span:
                has_invalid_mixed_model = 1
                has_mixed_model = 1
            if has_exons and has_introns and has_validated_intron and not has_exin_intron_span:
                has_valid_mixed_model = 1
                has_mixed_model = 1
            if not has_exin_intron_span:
                has_only_span_exin_model = 0

        if multi_gene:
            return

        if not len(molitem.mappings_record):
            return

        gene_ix = geneid2ix[transcript_model.geneid]

        if has_onlyexo_model and not has_onlyintron_model and not has_mixed_model:
            spliced[gene_ix, cell_bcidx] += 1
            return

        if has_only_span_exin_model:
            unspliced[gene_ix, cell_bcidx] += 1
            return

        if has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
            if all_intrononly_are_ir:
                intron_retained[gene_ix, cell_bcidx] += 1
            else:
                unspliced[gene_ix, cell_bcidx] += 1
            return

        # Non-validated intron-only → discard (stricter than Permissive)
        if has_onlyintron_model and not has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
            return

        if has_invalid_mixed_model and not has_valid_mixed_model and not has_onlyintron_model and not has_onlyexo_model and not has_only_span_exin_model:
            return

        if has_valid_mixed_model and not has_onlyintron_model and not has_onlyexo_model and not has_only_span_exin_model:
            unspliced[gene_ix, cell_bcidx] += 1
            return

        if has_onlyintron_model and has_onlyexo_model and not has_mixed_model:
            ambiguous[gene_ix, cell_bcidx] += 1
            return

        if has_onlyintron_model and not has_onlyexo_model and has_mixed_model:
            ambiguous[gene_ix, cell_bcidx] += 1
            return

        if not has_onlyintron_model and has_onlyexo_model and has_mixed_model:
            ambiguous[gene_ix, cell_bcidx] += 1
            return

        if has_onlyintron_model and has_onlyexo_model and has_mixed_model:
            ambiguous[gene_ix, cell_bcidx] += 1
            return


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

#: All logic classes provided by plantvelo, usable via --logic on the CLI
PLANT_LOGICS = {
    "PlantPermissive10X": PlantPermissive10X,
    "PlantValidated10X": PlantValidated10X,
}
