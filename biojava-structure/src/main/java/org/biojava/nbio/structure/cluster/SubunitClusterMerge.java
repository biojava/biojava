package org.biojava.nbio.structure.cluster;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsembleImpl;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class SubunitClusterMerge {
    private final SubunitCluster subunitCluster;
    private static final Logger logger = LoggerFactory.getLogger(SubunitCluster.class);

    /**
     * A constructor from a single SubunitCluster. To obtain a
     * SubunitCluster with multiple Subunits, initialize different
     * SubunitClusters and merge them.
     *
     * @param subunitCluster
     *            SubunitCluster
     */
    public SubunitClusterMerge(SubunitCluster subunitCluster) {
        this.subunitCluster = subunitCluster;
    }

    /**
     * Tells whether the other SubunitCluster contains exactly the same Subunit.
     * This is checked by String equality of their residue one-letter sequences.
     *
     * @param other SubunitCluster
     * @return true if the SubunitClusters are identical, false otherwise
     */
    public boolean isIdenticalTo(SubunitCluster other) {
        String thisSequence = subunitCluster.getSubunits().get(subunitCluster.getRepresentative())
                .getProteinSequenceString();
        String otherSequence = other.getSubunits().get(other.getRepresentative())
                .getProteinSequenceString();
        return thisSequence.equals(otherSequence);
    }

    /**
     * Tells whether the other SubunitCluster contains exactly the same Subunit.
     * This is checked by equality of their entity identifiers if they are present.
     *
     * @param other SubunitCluster
     * @return true if the SubunitClusters are identical, false otherwise
     */
    public boolean isIdenticalByEntityIdTo(SubunitCluster other) {
        Subunit thisSub = subunitCluster.getSubunits().get(subunitCluster.getRepresentative());
        Subunit otherSub = other.getSubunits().get(other.getRepresentative());
        String thisName = thisSub.getName();
        String otherName = otherSub.getName();

        Structure thisStruct = thisSub.getStructure();
        Structure otherStruct = otherSub.getStructure();
        if (thisStruct == null || otherStruct == null) {
            SubunitClusterMerge.logger.info("SubunitClusters {}-{} have no referenced structures. Ignoring identity check by entity id",
                    thisName,
                    otherName);
            return false;
        }
        if (thisStruct != otherStruct) {
            // different object references: will not cluster even if entity id is same
            return false;
        }
        Chain thisChain = thisStruct.getChain(thisName);
        Chain otherChain = otherStruct.getChain(otherName);
        if (thisChain == null || otherChain == null) {
            SubunitClusterMerge.logger.info("Can't determine entity ids of SubunitClusters {}-{}. Ignoring identity check by entity id",
                    thisName,
                    otherName);
            return false;
        }
        if (thisChain.getEntityInfo() == null || otherChain.getEntityInfo() == null) {
            SubunitClusterMerge.logger.info("Can't determine entity ids of SubunitClusters {}-{}. Ignoring identity check by entity id",
                    thisName,
                    otherName);
            return false;
        }
        int thisEntityId = thisChain.getEntityInfo().getMolId();
        int otherEntityId = otherChain.getEntityInfo().getMolId();
        return thisEntityId == otherEntityId;
    }

    /**
     * Merges the other SubunitCluster into this one if it contains exactly the
     * same Subunit. This is checked by {@link #isIdenticalTo(SubunitCluster)}.
     *
     * @param other SubunitCluster
     * @return true if the SubunitClusters were merged, false otherwise
     */
    public boolean mergeIdentical(SubunitCluster other) {

        if (!isIdenticalTo(other))
            return false;

        SubunitClusterMerge.logger.info("SubunitClusters {}-{} are identical in sequence",
                subunitCluster.getSubunits().get(subunitCluster.getRepresentative()).getName(),
                other.getSubunits().get(other.getRepresentative()).getName());

        subunitCluster.getSubunits().addAll(other.getSubunits());
        subunitCluster.getSubunitEQR().addAll(other.getSubunitEQR());

        return true;
    }

    /**
     * Merges the other SubunitCluster into this one if it contains exactly the
     * same Subunit. This is checked by comparing the entity identifiers of the subunits
     * if one can be found.
     * Thus this only makes sense when the subunits are complete chains of a
     * deposited PDB entry.
     *
     * @param other SubunitCluster
     * @return true if the SubunitClusters were merged, false otherwise
     */
    public boolean mergeIdenticalByEntityId(SubunitCluster other) {

        if (!isIdenticalByEntityIdTo(other))
            return false;

        Subunit thisSub = subunitCluster.getSubunits().get(subunitCluster.getRepresentative());
        Subunit otherSub = other.getSubunits().get(other.getRepresentative());
        String thisName = thisSub.getName();
        String otherName = otherSub.getName();

        SubunitClusterMerge.logger.info("SubunitClusters {}-{} belong to same entity. Assuming they are identical",
                thisName,
                otherName);

        List<Integer> thisAligned = new ArrayList<Integer>();
        List<Integer> otherAligned = new ArrayList<Integer>();

        // we've merged by entity id, we can assume structure, chain and entity are available (checked in isIdenticalByEntityIdTo())
        Structure thisStruct = thisSub.getStructure();
        Structure otherStruct = otherSub.getStructure();
        Chain thisChain = thisStruct.getChain(thisName);
        Chain otherChain = otherStruct.getChain(otherName);
        EntityInfo entityInfo = thisChain.getEntityInfo();

        // Extract the aligned residues of both Subunits
        for (int thisIndex = 0; thisIndex < thisSub.size(); thisIndex++) {

            Group g = thisSub.getRepresentativeAtoms()[thisIndex].getGroup();

            int seqresIndex = entityInfo.getAlignedResIndex(g, thisChain);

            if (seqresIndex == -1) {
                // this might mean that FileParsingParameters.setAlignSeqRes() wasn't set to true during parsing
                continue;
            }

            // note the seqresindex is 1-based
            Group otherG = otherChain.getSeqResGroups().get(seqresIndex - 1);

            int otherIndex = otherChain.getAtomGroups().indexOf(otherG);
            if (otherIndex == -1) {
                // skip residues that are unobserved in other sequence ("gaps" in the entity SEQRES alignment)
                continue;
            }

            // Only consider residues that are part of the SubunitCluster
            if (subunitCluster.getSubunitEQR().get(subunitCluster.getRepresentative()).contains(thisIndex)
                    && other.getSubunitEQR().get(other.getRepresentative()).contains(otherIndex)) {
                thisAligned.add(thisIndex);
                otherAligned.add(otherIndex);
            }
        }

        if (thisAligned.size() == 0 && otherAligned.size() == 0) {
            SubunitClusterMerge.logger.warn("No equivalent aligned atoms found between SubunitClusters {}-{} via entity SEQRES alignment. Is FileParsingParameters.setAlignSeqRes() set?", thisName, otherName);
        }

        updateEquivResidues(other, thisAligned, otherAligned);

        return true;
    }

    /**
     * Merges the other SubunitCluster into this one if their representatives
     * sequences are similar (according to the criteria in params).
     * <p>
     * The sequence alignment is performed using linear {@link SimpleGapPenalty} and
     * BLOSUM62 as scoring matrix.
     *
     * @param other  SubunitCluster
     * @param params SubunitClustererParameters, with information whether to use local
     *               or global alignment, sequence identity and coverage thresholds.
     *               Threshold values lower than 0.7 are not recommended.
     *               Use {@link #mergeStructure} for lower values.
     * @return true if the SubunitClusters were merged, false otherwise
     * @throws CompoundNotFoundException
     */

    public boolean mergeSequence(SubunitCluster other, SubunitClustererParameters params) throws CompoundNotFoundException {
        Alignments.PairwiseSequenceAlignerType alignerType = Alignments.PairwiseSequenceAlignerType.LOCAL;
        if (params.isUseGlobalMetrics()) {
            alignerType = Alignments.PairwiseSequenceAlignerType.GLOBAL;
        }
        return mergeSequence(other, params, alignerType
                , new SimpleGapPenalty(),
                SubstitutionMatrixHelper.getBlosum62());
    }

    /**
     * Merges the other SubunitCluster into this one if their representatives
     * sequences are similar (according to the criteria in params).
     * <p>
     * The sequence alignment is performed using linear {@link SimpleGapPenalty} and
     * BLOSUM62 as scoring matrix.
     *
     * @param other       SubunitCluster
     * @param params      {@link SubunitClustererParameters}, with information whether to use local
     *                    or global alignment, sequence identity and coverage thresholds.
     *                    Threshold values lower than 0.7 are not recommended.
     *                    Use {@link #mergeStructure} for lower values.
     * @param alignerType parameter for the sequence alignment algorithm
     * @param gapPenalty  parameter for the sequence alignment algorithm
     * @param subsMatrix  parameter for the sequence alignment algorithm
     * @return true if the SubunitClusters were merged, false otherwise
     * @throws CompoundNotFoundException
     */

    public boolean mergeSequence(SubunitCluster other, SubunitClustererParameters params,
                                 Alignments.PairwiseSequenceAlignerType alignerType,
                                 GapPenalty gapPenalty,
                                 SubstitutionMatrix<AminoAcidCompound> subsMatrix)
            throws CompoundNotFoundException {

        // Extract the protein sequences as BioJava alignment objects
        ProteinSequence thisSequence = subunitCluster.getSubunits().get(subunitCluster.getRepresentative())
                .getProteinSequence();
        ProteinSequence otherSequence = other.getSubunits()
                .get(other.getRepresentative()).getProteinSequence();

        // Perform the alignment with provided parameters
        PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> aligner = Alignments
                .getPairwiseAligner(thisSequence, otherSequence, alignerType,
                        gapPenalty, subsMatrix);

        double sequenceIdentity;
        if (params.isUseGlobalMetrics()) {
            sequenceIdentity = aligner.getPair().getPercentageOfIdentity(true);
        } else {
            sequenceIdentity = aligner.getPair().getPercentageOfIdentity(false);
        }

        if (sequenceIdentity < params.getSequenceIdentityThreshold())
            return false;

        double sequenceCoverage = 0;
        if (params.isUseSequenceCoverage()) {
            // Calculate real coverage (subtract gaps in both sequences)
            double gaps1 = aligner.getPair().getAlignedSequence(1)
                    .getNumGapPositions();
            double gaps2 = aligner.getPair().getAlignedSequence(2)
                    .getNumGapPositions();
            double lengthAlignment = aligner.getPair().getLength();
            double lengthThis = aligner.getQuery().getLength();
            double lengthOther = aligner.getTarget().getLength();
            sequenceCoverage = (lengthAlignment - gaps1 - gaps2)
                    / Math.max(lengthThis, lengthOther);

            if (sequenceCoverage < params.getSequenceCoverageThreshold())
                return false;
        }

        SubunitClusterMerge.logger.info(String.format("SubunitClusters %s-%s are similar in sequence "
                        + "with %.2f sequence identity and %.2f coverage",
                subunitCluster.getSubunits().get(subunitCluster.getRepresentative()).getName(),
                other.getSubunits().get(other.getRepresentative()).getName(),
                sequenceIdentity, sequenceCoverage));

        // If coverage and sequence identity sufficient, merge other and this
        List<Integer> thisAligned = new ArrayList<Integer>();
        List<Integer> otherAligned = new ArrayList<Integer>();

        // Extract the aligned residues of both Subunit
        for (int p = 1; p < aligner.getPair().getLength() + 1; p++) {

            // Skip gaps in any of the two sequences
            if (aligner.getPair().getAlignedSequence(1).isGap(p))
                continue;
            if (aligner.getPair().getAlignedSequence(2).isGap(p))
                continue;

            int thisIndex = aligner.getPair().getIndexInQueryAt(p) - 1;
            int otherIndex = aligner.getPair().getIndexInTargetAt(p) - 1;

            // Only consider residues that are part of the SubunitCluster
            if (subunitCluster.getSubunitEQR().get(subunitCluster.getRepresentative()).contains(thisIndex)
                    && other.getSubunitEQR().get(other.getRepresentative()).contains(otherIndex)) {
                thisAligned.add(thisIndex);
                otherAligned.add(otherIndex);
            }
        }

        updateEquivResidues(other, thisAligned, otherAligned);

        subunitCluster.setMethod(SubunitClustererMethod.SEQUENCE);
        subunitCluster.setPseudoStoichiometric(!params.isHighConfidenceScores(sequenceIdentity, sequenceCoverage));

        return true;
    }

    /**
     * Merges the other SubunitCluster into this one if their representative
     * Atoms are structurally similar (according to the criteria in params).
     * <p>
     *
     * @param other  SubunitCluster
     * @param params {@link SubunitClustererParameters}, with information on what alignment
     *               algorithm to use, RMSD/TMScore and structure coverage thresholds.
     * @return true if the SubunitClusters were merged, false otherwise
     * @throws StructureException
     */

    public boolean mergeStructure(SubunitCluster other, SubunitClustererParameters params) throws StructureException {

        StructureAlignment aligner = StructureAlignmentFactory.getAlgorithm(params.getSuperpositionAlgorithm());
        ConfigStrucAligParams aligner_params = aligner.getParameters();

        Method setOptimizeAlignment = null;
        try {
            setOptimizeAlignment = aligner_params.getClass().getMethod("setOptimizeAlignment", boolean.class);
        } catch (NoSuchMethodException e) {
            //alignment algorithm does not have an optimization switch, moving on
        }
        if (setOptimizeAlignment != null) {
            try {
                setOptimizeAlignment.invoke(aligner_params, params.isOptimizeAlignment());
            } catch (IllegalAccessException | InvocationTargetException e) {
                SubunitClusterMerge.logger.warn("Could not set alignment optimisation switch");
            }
        }

        AFPChain afp = aligner.align(subunitCluster.getSubunits().get(subunitCluster.getRepresentative())
                        .getRepresentativeAtoms(),
                other.getSubunits().get(other.getRepresentative())
                        .getRepresentativeAtoms());

        // Convert AFPChain to MultipleAlignment for convenience
        MultipleAlignment msa = new MultipleAlignmentEnsembleImpl(
                afp,
                subunitCluster.getSubunits().get(subunitCluster.getRepresentative()).getRepresentativeAtoms(),
                other.getSubunits().get(other.getRepresentative())
                        .getRepresentativeAtoms(), false)
                .getMultipleAlignment(0);

        double structureCoverage = Math.min(msa.getCoverages().get(0), msa
                .getCoverages().get(1));

        if (params.isUseStructureCoverage() && structureCoverage < params.getStructureCoverageThreshold()) {
            return false;
        }

        double rmsd = afp.getTotalRmsdOpt();
        if (params.isUseRMSD() && rmsd > params.getRMSDThreshold()) {
            return false;
        }

        double tmScore = afp.getTMScore();
        if (params.isUseTMScore() && tmScore < params.getTMThreshold()) {
            return false;
        }

        SubunitClusterMerge.logger.info(String.format("SubunitClusters are structurally similar with "
                + "%.2f RMSD %.2f coverage", rmsd, structureCoverage));

        // Merge clusters
        List<List<Integer>> alignedRes = msa.getBlock(0).getAlignRes();
        List<Integer> thisAligned = new ArrayList<Integer>();
        List<Integer> otherAligned = new ArrayList<Integer>();

        // Extract the aligned residues of both Subunit
        for (int p = 0; p < msa.length(); p++) {

            // Skip gaps in any of the two sequences
            if (alignedRes.get(0).get(p) == null)
                continue;
            if (alignedRes.get(1).get(p) == null)
                continue;

            int thisIndex = alignedRes.get(0).get(p);
            int otherIndex = alignedRes.get(1).get(p);

            // Only consider residues that are part of the SubunitCluster
            if (subunitCluster.getSubunitEQR().get(subunitCluster.getRepresentative()).contains(thisIndex)
                    && other.getSubunitEQR().get(other.getRepresentative()).contains(
                    otherIndex)) {
                thisAligned.add(thisIndex);
                otherAligned.add(otherIndex);
            }
        }

        updateEquivResidues(other, thisAligned, otherAligned);

        subunitCluster.setMethod(SubunitClustererMethod.STRUCTURE);
        subunitCluster.setPseudoStoichiometric(true);

        return true;
    }

    void updateEquivResidues(SubunitCluster other, List<Integer> thisAligned, List<Integer> otherAligned) {
        // Do a List intersection to find out which EQR columns to remove
        List<Integer> thisRemove = new ArrayList<Integer>();
        List<Integer> otherRemove = new ArrayList<Integer>();

        for (int t = 0; t < subunitCluster.getSubunitEQR().get(subunitCluster.getRepresentative()).size(); t++) {
            // If the index is aligned do nothing, otherwise mark as removing
            if (!thisAligned.contains(subunitCluster.getSubunitEQR().get(subunitCluster.getRepresentative()).get(t)))
                thisRemove.add(t);
        }

        for (int t = 0; t < other.getSubunitEQR().get(other.getRepresentative()).size(); t++) {
            // If the index is aligned do nothing, otherwise mark as removing
            if (!otherAligned.contains(other.getSubunitEQR().get(other.getRepresentative()).get(t)))
                otherRemove.add(t);
        }
        // Now remove unaligned columns, from end to start
        Collections.sort(thisRemove);
        Collections.reverse(thisRemove);
        Collections.sort(otherRemove);
        Collections.reverse(otherRemove);

        for (int t = 0; t < thisRemove.size(); t++) {
            for (List<Integer> eqr : subunitCluster.getSubunitEQR()) {
                int column = thisRemove.get(t);
                eqr.remove(column);
            }
        }

        for (int t = 0; t < otherRemove.size(); t++) {
            for (List<Integer> eqr : other.getSubunitEQR()) {
                int column = otherRemove.get(t);
                eqr.remove(column);
            }
        }

        // The representative is the longest sequence
        if (subunitCluster.getSubunits().get(subunitCluster.getRepresentative()).size() < other.getSubunits().get(other.getRepresentative()).size())
            subunitCluster.setRepresentative(other.getRepresentative() + subunitCluster.getSubunits().size());

        subunitCluster.getSubunits().addAll(other.getSubunits());
        subunitCluster.getSubunitEQR().addAll(other.getSubunitEQR());

    }
}