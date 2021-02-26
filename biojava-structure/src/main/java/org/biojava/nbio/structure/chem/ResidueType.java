package org.biojava.nbio.structure.chem;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * Enumerates the possible classifications of residues. These are generally more specific than PolymerTypes
 * This information is derived from the mmcif dictionary.
 * @author mulvaney
 * @author Andreas Prlic
 * @see <a href="http://mmcif.rcsb.org/dictionaries/mmcif_pdbx.dic/Items/_chem_comp.type.html">link into mmCIF dictionary</a>
 * @since 1.7
 */
public enum ResidueType implements Serializable {
    atomn(null, "null"), // present in db for _chem_comp.id_ = 'CFL' but not enumerated in dictionary
    // Peptides
    dPeptideLinking(PolymerType.dpeptide, "D-peptide linking"),
    lPeptideLinking(PolymerType.peptide, "L-peptide linking"),
    glycine(PolymerType.peptide,"PEPTIDE LINKING"),
    peptideLike(PolymerType.otherPolymer, "peptide-like"),
    dPeptideAminoTerminus(PolymerType.dpeptide, "D-peptide NH3 amino terminus"),
    lPeptideAminoTerminus(PolymerType.peptide, "L-peptide NH3 amino terminus"),
    dPeptideCarboxyTerminus(PolymerType.dpeptide, "D-peptide COOH carboxy terminus"),
    lPeptideCarboxyTerminus(PolymerType.peptide, "L-peptide COOH carboxy terminus"),
    // Nucleotides
    dnaLinking(PolymerType.dna, "DNA linking"),
    rnaLinking(PolymerType.rna, "RNA linking"),
    dna3PrimeTerminus(PolymerType.dna, "DNA OH 3 prime terminus"),
    rna3PrimeTerminus(PolymerType.rna, "RNA OH 3 prime terminus"),
    dna5PrimeTerminus(PolymerType.dna, "DNA OH 5 prime terminus"),
    rna5PrimeTerminus(PolymerType.rna, "RNA OH 5 prime terminus"),
    // Sugars
    dSaccharide(PolymerType.polysaccharide, "D-saccharide"),
    dSaccharide14and14linking(PolymerType.polysaccharide, "D-saccharide 1,4 and 1,4 linking"),
    dSaccharide14and16linking(PolymerType.polysaccharide, "D-saccharide 1,4 and 1,6 linking"),
    lSaccharide(PolymerType.lpolysaccharide, "L-saccharide"),
    lSaccharide14and14linking(PolymerType.lpolysaccharide, "L-saccharide 1,4 and 1,4 linking"),
    lSaccharide14and16linking(PolymerType.lpolysaccharide, "L-saccharide 1,4 and 1,6 linking"),
    saccharide(PolymerType.polysaccharide, "saccharide"),
    // Iso-peptides
    dBetaPeptideCGammaLinking(PolymerType.dpeptide,"D-beta-peptide, C-gamma linking"),
    dGammaPeptideCDeltaLinking(PolymerType.dpeptide,"D-gamma-peptide, C-delta linking"),
    lBetaPeptideCGammaLinking(PolymerType.peptide,"L-beta-peptide, C-gamma linking"),
    lGammaPeptideCDeltaLinking(PolymerType.peptide,"L-gamma-peptide, C-delta linking"),
    // L nucleotides. As of 2015-04, these are only found in D-DNA hybrids, so they don't have their own PolymerType
    lDNALinking(PolymerType.dna,"L-DNA linking"),
    lRNALinking(PolymerType.dna,"L-RNA linking"),
    // Other
    nonPolymer(null, "non-polymer"),
    otherChemComp(null, "other");

    static Map<String, ResidueType> lookupTable = new HashMap<>();

    static {
        for (ResidueType residueType : ResidueType.values() ) {
            lookupTable.put(residueType.chem_comp_type, residueType);
            lookupTable.put(residueType.chem_comp_type.toLowerCase(), residueType);
        }
    }

    ResidueType(PolymerType polymerType, String chem_comp_type) {
        this.polymerType = polymerType;
        this.chem_comp_type = chem_comp_type;
    }

    /**
     * The associated {@link PolymerType}
     */
    public final PolymerType polymerType;

    /**
     * Gets the associated PolymerType, which are less specific
     * @return
     */
    public PolymerType getPolymerType() {
        return polymerType;
    }

    /**
     * String value of the type
     */
    public final String chem_comp_type;

    /** Get ResidueType by chem_comp_type
     *
     * @param chem_comp_type e.g. L-peptide linking
     * @return
     */
    public static ResidueType getResidueTypeFromString(String chem_comp_type) {
        if (chem_comp_type == null) {
            return null;
        }

        // Almost all calls to this method are for L-peptide linking. Use this knowledge for a shortcut.
        if (chem_comp_type.equalsIgnoreCase(lPeptideLinking.chem_comp_type)) {
            return lPeptideLinking;
        }

        ResidueType lookedUpResidueType = lookupTable.get(chem_comp_type);
        if (lookedUpResidueType != null) {
            return lookedUpResidueType;
        }

        /*
         * Unfortunately it can be guaranteed that chem_comp_type case sensitivity is preserved.
         * E.g. mmtf has it all upper-case. As such we need to do a second check
         */
        lookedUpResidueType = lookupTable.get(chem_comp_type.toLowerCase());
        if (lookedUpResidueType != null) {
            return lookedUpResidueType;
        }

        // preserving previous behaviour. Not sure if this is really necessary?
        for (ResidueType residueType : ResidueType.values()) {
            if(residueType.chem_comp_type.equalsIgnoreCase(chem_comp_type)) {
                return residueType;
            }

            if (residueType.chem_comp_type.toLowerCase().startsWith(chem_comp_type.toLowerCase())) {
                return residueType;
            }
            if (chem_comp_type.toLowerCase().startsWith(residueType.chem_comp_type.toLowerCase())) {
                return residueType;
            }
        }
        return null;
    }
}
