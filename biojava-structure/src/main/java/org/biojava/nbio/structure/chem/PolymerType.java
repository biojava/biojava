package org.biojava.nbio.structure.chem;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Enumerates the classification of polymers.
 * This information is derived from the mmcif dictionary
 * @author mulvaney
 * @author Andreas Prlic
 * @see <a href="http://mmcif.rcsb.org/dictionaries/mmcif_pdbx.dic/Items/_entity_poly.type.html">link into mmCIF dictionary</a>
 * @since 1.7
 */
public enum PolymerType implements Serializable {
    /**
     * polypeptide(L)
     */
    peptide("polypeptide(L)"),
    /**
     * polypeptide(D)
     */
    dpeptide("polypeptide(D)"),
    /**
     * polydeoxyribonucleotide
     */
    dna("polydeoxyribonucleotide"),
    /**
     * polyribonucleotide
     */
    rna("polyribonucleotide"),
    /**
     * polydeoxyribonucleotide/polyribonucleotide hybrid
     */
    dnarna("polydeoxyribonucleotide/polyribonucleotide hybrid"),
    /**
     * polysaccharide(D)
     */
    polysaccharide("polysaccharide(D)"),
    /**
     * polysaccharide(L)
     */
    lpolysaccharide("polysaccharide(L)"),
    /**
     * other
     */
    otherPolymer("other"),
    /**
     * cyclic peptides
     */
    cyclicPeptide("cyclic-pseudo-peptide"),
    /**
     * Peptide nucleic acids
     */
    peptideNucleicAcid("peptide nucleic acid"),
    /**
     * if all else fails...
     */
    unknown(null);

    static Map<String, PolymerType> lookupTable = new HashMap<>();

    static {
        for (PolymerType polymerType : PolymerType.values()) {
            if (polymerType == unknown) {
                continue;
            }

            lookupTable.put(polymerType.entity_poly_type,polymerType);
            lookupTable.put(polymerType.entity_poly_type.toLowerCase(), polymerType);
        }
    }

    public final String entity_poly_type;

    PolymerType(String entity_poly_type) {
        this.entity_poly_type = entity_poly_type;
    }

    public static PolymerType polymerTypeFromString(String polymerTypeString) {
        if (polymerTypeString.equalsIgnoreCase(peptide.entity_poly_type)) {
            return peptide;
        }

        PolymerType lookedUpPolymerType = lookupTable.get(polymerTypeString);
        if (lookedUpPolymerType != null) {
            return lookedUpPolymerType;
        }

        lookedUpPolymerType = lookupTable.get(polymerTypeString.toLowerCase());
        if (lookedUpPolymerType != null) {
            return lookedUpPolymerType;
        }

        for (PolymerType polymerType : PolymerType.values()) {
            if (polymerTypeString.equals(polymerType.entity_poly_type)) {
                return polymerType;
            }
        }

        return unknown;
    }

    /**
     * Convenience <tt>Set</tt> of polymer types classified as protein.  This only contains {@link #peptide}
     */
    public static final Set<PolymerType> PROTEIN_ONLY;

    /**
     * Convenience <tt>Set</tt> of polymer types classified as DNA.  This only contains {@link #dna}
     */
    public static final Set<PolymerType> DNA_ONLY;

    /**
     * Convenience <tt>Set</tt> of polymer types classified as RNA.  This only contains {@link #rna}
     */
    public static final Set<PolymerType> RNA_ONLY;

    /**
     * Convenience <tt>Set</tt> of polymer types classified as DNA.  This contains:
     * <ul>
     * <li>{@link #dna}</li>
     * <li>{@link #rna}</li>
     * <li>{@link #dnarna}</li>
     * </ul>
     */
    public static final Set<PolymerType> POLYNUCLEOTIDE_ONLY;

    /**
     * Convenience <tt>Set</tt> of all polymer types.
     */
    public static final Set<PolymerType> ALL_POLYMER_TYPES;

    static {
        Set<PolymerType> tmp;

        tmp = new HashSet<>();
        tmp.add(peptide);
        PROTEIN_ONLY = Collections.unmodifiableSet(tmp);

        tmp = new HashSet<>();
        tmp.add(dna);
        DNA_ONLY = Collections.unmodifiableSet(tmp);

        tmp = new HashSet<>();
        tmp.add(rna);
        RNA_ONLY = Collections.unmodifiableSet(tmp);

        tmp = new HashSet<>();
        tmp.add(dna);
        tmp.add(rna);
        tmp.add(dnarna);
        POLYNUCLEOTIDE_ONLY = Collections.unmodifiableSet(tmp);

        ALL_POLYMER_TYPES = Collections.unmodifiableSet(new HashSet<>(Arrays.asList(values())));
    }
}
