/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.structure;

import org.biojava.nbio.structure.io.mmcif.chem.PolymerType;
import org.biojava.nbio.structure.io.mmcif.chem.ResidueType;

import java.util.*;


/**
 * This contains basic categories for Group types. It reflects the categorization
 * used in old PDB file (e.g. for storing whether a residue is composed of
 * ATOM or HETATM records. It is less specific than the mmCIF/PDBx-defined
 * ResidueType enum, which may be more suitable for future applications.
 *
 * @author Andreas Prlic
 * @author Spencer Bliven
 * @since 1.7
 * @see ResidueType
 *
 */
public enum GroupType {


	/**
	 * The type for amino acids (L-peptides)
	 */
	AMINOACID("amino",GroupType.matchPolymerTypes(PolymerType.PROTEIN_ONLY)),

	/**
	 * The type for nucleotide groups (dna and rna)
	 */
	NUCLEOTIDE("nucleotide",GroupType.matchPolymerTypes(PolymerType.POLYNUCLEOTIDE_ONLY)),

	/**
	 * The type for hetero groups (everything else)
	 */
	HETATM("hetatm",GroupType.getHetatmTypes());

	private final String name;
	private final Set<ResidueType> types;
	private GroupType(String name,Set<ResidueType> types) {
		this.name = name;
		this.types = types;
	}

	/**
	 * The 3-letter codes used in the PDB to identify water molecules
	 * @see Group#isWater()
	 */
	public static final List<String> WATERNAMES = Arrays.asList("HOH", "DOD", "WAT");

	/**
	 * @return The name of this GroupType. One of "amino", "nucleotide", or "hetatm"
	 */
	@Override
	public String toString() {
		return name;
	}

	/**
	 * Get a set of ResidueTypes loosely equivalent to this GroupType.
	 *
	 * <p>Because mmCIF and PDB handle modified residues differently, some
	 * Groups may have a well-defined ResidueType yet still be HETATMs.
	 * @return A Set of ResidueTypes commonly classified as this GroupType
	 */
	public Set<ResidueType> getResidueTypes() {
		return types;
	}

	/**
	 * Get ResidueTypes whose polymerType is contained in a certain set.
	 * This is used for defining the AMINOACID and NUCLEOTIDE sets.
	 * @param allowedTypes
	 * @return
	 */
	private static Set<ResidueType> matchPolymerTypes(Set<PolymerType> allowedTypes) {
		Set<ResidueType> matched = new HashSet<ResidueType>();
		for(ResidueType restype : ResidueType.values()) {
			if(allowedTypes.contains(restype.polymerType)) {
				matched.add(restype);
			}
		}
		return Collections.unmodifiableSet(matched);
	}

	/**
	 * Bundles everything not in AMINOACID or NUCLEOTIDE into the HETATM entry
	 * @return
	 */
	private static Set<ResidueType> getHetatmTypes() {
		Set<ResidueType> unmatched = new HashSet<ResidueType>();
		for(ResidueType restype : ResidueType.values()) {
			if(!AMINOACID.getResidueTypes().contains(restype) &&
					!NUCLEOTIDE.getResidueTypes().contains(restype) ) {
				unmatched.add(restype);
			}
		}
		return Collections.unmodifiableSet(unmatched);
	}
}
