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
 *
 */
package org.biojava.nbio.structure.io.mmcif.chem;

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


	static Map<String,ResidueType> lookupTable = new HashMap<>();

	static {

		for (ResidueType rt : ResidueType.values() ) {
				lookupTable.put(rt.chem_comp_type,rt);
				lookupTable.put(rt.chem_comp_type.toLowerCase(),rt);
		}
	}
	ResidueType(PolymerType pt, String chem_comp_type)
	{
		this.polymerType = pt;
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
	public PolymerType getPolymerType() {return polymerType;}

	/**
	 * String value of the type
	 */
	public final String chem_comp_type;

	/** Get ResidueType by chem_comp_type
	 *
	 * @param chem_comp_type e.g. L-peptide linking
	 * @return
	 */
	public static ResidueType getResidueTypeFromString(String chem_comp_type)
	{

		// Almost all calls to this method are for L-peptide linking. Use this knowledge for a shortcut.

		if ( chem_comp_type.equalsIgnoreCase(lPeptideLinking.chem_comp_type) )
			return lPeptideLinking;

		ResidueType rtype = lookupTable.get(chem_comp_type);
		if ( rtype != null)
			return rtype;

		/** Unfortunately it can be guaranteed that chem_comp_type case sensitivity is preserved.
		 * E.g. mmtf has it all upper-case. As such we need to do a second check
		 */
		rtype = lookupTable.get(chem_comp_type.toLowerCase());
		if ( rtype != null)
			return rtype;



		// preserving previous behaviour. Not sure if this is really necessary?
		for(ResidueType rt : ResidueType.values())
		{
			if(rt.chem_comp_type.equalsIgnoreCase(chem_comp_type))
			{
				return rt;
			}
			if ( rt.chem_comp_type.startsWith(chem_comp_type))
				return rt;
			if ( chem_comp_type.startsWith(rt.chem_comp_type))
				return rt;
		}
		return null;
	}
}
