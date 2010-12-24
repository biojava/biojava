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
package org.biojava.bio.structure.io.mmcif.chem;

import java.io.Serializable;


/**
 * Enumerates the possible classifications of residues.
 * This information is derived from the mmcif dictionary.
 * @author mulvaney
 * @author Andreas Prlic
 * @see <a href="http://mmcif.rcsb.org/dictionaries/mmcif_pdbx.dic/Items/_chem_comp.type.html">link into mmCIF dictionary</a>
 * @since 1.7
 */

public enum ResidueType implements Serializable {

   atomn(null, "null"), // present in db for _chem_comp.id_ = 'CFL' but not enumerated in dictionary
   dPeptideLinking(PolymerType.dpeptide, "D-peptide linking"),
   lPeptideLinking(PolymerType.peptide, "L-peptide linking"),
   glycine(PolymerType.peptide,"PEPTIDE LINKING"),
   dPeptideAminoTerminus(PolymerType.dpeptide, "D-peptide NH3 amino terminus"),
   lPeptideAminoTerminus(PolymerType.peptide, "L-peptide NH3 amino terminus"),
   dPeptideCarboxyTerminus(PolymerType.dpeptide, "D-peptide COOH carboxy terminus"),
   lPeptideCarboxyTerminus(PolymerType.peptide, "L-peptide COOH carboxy terminus"),
   dnaLinking(PolymerType.dna, "DNA linking"),
   rnaLinking(PolymerType.rna, "RNA linking"),
   dna3PrimeTerminus(PolymerType.dna, "DNA OH 3 prime terminus"),
   rna3PrimeTerminus(PolymerType.rna, "RNA OH 3 prime terminus"),
   dna5PrimeTerminus(PolymerType.dna, "DNA OH 5 prime terminus"),
   rna5PrimeTerminus(PolymerType.rna, "RNA OH 5 prime terminus"),
   dSaccharide(PolymerType.polysaccharide, "D-saccharide"),
   dSaccharide14and14linking(PolymerType.polysaccharide, "D-saccharide 1,4 and 1,4 linking"),
   dSaccharide14and16linking(PolymerType.polysaccharide, "D-saccharide 1,4 and 1,6 linking"),
   lSaccharide(PolymerType.lpolysaccharide, "L-saccharide"),
   lSaccharide14and14linking(PolymerType.lpolysaccharide, "L-saccharide 1,4 and 1,4 linking"),
   lSaccharide14and16linking(PolymerType.lpolysaccharide, "L-saccharide 1,4 and 1,6 linking"),
   saccharide(PolymerType.polysaccharide, "saccharide"),
   nonPolymer(null, "non-polymer"),
   otherChemComp(null, "other");

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
    * String value of the type
    */
   public final String chem_comp_type;

   public static ResidueType getResidueTypeFromString(String chem_comp_type)
   {
	   
	   chem_comp_type = chem_comp_type.replaceAll("'", "");
	   chem_comp_type = chem_comp_type.replaceAll("\"", "");
	   
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
