/* This class is based on the original FATCAT implementation by
 * <pre>
 * Yuzhen Ye & Adam Godzik (2003)
 * Flexible structure alignment by chaining aligned fragment pairs allowing twists.
 * Bioinformatics vol.19 suppl. 2. ii246-ii255.
 * http://www.ncbi.nlm.nih.gov/pubmed/14534198
 * </pre>
 * 
 * Thanks to Yuzhen Ye and A. Godzik for granting permission to freely use and redistribute this code.
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
 *
 * Created by Andreas Prlic - RCSB PDB 
 * 
 */
package org.biojava.bio.structure.align.fatcat;

import org.biojava.bio.structure.Atom;

import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.bio.structure.align.fatcat.calc.FatCatParameters;

import org.biojava.bio.structure.align.model.AFPChain;


public class FatCatRigid extends FatCat implements StructureAlignment{

	public static String algorithmName = "jFatCat_rigid";

	FatCatParameters params;

	public FatCatRigid(){
		super();
		params = new FatCatParameters();
	    params.setMaxTra(0);		
	}

	public AFPChain align(Atom[] ca1, Atom[] ca2) throws StructureException {

		AFPChain afpChain = alignRigid(ca1, ca2, params);
		afpChain.setAlgorithmName(algorithmName);
		afpChain.setVersion(VERSION+"");
		return afpChain;
	}

	public AFPChain align(Atom[] ca1, Atom[] ca2, Object param)
	throws StructureException {

		if ( ! (param instanceof FatCatParameters)){
			throw new IllegalArgumentException("FatCat algorithm needs FatCatParameters object as argument.");
		}

		params = (FatCatParameters) param;		

		AFPChain afpChain= alignRigid(ca1, ca2, params);
		afpChain.setAlgorithmName(algorithmName);
		afpChain.setVersion(VERSION+"");
		return afpChain;
	}

	public String getAlgorithmName() {

		return algorithmName;
	}

	public ConfigStrucAligParams getParameters() {

		return params;
	}

	public String getVersion(){
		return VERSION+"";
	}

//	public StructureAlignmentJmol display(AFPChain afpChain, Atom[] ca1,
//			Atom[] ca2, List<Group> hetatms, List<Group> nucs,
//			List<Group> hetatms2, List<Group> nucs2) throws StructureException {
//
//		StructureAlignmentJmol gui =  super.display(afpChain, ca1, ca2, hetatms, nucs, hetatms2, nucs2);
//		gui.setTitle(getAlgorithmName() + " : " + afpChain.getName1() + " vs. " + afpChain.getName2());
//		return gui;
//	}

	public void setParameters(ConfigStrucAligParams parameters) {
		if (! (parameters instanceof FatCatParameters)){
			throw new IllegalArgumentException("Provided parameters are not of type FatCatParameters!");
		}
		params = (FatCatParameters) parameters;
	}



}
