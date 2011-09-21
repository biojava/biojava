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
 * Author: Andreas Prlic 
 *
 */
package org.biojava.bio.structure.domain;

import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.ResidueNumber;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.domain.pdp.ClusterDomains;
import org.biojava.bio.structure.domain.pdp.CutDomain;
import org.biojava.bio.structure.domain.pdp.CutSites;
import org.biojava.bio.structure.domain.pdp.Domain;
import org.biojava.bio.structure.domain.pdp.GetDistanceMatrix;
import org.biojava.bio.structure.domain.pdp.PDPDistanceMatrix;
import org.biojava.bio.structure.domain.pdp.ShortSegmentRemover;


/** Protein Domain Parser is a an algorithm that attempts at assigning domains for 3D protein structures.
 * Since domains in proteins are difficult to define, results detected by automated algorithms have to be taken with a grain of salt.
 *
 *   see<pre>
J Mol Biol. 2004 Jun 4;339(3):647-78.

Toward consistent assignment of structural domains in proteins.

Veretnik S, Bourne PE, Alexandrov NN, Shindyalov IN.
</pre>

 * This implementation is based on a Java port of the PDP algorithm, as described in:
 *   
 * 
 * @author Andreas Prlic
 * @since 3.0.2
 */
public class LocalProteinDomainParser {


	/** make sure this class can only get accessed via the static method calls
	 * 
	 */
	private LocalProteinDomainParser(){
		
	}
	
	/** Suggest domains for a protein structure
	 * 
	 * @param s the protein structure
	 * @return a list of possible domains
	 * @throws StructureException
	 */
	public static List<Domain> suggestDomains(Structure s) throws StructureException{

		Atom[] ca = StructureTools.getAtomCAArray(s);

		return suggestDomains(ca);
	} 


	/** Suggest domains for a set of Calpha atoms
	 * 
	 * @param ca an array of Calpha atoms
	 * @return a list of possible domains
	 * @throws StructureException
	 */
	public static List<Domain> suggestDomains(Atom[] ca) throws StructureException{

		GetDistanceMatrix distMaxCalculator = new GetDistanceMatrix();

		PDPDistanceMatrix pdpMatrix = distMaxCalculator.getDistanceMatrix(ca);



		Domain dom = new Domain();
		Chain c = ca[0].getGroup().getChain();
		dom.setId("D"+c.getParent().getPDBCode()+c.getId()+"1");
		dom.setSize(ca.length);
		dom.setNseg(1);
		dom.getSegmentAtPos(0).setFrom(0);
		dom.getSegmentAtPos(0).setTo(ca.length-1);

		CutSites cutSites = new CutSites();

		// Do the initial splitting
		CutDomain cutDomain = new CutDomain(ca,pdpMatrix);
		cutDomain.cutDomain(dom, cutSites, pdpMatrix);
		List<Domain> domains =  cutDomain.getDomains();

		
		//
		domains = ClusterDomains.cluster(domains, pdpMatrix);

		ShortSegmentRemover.cleanup(domains);


		return domains;

	}

	

}
