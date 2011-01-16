package org.biojava.bio.structure.domain;

import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.domain.pdp.ClusterDomains;
import org.biojava.bio.structure.domain.pdp.CutDomain;
import org.biojava.bio.structure.domain.pdp.CutSites;
import org.biojava.bio.structure.domain.pdp.Domain;
import org.biojava.bio.structure.domain.pdp.GetDistanceMatrix;
import org.biojava.bio.structure.domain.pdp.PDPDistanceMatrix;
import org.biojava.bio.structure.domain.pdp.Segment;
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
 *
 */
public class ProteinDomainParser {


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
		dom.size = ca.length;
		dom.nseg=1;
		dom.getSegment(0).setFrom(0);
		dom.getSegment(0).setTo(ca.length-1);

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

	public static final void listdomains(List<Domain> domains){

		int i = -1;
		for ( Domain dom : domains){
			i++;
			System.out.println("DOMAIN:" + i + " size:" + dom.size + " " +  dom.score);
			List<Segment> segments = dom.getSegments();

			for ( Segment s : segments){
				System.out.println("   Segment: " + s);

			}
		}
	}

}
