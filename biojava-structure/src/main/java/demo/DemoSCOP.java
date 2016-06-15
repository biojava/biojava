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
 * Created on Jun 30, 2010
 * Author: Andreas Prlic
 *
 */

package demo;


import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.ce.CeMain;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.AfpChainWriter;
import org.biojava.nbio.structure.align.util.AFPChainScorer;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.scop.*;

import java.util.List;

/** A class demonstrating the use of the SCOP parsing tools
 *
 * @author Andreas Prlic
 *
 */
public class DemoSCOP
{
	public static void main(String[] args){

		DemoSCOP demo = new DemoSCOP();


		// this creates a local copy of SCOP
		ScopDatabase scop = new ScopInstallation();

		// an alternative would be this one, which fetches data dynamically
		//ScopDatabase scop = new RemoteScopInstallation();

		ScopFactory.setScopDatabase(scop);

		demo.getCategories();
		demo.printDomainsForPDB();
		demo.traverseHierarchy();
		demo.alignSuperfamily();
	}

	/** Traverse throught the SCOP hierarchy
	 *
	 */
	public void traverseHierarchy()
	{
		String pdbId = "4HHB";
		// download SCOP if required and load into memory
		ScopDatabase scop = ScopFactory.getSCOP();

		List<ScopDomain> domains = scop.getDomainsForPDB(pdbId);

		// show the hierachy for the first domain:

		ScopNode node = scop.getScopNode(domains.get(0).getSunid());

		while (node != null){

			System.out.println("This node: sunid:" + node.getSunid() );
			System.out.println(scop.getScopDescriptionBySunid(node.getSunid()));
			node = scop.getScopNode(node.getParentSunid());
		}

	}

	/** Get various categories
	 *
	 */
	public void getCategories(){
		// download SCOP if required and load into memory
		ScopDatabase scop = ScopFactory.getSCOP();
		List<ScopDescription> superfams = scop.getByCategory(ScopCategory.Superfamily);

		System.out.println("Total nr. of superfamilies:" + superfams.size());

		List<ScopDescription> folds = scop.getByCategory(ScopCategory.Fold);
		System.out.println("Total nr. of folds:" + folds.size());

	}

	public void alignSuperfamily(){
		// download SCOP if required and load into memory
		ScopDatabase scop = ScopFactory.getSCOP();
		List<ScopDescription> superfams = scop.getByCategory(ScopCategory.Superfamily);

		System.out.println("Total nr. of superfamilies:" + superfams.size());


		// configure where to load PDB files from and
		// what information to load
		AtomCache cache = new AtomCache();
		FileParsingParameters fileparams = new FileParsingParameters() ;
		fileparams.setAlignSeqRes(false);
		fileparams.setParseSecStruc(false);
		cache.setFileParsingParams(fileparams);

		// get tge first superfamily
		ScopDescription superfam1 = superfams.get(0);
		System.out.println("First superfamily: " + superfam1);

		ScopNode node = scop.getScopNode(superfam1.getSunID());
		System.out.println("scopNode for first superfamily:" + node);

		List<ScopDomain> doms4superfam1 = scop.getScopDomainsBySunid(superfam1.getSunID());
		ScopDomain dom1 = doms4superfam1.get(0);

		// align the first domain against all others members of this superfamily
		for ( int i = 1 ; i < doms4superfam1.size() ; i ++){

			ScopDomain dom2 = doms4superfam1.get(i);

			try {
				Structure s1 = cache.getStructureForDomain(dom1);
				Structure s2 = cache.getStructureForDomain(dom2);

				Atom[] ca1 = StructureTools.getAtomCAArray(s1);
				Atom[] ca2 = StructureTools.getAtomCAArray(s2);
				StructureAlignment ce = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
				AFPChain afpChain = ce.align(ca1, ca2);

				//System.out.println(afpChain.toCE(ca1, ca2));

				//StructureAlignmentDisplay.display(afpChain, ca1, ca2);

				System.out.println(dom1.getScopId() + " vs. " + dom2.getScopId()+ " :" + afpChain.getProbability());
				double tmScore = AFPChainScorer.getTMScore(afpChain, ca1, ca2);
				afpChain.setTMScore(tmScore);
				System.out.println(AfpChainWriter.toScoresList(afpChain));

			} catch (Exception e){
				e.printStackTrace();
			}
		}

	}

	public void printDomainsForPDB(){
		String pdbId = "4HHB";

		// download SCOP if required and load into memory
		ScopDatabase scop = ScopFactory.getSCOP();

		List<ScopDomain> domains = scop.getDomainsForPDB(pdbId);

		System.out.println(domains);

	}


}
