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
 * created at 14 Oct 2013
 * Author: ap3
 */

package demo;

import org.biojava.nbio.structure.align.client.StructureName;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.cath.CathDatabase;
import org.biojava.nbio.structure.cath.CathDomain;
import org.biojava.nbio.structure.cath.CathInstallation;
import org.biojava.nbio.structure.cath.CathSegment;
import org.biojava.nbio.structure.gui.BiojavaJmol;
import org.biojava.nbio.structure.StructureIO;

import java.util.List;

public class DemoShowCATHDomain {


	private static final String DEFAULT_SCRIPT ="select * ; cartoon on; spacefill off; wireframe off; select ligands; wireframe on; spacefill on;";
	private static final String[] colors = new String[]{"red","green","blue","yellow"};

	public static void main(String[] args){

		UserConfiguration config = new UserConfiguration();
		config.setPdbFilePath("/tmp/");

		String pdbID = "1DAN";

		CathDatabase cath = new CathInstallation(config.getPdbFilePath());

		List<CathDomain> domains = cath.getDomainsForPdb(pdbID);

		try {

			// show the structure in 3D
			BiojavaJmol jmol = new BiojavaJmol();
			jmol.setStructure(StructureIO.getStructure(pdbID));
			jmol.evalString(DEFAULT_SCRIPT);

			System.out.println("got " + domains.size() + " domains");

			// now color the domains on the structure
			int colorpos = -1;

			for ( CathDomain domain : domains){

				colorpos++;

				showDomain(jmol, domain,colorpos);
			}


		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}



	private static void showDomain(BiojavaJmol jmol, CathDomain domain, int colorpos) {
		List<CathSegment> segments = domain.getSegments();

		StructureName key = new StructureName(domain.getDomainName());
		String chainId = key.getChainId();

		String color = colors[colorpos];

		System.out.println(" * domain " + domain.getDomainName() + " has # segments: " + domain.getSegments().size() + " color: " + color);

		for ( CathSegment segment : segments){
			System.out.println("   * " + segment);
			String start = segment.getStart();

			String stop = segment.getStop();

			String script = "select " + start + "-" + stop+":"+chainId + "; color " + color +";";

			jmol.evalString(script );
		}

	}

}
