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
package demo;


import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.cath.CathDatabase;
import org.biojava.nbio.structure.cath.CathDomain;
import org.biojava.nbio.structure.cath.CathInstallation;

/** An example for how to access CATH data.
 *
 * run with -Xmx512M
 *
 * @author Andreas Prlic
 *
 */
public class DemoCATH {
	public static void main(String[] args){

		AtomCache cache  = new AtomCache();

		CathDatabase database = new CathInstallation(cache.getPath());

		String domainID = "1cdgA01";

		CathDomain cathDomain = database.getDescriptionByCathId(domainID);

		System.out.println(cathDomain);


		Structure cathDomainStructure;
		try {
			cathDomainStructure = cache.getStructure(domainID);
			System.out.println(cathDomainStructure);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}




	}
}
