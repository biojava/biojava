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

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.scop.ScopDatabase;
import org.biojava.nbio.structure.scop.ScopFactory;

/** A demo for how to use the Berkeley version of SCOP instead of the default UK-SCOP
 *
 * @author andreas
 * @since 3.0.4
 *
 */
public class DemoBerkeleyScop {

	public static void main(String[]args){


		//ScopDatabase berkeley = new BerkeleyScopInstallation();

		//ScopFactory.setScopDatabase(berkeley);

		AtomCache cache = new AtomCache();
		// whenever you want to get access to SCOP now request it like this:
		ScopDatabase scop = ScopFactory.getSCOP("1.75");
		ScopFactory.setScopDatabase(scop);


		System.out.println(cache.getPath());
		System.out.println(cache.getCachePath());
		// ... and do something with it


		// eg. you can run all the demos that work for the UK - SCOP (currently at version 1.75)
		// this demo no automatically picks up the Berkeley version (currently 1.75A)
		DemoSCOP scopDemo = new DemoSCOP();

		scopDemo.getCategories();
		scopDemo.printDomainsForPDB();
		scopDemo.traverseHierarchy();
		scopDemo.alignSuperfamily();


	}
}
