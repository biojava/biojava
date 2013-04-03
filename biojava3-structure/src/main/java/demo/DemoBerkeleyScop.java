package demo;

import org.biojava.bio.structure.scop.BerkeleyScopInstallation;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopFactory;

/** A demo for how to use the Berkeley version of SCOP instead of the default UK-SCOP
 * 
 * @author andreas
 * @since 3.0.4
 *
 */
public class DemoBerkeleyScop {
	public static void main(String[]args){


		ScopDatabase berkeley = new BerkeleyScopInstallation();

		ScopFactory.setScopDatabase(berkeley);

		// whenever you want to get access to SCOP now request it like this:
		ScopDatabase scop = ScopFactory.getSCOP();
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
