package demo;

import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopDatabase;
import org.biojava.bio.structure.scop.ScopFactory;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/** A demo for how to use the Berkeley version of SCOP instead of the default UK-SCOP
 * 
 * @author andreas
 * @since 3.0.4
 *
 */
public class DemoBerkeleyScop {
	
	private static final Logger logger = LoggerFactory.getLogger(DemoBerkeleyScop.class);

	public static void main(String[]args){


		//ScopDatabase berkeley = new BerkeleyScopInstallation();

		//ScopFactory.setScopDatabase(berkeley);

		AtomCache cache = new AtomCache();
		// whenever you want to get access to SCOP now request it like this:
		ScopDatabase scop = ScopFactory.getSCOP("1.75");
		ScopFactory.setScopDatabase(scop);
		
		
		logger.info("Atom Cache Path: {}", cache.getPath());
		logger.info("{}", cache.getCachePath());
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