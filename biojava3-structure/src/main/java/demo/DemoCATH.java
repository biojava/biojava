package demo;


import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.cath.CathDatabase;
import org.biojava.bio.structure.cath.CathDomain;
import org.biojava.bio.structure.cath.CathInstallation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/** An example for how to access CATH data.
 * 
 * run with -Xmx512M
 * 
 * @author Andreas Prlic
 *
 */
public class DemoCATH {
	
	private static final Logger logger = LoggerFactory.getLogger(DemoCATH.class);

	public static void main(String[] args){
				
		AtomCache cache  = new AtomCache();
				
		CathDatabase database = new CathInstallation(cache.getPath());
		
		String domainID = "1cdgA01";
		
		CathDomain cathDomain = database.getDescriptionByCathId(domainID);
			
		logger.info("{}", cathDomain);
		
		Structure cathDomainStructure;
		try {
			cathDomainStructure = cache.getStructure(domainID);
			logger.info("{}", cathDomainStructure);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			logger.error("Exception: ", e);
		}		
	}
}