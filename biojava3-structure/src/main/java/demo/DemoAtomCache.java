package demo;


import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava3.structure.StructureIO;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/** Example of how to load PDB files using the AtomCache class.
 * 
 * @author Andreas Prlic
 *
 */
public class DemoAtomCache {
	
	private static final Logger logger = LoggerFactory.getLogger(DemoAtomCache.class);

	public static void main(String[] args){
		demoAtomCache();
		demoStructureIO();

	}

	@SuppressWarnings("unused")
	private static void demoStructureIO()  {


		try {
			Structure s1 = StructureIO.getStructure("4hhb");

			Structure bioAssembly = StructureIO.getBiologicalAssembly("1stp",1);
			
			// do something with them...
		} catch (Exception e){
			logger.error("Exception: ", e);
		}

	}

	private static void demoAtomCache() {
		AtomCache cache = new AtomCache();

		FileParsingParameters params = cache.getFileParsingParams();

		params.setLoadChemCompInfo(false);
		params.setAlignSeqRes(true);
		params.setHeaderOnly(false);
		params.setParseCAOnly(false);
		params.setParseSecStruc(false);

		String[] pdbIDs = new String[]{"4hhb", "1cdg","5pti","1gav", "WRONGID" };

		for (String pdbID : pdbIDs){

			try {
				Structure s = cache.getStructure(pdbID);
				if ( s == null) {
					logger.warn("could not find structure {}", pdbID);
					continue;
				}
				// do something with the structure
				logger.info("do something with: {}", s);

			} catch (Exception e){
				// something crazy happened...
				logger.error("Can't load structure {} reason: ", pdbID, e);
				//e.printStackTrace();
			}
		}
	}
}