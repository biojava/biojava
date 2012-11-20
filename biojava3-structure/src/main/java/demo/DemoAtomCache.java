package demo;


import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava3.structure.StructureIO;


/** Example of how to load PDB files using the AtomCache class.
 * 
 * @author Andreas Prlic
 *
 */
public class DemoAtomCache {
	public static void main(String[] args){
		demoAtomCache();
		demoStructureIO();

	}

	private static void demoStructureIO()  {


		try {
			Structure s1 = StructureIO.getStructure("4hhb");

			Structure bioAssembly = StructureIO.getBiologicalAssembly("1stp",1);
			
			// do something with them...
		} catch (Exception e){
			e.printStackTrace();
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
					System.out.println("could not find structure " + pdbID);
					continue;
				}
				// do something with the structure
				System.out.println(s);

			} catch (Exception e){
				// something crazy happened...
				System.err.println("Can't load structure " + pdbID + " reason: " + e.getMessage());
				e.printStackTrace();
			}
		}

	}
}
