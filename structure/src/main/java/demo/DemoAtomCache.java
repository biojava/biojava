package demo;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.io.FileParsingParameters;

public class DemoAtomCache {
	public static void main(String[] args){
		
		boolean splitFileOrganisation = true;
		AtomCache cache = new AtomCache("/tmp",splitFileOrganisation);
		
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
