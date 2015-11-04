package org.biojava.nbio.structure.io;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;
import org.junit.Assert;
import org.junit.Test;

/**
 * Created by edlunde-dnastar 
 * @since 10/30/2015.
 */
public class TestParseMmCIFLigands {
	
	int HEM_COUNT_4HHB = 172;	//Number of atoms in HEM groups of 4HHB (manually determined from CIF file)
	
	@Test
	public void testLigandConnections()throws IOException, StructureException {
		AtomCache cache = new AtomCache();

		StructureIO.setAtomCache(cache);

		cache.setUseMmCif(true);
		ChemCompGroupFactory.setChemCompProvider(new DownloadChemCompProvider());

		FileParsingParameters params = cache.getFileParsingParams();
		params.setCreateLigandConects(true);
		StructureIO.setAtomCache(cache);

		Structure sCif = StructureIO.getStructure("4HHB");
		List<Map<String, Integer>> conects = sCif.getConnections();
		
		assertNotNull(conects);
		Assert.assertFalse(conects.isEmpty());
		//Verify that we have all HEM atoms from the CIF file.
		Assert.assertTrue(countUniqueAtomsInConnectionsList(conects) == HEM_COUNT_4HHB);
	}
	
	private int countUniqueAtomsInConnectionsList(List<Map<String, Integer>> conects){

		HashSet<Integer> uniqueIDs = new HashSet<Integer>();
		for ( Map<String, Integer> c : conects){
			uniqueIDs.add(c.get("atomserial"));
			for (int i = 1; i < 10; i ++ ){	//Look for each bond. Max bound for 4hhb/Fe is 4, but I'll err on the side of caution.
				if ( c.get("bond" + i) != null ){
					uniqueIDs.add(c.get("bond"+ i));
				}else{
					break;
				}
			}
		}
		return uniqueIDs.size();
	}
}