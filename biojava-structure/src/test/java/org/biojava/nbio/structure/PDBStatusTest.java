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
/**
 * 
 */
package org.biojava.nbio.structure;

import junit.framework.TestCase;
import org.biojava.nbio.structure.PDBStatus.Status;

import java.lang.reflect.Method;
import java.util.*;

/**
 * @author Spencer Bliven <sbliven@ucsd.edu>
 *
 */
public class PDBStatusTest extends TestCase {

	/**
	 * Test {@link PDBStatus#getStatus(String)}.
	 * 
	 * <p>Uses the following PDBs:<br/>
	 * <pre>1HHB    OBSOLETE	replacedBy=4HHB
	 *4HHB    CURRENT	replaces=1HHB
	 *3HHB    CURRENT	replaces=1HHB
	 *</pre>
	 */
	public void testGetStatus() {
		assertEquals(Status.OBSOLETE, PDBStatus.getStatus("1HHB"));
		assertEquals(Status.CURRENT, PDBStatus.getStatus("3HHB"));
		assertEquals(Status.CURRENT, PDBStatus.getStatus("4HHB"));
	}
	
	public void testGetReplacement() {
		assertFalse(Arrays.asList("YES").equals(Arrays.asList("NO"))); //check for deep equals
		
		// 1CMW is replacedBy NONE
		assertEquals(Arrays.asList(),PDBStatus.getReplacement("1CMW", true, false));
		assertEquals(Arrays.asList("1CMW"),PDBStatus.getReplacement("1CMW", true, true));
		
		// 1HHB is replacedBy 2-4HHB
		assertEquals(Arrays.asList("3HHB"),PDBStatus.getReplacement("3HHB",false,false));
		assertEquals(Arrays.asList("3HHB"),PDBStatus.getReplacement("3HHB",false,true));
		assertEquals(Arrays.asList("4HHB","3HHB","2HHB"),PDBStatus.getReplacement("1HHB",false,false));
		assertEquals(Arrays.asList("4HHB","3HHB","2HHB","1HHB"),PDBStatus.getReplacement("1HHB",false,true));
		
		// 1CAT is replacedBy 3CAT is replacedBy 7-8CAT
		assertEquals(Arrays.asList("8CAT","7CAT","3CAT","1CAT"),PDBStatus.getReplacement("1CAT",true,true));
		assertEquals(Arrays.asList("8CAT","7CAT"),PDBStatus.getReplacement("1CAT",true,false));
		assertEquals(Arrays.asList("8CAT","7CAT","3CAT"),PDBStatus.getReplacement("3CAT",true,true));
		assertEquals(Arrays.asList("8CAT","7CAT"),PDBStatus.getReplacement("3CAT",true,false));
	}
	

	public void testGetCurrent() {
		assertEquals("4HHB",PDBStatus.getCurrent("1HHB"));
		assertEquals("3HHB",PDBStatus.getCurrent("3HHB"));
		assertEquals(null, PDBStatus.getCurrent("1CMW"));
		assertEquals("3ENI",PDBStatus.getCurrent("1KSA"));
		assertEquals("8CAT",PDBStatus.getCurrent("1CAT"));
		assertEquals("8CAT",PDBStatus.getCurrent("3CAT"));
		assertEquals("7CAT",PDBStatus.getCurrent("7CAT"));
	}
	
	public void testGetReplaces() {
		assertEquals(new ArrayList<String>(), Arrays.asList(new String[] {}));
		
		assertEquals(Arrays.asList("1HHB"),PDBStatus.getReplaces("4HHB",false));
		assertEquals(Arrays.asList("1HHB"),PDBStatus.getReplaces("3HHB",false));
		assertEquals(Arrays.asList(), PDBStatus.getReplaces("1HHB", false));
		assertEquals(Arrays.asList("1M50","1KSA"),PDBStatus.getReplaces("3ENI",false));
		assertEquals(Arrays.asList("1M50","1KSA"),PDBStatus.getReplaces("3ENI",true));
		assertEquals(Arrays.asList("3CAT"),PDBStatus.getReplaces("8CAT",false));
		assertEquals(Arrays.asList("3CAT","1CAT"),PDBStatus.getReplaces("8CAT",true));
		
	}
	
	/**
	 * Tests a helper method for merging that was giving me problems
	 */
	public void testMergeReversed() {
		try {
			Method mergeReversed = PDBStatus.class.getDeclaredMethod("mergeReversed",
					List.class,List.class);
			mergeReversed.setAccessible(true);


			List<String> a,b;

			b = Arrays.asList("F","A");
			a = new LinkedList<String>();
			mergeReversed.invoke(null, a,b);
			assertEquals(Arrays.asList("F","A"),a);

			a = new LinkedList<String>();
			a.add("B");
			mergeReversed.invoke(null, a,b);
			assertEquals(Arrays.asList("F","B","A"),a);

			a = new LinkedList<String>();
			a.add("G");
			mergeReversed.invoke(null, a,b);
			assertEquals(Arrays.asList("G","F","A"),a);

			a = new LinkedList<String>();
			a.add("1");
			mergeReversed.invoke(null, a,b);
			assertEquals(Arrays.asList("F","A", "1"),a);

			a = new LinkedList<String>();
			a.add("G");
			a.add("1");
			mergeReversed.invoke(null, a,b);
			assertEquals(Arrays.asList("G","F","A", "1"),a);
			
			b = Arrays.asList();
			mergeReversed.invoke(null, a,b);
			assertEquals(Arrays.asList("G","F","A", "1"),a);

			b = Arrays.asList("G","D","C","A");
			a = new LinkedList<String>();
			a.add("F");
			a.add("B");
			a.add("1");
			mergeReversed.invoke(null, a,b);
			assertEquals(Arrays.asList("G","F","D","C","B","A", "1"),a);

		} catch(Exception e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
	
	/**
	 * Test low-level connectivity to the PDB
	 */
	@SuppressWarnings("unchecked")
	public void testGetStatusIdRecords() {
		try {
		Method getStatusIdRecords = PDBStatus.class.getDeclaredMethod("getStatusIdRecords",
				String[].class);
		getStatusIdRecords.setAccessible(true);

		
			List<Map<String,String>> attrsList;
			String[] pdbIds;
			Map<String,String> attrs;
			
			// Test invocation with a single ID
			pdbIds = new String[] {"1HHB"};
			attrsList = (List<Map<String,String>>) getStatusIdRecords.invoke(null, (Object) pdbIds);
			assertEquals("Wrong number of records.",1, attrsList.size());
			attrs = attrsList.get(0);
			assertEquals("Wrong number of attributes",3,attrs.size());
			assertEquals("Wrong structureId","1HHB",attrs.get("structureId"));
			assertEquals("Wrong status","OBSOLETE",attrs.get("status"));
			assertEquals("Wrong replacedBy","4HHB 3HHB 2HHB",attrs.get("replacedBy"));
			
			// Test with multiple IDs
			pdbIds = new String[] {"1HHB","4HHB"};
			attrsList = (List<Map<String,String>>) getStatusIdRecords.invoke(null, (Object) pdbIds);
			assertEquals("Wrong number of records.",2, attrsList.size());
			attrs = attrsList.get(1);
			assertEquals("Wrong number of attributes",3,attrs.size());
			assertEquals("Wrong structureId","4HHB",attrs.get("structureId"));
			assertEquals("Wrong status","CURRENT",attrs.get("status"));
			assertEquals("Wrong replaces","1HHB",attrs.get("replaces"));
			attrs = attrsList.get(0);
			assertEquals("Wrong number of attributes",3,attrs.size());
			assertEquals("Wrong structureId","1HHB",attrs.get("structureId"));
			assertEquals("Wrong status","OBSOLETE",attrs.get("status"));
			assertEquals("Wrong replacedBy","4HHB 3HHB 2HHB",attrs.get("replacedBy"));
			
			// Test invocation with a single ID
			pdbIds = new String[] {"3ENI"};
			attrsList = (List<Map<String,String>>) getStatusIdRecords.invoke(null, (Object) pdbIds);
			assertEquals("Wrong number of records.",1, attrsList.size());
			attrs = attrsList.get(0);
			assertEquals("Wrong number of attributes",3,attrs.size());
			assertEquals("Wrong structureId","3ENI",attrs.get("structureId"));
			assertEquals("Wrong status","CURRENT",attrs.get("status"));
			assertEquals("Wrong replacedBy","1M50 1KSA",attrs.get("replaces"));
			
			
		} catch(Exception e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}
	
}
