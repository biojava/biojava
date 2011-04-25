/**
 * 
 */
package org.biojava.bio.structure;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.biojava.bio.structure.PDBStatus;
import org.biojava.bio.structure.PDBStatus.Status;

import junit.framework.TestCase;

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
		//TODO should test a case with a longer/more complicated version tree
		assertEquals("4HHB",PDBStatus.getReplacement("1HHB",false));
		assertEquals("3HHB",PDBStatus.getReplacement("3HHB",false));
	}
	

	public void testGetCurrent() {
		//TODO should test a case with a longer/more complicated version tree
		assertEquals("4HHB",PDBStatus.getCurrent("1HHB"));
		assertEquals("3HHB",PDBStatus.getCurrent("3HHB"));
		assertEquals(null, PDBStatus.getCurrent("1CMW"));
		assertEquals("3ENI",PDBStatus.getCurrent("1KSA"));
	}
	
	public void testGetReplaces() {
		//TODO should test a case with a longer/more complicated version tree
		assertEquals(new ArrayList<String>(), Arrays.asList(new String[] {}));
		
		assertEquals(Arrays.asList("1HHB"),PDBStatus.getReplaces("4HHB"));
		assertEquals(Arrays.asList("1HHB"),PDBStatus.getReplaces("3HHB"));
		assertEquals(Arrays.asList(), PDBStatus.getReplaces("1HHB"));
		assertEquals(Arrays.asList("1M50","1KSA"),PDBStatus.getReplaces("3ENI",false));
		assertEquals(Arrays.asList("1M50","1KSA"),PDBStatus.getReplaces("3ENI",true));
	}
	
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
			assertEquals("Wrong replacedBy","4HHB",attrs.get("replacedBy"));
			
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
			assertEquals("Wrong replacedBy","4HHB",attrs.get("replacedBy"));
			
		} catch(Exception e) {
			e.printStackTrace();
			fail();
		}
	}
	
}
