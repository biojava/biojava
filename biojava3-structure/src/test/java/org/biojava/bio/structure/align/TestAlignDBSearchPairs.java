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
 * Created on Mar 1, 2010
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align;



import junit.framework.TestCase;

public class TestAlignDBSearchPairs extends TestCase
{

	
	public void testNothing(){
		
	}
	
	
	// speedup... nothing new being tested here, so disabling for now
//   public void testParsePairs(){
//
//      String tmpDir = System.getProperty("java.io.tmpdir");
//      
//      AtomCache cache = new AtomCache(tmpDir,true);
//
//      InputStream inStream = this.getClass().getResourceAsStream("/db_search.pairs");
//      assertNotNull(inStream);
//
//      BufferedReader is = new BufferedReader (new InputStreamReader(inStream)) ;
//      try {
//         StructureAlignment algorithm =  StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
//         String line = null;
//         while ( (line = is.readLine()) != null){
//            if ( line.startsWith("#"))
//               continue;
//           // System.out.println("aligning: " + line);
//            String[] spl = line.split(" ");
//            String pdb1 = spl[0];
//            String pdb2 = spl[1];
//
//
//            Structure structure1 = cache.getStructure(pdb1);
//            Structure structure2 = cache.getStructure(pdb2);
//
//            Atom[] ca1;
//            Atom[] ca2;
//
//
//            ca1 = StructureTools.getAtomCAArray(structure1);
//            ca2 = StructureTools.getAtomCAArray(structure2);
//
//            algorithm.align(ca1,ca2);
//
//         }
//      } catch (Exception e){
//         e.printStackTrace();
//         fail(e.getMessage());
//      }
//   }
}
