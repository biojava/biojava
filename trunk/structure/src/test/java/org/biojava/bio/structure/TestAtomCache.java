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

package org.biojava.bio.structure;

import org.biojava.bio.structure.align.util.AtomCache;

import junit.framework.TestCase;

public class TestAtomCache extends TestCase
{


   public void testAtomCacheNameParsing(){

      String name1= "4hhb";
      
      String name2 = "4hhb.C";
      String chainId2 = "C";
      
      String name3 = "4hhb:1";
      String chainId3 = "B";
      
      String name4 = "4hhb:A:10-20,B:10-20,C:10-20";
      String name5 = "4hhb:(A:10-20,A:30-40)";
      
      String tmpDir = System.getProperty("java.io.tmpdir");
      
      boolean isSplit = true;
      AtomCache cache = new AtomCache(tmpDir,isSplit);

      try {
         Structure s = cache.getStructure(name1);
         assertNotNull(s);
         assertTrue(s.getChains().size() == 4);
         s = cache.getStructure(name2);
         
         assertTrue(s.getChains().size() == 1);
         Chain c = s.getChainByPDB(chainId2);
         assertEquals(c.getName(),chainId2);
         
         s = cache.getStructure(name3);
         assertNotNull(s);
         assertTrue(s.getChains().size() == 1);
         
         c = s.getChainByPDB(chainId3);
         assertEquals(c.getName(),chainId3);
         
         s = cache.getStructure(name4);
         assertNotNull(s);
  
         assertEquals(s.getChains().size(), 3);
         
         c = s.getChainByPDB("B");
         assertEquals(c.getAtomLength(),11);
         
         s =cache.getStructure(name5);
         assertNotNull(s);
         
         assertEquals(s.getChains().size(),1 );
         c = s.getChainByPDB("A");
         assertEquals(c.getAtomLength(),22);
         
         
      } catch (Exception e){
         e.printStackTrace();
         fail(e.getMessage());
      }

   }
}
