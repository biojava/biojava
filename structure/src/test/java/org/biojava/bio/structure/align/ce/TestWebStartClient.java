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
 * Created on May 18, 2010
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.ce;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.TmpAtomCache;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.client.JFatCatClient;

import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.xml.AFPChainFlipper;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.align.xml.AFPChainXMLParser;

import junit.framework.TestCase;


public class TestWebStartClient extends TestCase
{

   public void testCPAlignment(){
      
      String name1="1cdg.A";
      String name2="1tim.A";
      
      
      
      AtomCache cache = TmpAtomCache.cache;
      try {
        
         Atom[] ca1 = cache.getAtoms(name1);
         Atom[] ca2 = cache.getAtoms(name2);
         
         String serverLocation = "http://beta.rcsb.org/pdb/rest/";
         AFPChain afpServer = JFatCatClient.getAFPChainFromServer(serverLocation,CeCPMain.algorithmName, name1, name2, ca1, ca2, 5000);
         assertNotNull(afpServer); 
         
         assertTrue(afpServer.getAlgorithmName().equals(CeCPMain.algorithmName));
         assertTrue(afpServer.getBlockNum() > 1);
                 
         StructureAlignment alignment = StructureAlignmentFactory.getAlgorithm(CeCPMain.algorithmName);
        
         
         AFPChain afpChain = alignment.align(ca1,ca2);
         afpChain.setName1(name1);  
         afpChain.setName2(name2);
         
         assertNotNull(afpChain);
         assertNotNull(afpChain.getAlgorithmName());
         assertTrue(afpChain.getAlgorithmName().equals(CeCPMain.algorithmName));
         
         String xml = AFPChainXMLConverter.toXML(afpChain,ca1,ca2);
                
         String xml2 = AFPChainXMLConverter.toXML(afpServer, ca1, ca2);
         assertEquals(xml,xml2);
         
         AFPChain afpFlip = AFPChainFlipper.flipChain(afpChain);
         String xml3 = AFPChainXMLConverter.toXML(afpFlip, ca2, ca1);
         AFPChain afpDoubleFlip = AFPChainXMLParser.fromXML(xml3, ca2, ca1);
         String xml5 = AFPChainXMLConverter.toXML(afpDoubleFlip, ca1, ca2);
         
         assertEquals(xml,xml5);
         
         
         
      } catch (Exception e){
         fail (e.getMessage());
      }
   }
}
