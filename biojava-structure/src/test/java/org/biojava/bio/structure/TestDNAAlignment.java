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
 * Created on May 11, 2010
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure;

import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;

import junit.framework.TestCase;

/** make sure DNA alignments fail gracefully
 * 
 * @author Andreas Prlic
 *
 */
public class TestDNAAlignment extends TestCase
{

   public void test1(){
      String name1="1l3s.A";
      String name2="1t7p.P";

      AtomCache cache = new AtomCache();
      try {
         Atom[] ca1 = cache.getAtoms(name1);
         Atom[] ca2 = cache.getAtoms(name2);
         CeMain ce = new CeMain();
         AFPChain afpChain = ce.align(ca1,ca2);
         assertNotNull(afpChain);
         
        String txt = afpChain.toFatcat(ca1, ca2);
        
        assertNotNull(txt);
        
      } catch (Exception e){
         e.printStackTrace();
         fail(e.getMessage());
      }
   }
}
