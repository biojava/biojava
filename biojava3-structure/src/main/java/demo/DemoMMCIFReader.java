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
 * Created on May 17, 2010
 * Author: Andreas Prlic 
 *
 */

package demo;

import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.GroupType;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.io.MMCIFFileReader;
import org.biojava.bio.structure.io.StructureIOFile;

/** An example of how to read MMcif files
 * 
 * @author Andreas Prlic
 * 
 */
public class DemoMMCIFReader
{

   public static void main(String[] args){
      String pdbId = "193D";
      
      StructureIOFile pdbreader = new MMCIFFileReader();
      
      try {
         pdbreader.setAutoFetch(true);
          Structure s = pdbreader.getStructureById(pdbId);
          System.out.println(s);
       
          Chain c = s.getChainByPDB("A");

          System.out.println(c.getSeqResSequence());
          System.out.println(c.getAtomSequence());
          System.out.println(c.getAtomGroups(GroupType.HETATM));
          Chain d = s.getChainByPDB("B");
          System.out.println(d.getSeqResSequence());
          System.out.println(d.getAtomSequence());
        
          
          for (Group g : d.getAtomGroups(GroupType.HETATM)){
             System.out.println(g.getResidueNumber() +  " " +  g.getPDBName() + " " + g);
          }
      } catch (Exception e) {
          e.printStackTrace();
      }


   }
}
