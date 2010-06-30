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
 * Created on Jun 30, 2010
 * Author: ap3 
 *
 */

package demo;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeMain;

import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.scop.ScopCategory;
import org.biojava.bio.structure.scop.ScopDescription;
import org.biojava.bio.structure.scop.ScopDomain;
import org.biojava.bio.structure.scop.ScopInstallation;
import org.biojava.bio.structure.scop.ScopNode;

public class DemoSCOP
{
   public static void main(String[] args){

      String cacheLocation = "/tmp/";

      ScopInstallation scop = new ScopInstallation(cacheLocation);

      List<ScopDomain> domains = scop.getDomainsForPDB("4HHB");

      System.out.println(domains);

      List<ScopDescription> superfams = scop.getByCategory(ScopCategory.Superfamily);

      System.out.println("Total nr. of superfamilies:" + superfams.size());


      ScopDescription superfam1 = superfams.get(0);
      List<ScopDomain> doms4superfam1 = scop.getScopDomainsBySunid(superfam1.getSunID());
      ScopDomain dom1 = doms4superfam1.get(0);
      
      AtomCache cache = new AtomCache(cacheLocation, true);
      
      for ( int i = 1 ; i < superfams.size() ; i ++){
         ScopDescription superfam = superfams.get(i);

         ScopNode node = scop.getScopNode(superfam.getSunID());
         
         List<ScopDomain> doms = scop.getScopDomainsBySunid(superfam.getSunID());

         ScopDomain dom2 = doms.get(0);
         align(dom1,dom2,cache);
      }

    

   }
   
   public static final void align(ScopDomain dom1, ScopDomain dom2, AtomCache cache){

      try {
         Structure s1 = cache.getStructureForDomain(dom1);
         Structure s2 = cache.getStructureForDomain(dom2);
         
         Atom[] ca1 = StructureTools.getAtomCAArray(s1);
         Atom[] ca2 = StructureTools.getAtomCAArray(s2);
         StructureAlignment ce = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
         AFPChain afpChain = ce.align(ca1, ca2);
         
         //System.out.println(afpChain.toCE(ca1, ca2));
         
         //StructureAlignmentDisplay.display(afpChain, ca1, ca2);
         
         System.out.println(dom1.getScopId() + " vs. " + dom2.getScopId()+ " :" + afpChain.getProbability());
      } catch (Exception e){
         e.printStackTrace();
      }

   }
}
