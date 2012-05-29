/*
 *                    PDB web development code
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
 *
 * Created on Jul 23, 2009
 * Created by ap3
 *
 */

package org.biojava.bio.structure.align.fatcat;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.fatcat.calc.FatCatParameters;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.align.xml.AFPChainXMLParser;
import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileReader;



public class MyTestHelper
{

   // 
   public static final String pdbPath = System.getProperty("java.io.tmpdir");

   public static String compareAlignment(String pdb1, String chain1, String pdb2, String chain2, String originalOutput, boolean doRigid){

      AFPChain afpChain = null;

      PDBFileReader pdbpars = new PDBFileReader();
      pdbpars.setPath(pdbPath);
      pdbpars.setAutoFetch(true);
      
      FileParsingParameters params = new FileParsingParameters();
      params.setAlignSeqRes(true);
      params.setLoadChemCompInfo(false);
      params.setParseCAOnly(true);
      pdbpars.setFileParsingParameters(params);
      
      Structure structure1;
      Structure structure2;

      //                   default:      new:
      // 1buz - 1ali : time: 8.3s eqr 68 rmsd 3.1 score 161 | time 6.4 eqr 58 rmsd 3.0 scre 168 | rigid: identical, flexible: not significant alignment, 
      // 5pti - 1tap : time: 6.2s eqr 48 rmsd 2.67 score 164 | time 5.2 eqr 49 rmsd 2.9 score 151 | rigid: 
      // 1cdg - 8tim
      // 1jbe - 1ord : identical with fatcat
      // 1nbw.A - 1kid : rigid: identical, flexible: not identical, alignment not significant.
      // 1t4y - 1rp5
      // 1a64.A - 1hng.B
      // 1zzw - 1bw6

      try {
         structure1 = pdbpars.getStructureById(pdb1);
         structure2 = pdbpars.getStructureById(pdb2);

         //structure1 = pdbpars.getStructureById("1cdg");
         Chain c1 = structure1.getChainByPDB(chain1);

         //structure2 = pdbpars.getStructureById("2aaa");         
         Chain c2 = structure2.getChainByPDB(chain2);

         Structure s3 = new StructureImpl();
         s3.addChain(c1);

         Structure s4 = new StructureImpl();
         s4.addChain(c2);

         Atom[] ca1 = StructureTools.getAtomCAArray(s3);
         Atom[] ca2 = StructureTools.getAtomCAArray(s4);

         // keep an independent copy of them for tests further down..
         Atom[] ca3 = new Atom[ca2.length];
         for (int i = 0 ; i < ca2.length; i++){
            Group g = (Group)ca2[i].getGroup().clone();
            g.setChain(ca2[i].getGroup().getChain());
            ca3[i] = g.getAtom(StructureTools.caAtomName);
         }

         StructureAlignment fatCat;

         FatCatParameters fparams = new FatCatParameters();

         if ( doRigid)
            fatCat = new FatCatRigid();            
         else 
            fatCat = new FatCatFlexible();

         afpChain = fatCat.align(ca1, ca2, fparams);

         afpChain.setName1(pdb1+chain1);
         afpChain.setName2(pdb2+chain2);


         // TEST THE XML SERIALIZATION AND DE_SERIALIZATION!
         String result = afpChain.toFatcat(ca1, ca2);

         String xml = AFPChainXMLConverter.toXML(afpChain,ca1,ca2);

         //System.out.println(xml);

         AFPChain newChain = AFPChainXMLParser.fromXML (xml, ca1, ca3);

         String resultSerialized = newChain.toFatcat(ca1, ca3);
         //System.out.println("*** RESULT2 "+result2);

         if ( ! result.equals(resultSerialized)) {
            System.out.println("not identical toFatCat()!!!");
            System.out.println(xml);
            System.out.println(result);
            System.out.println("***");
            System.out.println(resultSerialized);
            throw new StructureException("the JFatCat alignment image does not look identical after XML serialization.");        	 
         }

         if ( ! afpChain.toString().equals(newChain.toString())){
            System.err.println("not identical toStrings!!!");
            System.err.println(afpChain.toString());
            System.err.println(newChain.toString());
            throw new StructureException("The AFPChain.toString() does not look identical after XML serialization.");
         }



         if ( result.equals(originalOutput)){
            return "";
         } else {
            return result;
         }

      } catch ( Exception e){
         e.printStackTrace();
         return e.getMessage();
      }

   }

}
