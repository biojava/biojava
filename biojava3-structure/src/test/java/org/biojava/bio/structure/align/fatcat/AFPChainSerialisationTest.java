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

import java.io.IOException;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.fatcat.calc.FatCatParameters;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.align.xml.AFPChainXMLParser;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

public class AFPChainSerialisationTest 
extends TestCase
{



   /**
    * @return the suite of tests being tested
    */
   public static Test suite()
   {
      return new TestSuite( AFPChainSerialisationTest.class );
   }

   public AFPChainSerialisationTest(String testName){
      super(testName);
   }

//   public void testSerialization1P80_2IUF(){
//      String name1 = "1P80.D";
//      String name2 = "2IUF.E";
//
//      AtomCache cache =    TmpAtomCache.cache;
//      try {
//         
//         Atom[] ca1 = cache.getAtoms(name1);         
//         Atom[] ca2 = cache.getAtoms(name2);
//         
//         testAlignment(name1,name2 ,ca1,ca2,false);
//         testAlignment(name1,name2,ca1,ca2,true);
//         
//      } catch (Exception e){
//         e.printStackTrace();
//         fail(e.getMessage());
//      }
//
//   }

   public void testSerialization1a21_1hwg(){


      try {

         Structure s1 = getStructure("1a21", "A");			
         Structure s2 = getStructure("1hwg","C");


         Atom[] ca1 = StructureTools.getAtomCAArray(s1);
         Atom[] ca2 = StructureTools.getAtomCAArray(s2);

         testAlignment("1a21.A","1hwg.C",ca1,ca2,false);
         testAlignment("1a21.A","1hwg.C",ca1,ca2,true);

      } catch (Exception e) {
         e.printStackTrace();
         fail(e.getMessage());
      }
   }

   
   // slow...
//   public void testSerialization1cdg_1tim(){
//
//      try {
//
//         Structure s1 = getStructure("1cdg", "A");		
//         Structure s2 = getStructure("1tim","A");
//
//         Atom[] ca1 = StructureTools.getAtomCAArray(s1);
//         Atom[] ca2 = StructureTools.getAtomCAArray(s2);
//
//         testAlignment("1cdg.A","1tim.A",ca1,ca2,false);
//         testAlignment("1cdg.A","1tim.A",ca1,ca2,true);
//
//      } catch (Exception e) {
//         e.printStackTrace();
//         fail(e.getMessage());
//      }
//   }
//   public void testSerialization4hhb(){
//
//      try {
//
//         Structure s1 = getStructure("4hhb", "A");		
//         Structure s2 = getStructure("4hhb", "B");
//
//         Atom[] ca1 = StructureTools.getAtomCAArray(s1);
//         Atom[] ca2 = StructureTools.getAtomCAArray(s2);
//
//         testAlignment("4hhb.A","4hhb.B",ca1,ca2,false);
//         testAlignment("4hhb.A","4hhb.B",ca1,ca2,true);
//
//      } catch (Exception e) {
//         e.printStackTrace();
//         fail(e.getMessage());
//      }
//   }


   private Structure getStructure(String pdbId, String chainId) throws IOException, StructureException{
	   AtomCache cache = new AtomCache();
      Structure structure1 = cache.getStructure(pdbId);

      Chain c = structure1.getChainByPDB(chainId);

      Structure s = new StructureImpl();
      s.addChain(c);

      return s;

   }

   private void testAlignment(String name1, String name2, Atom[] ca1, Atom[] ca2, boolean doRigid) throws StructureException,IOException{


      Atom[] ca3 = StructureTools.cloneCAArray(ca2);


      AFPChain afpChain = doAlign(name1, name2, ca1,ca2,doRigid);



      String fatcat = afpChain.toFatcat(ca1, ca2);
      String xml 	  = AFPChainXMLConverter.toXML(afpChain,ca1,ca2);


      //System.out.println(xml);
      AFPChain newChain = AFPChainXMLParser.fromXML (xml, ca1, ca3);

      // test blockNum and optLen arrays 
      int blockNum = afpChain.getBlockNum();
      int[] optLen = afpChain.getOptLen();

      assertTrue("The nr of aligned blocks is not the same! " + blockNum + " " + newChain.getBlockNum() , blockNum == newChain.getBlockNum());



      for ( int i =0 ; i < blockNum ; i++){
         int newLenI = newChain.getOptLen()[i];
         assertTrue("The values in the optLen field don't match! pos:" + i + " orig:" + optLen[i] + " new:" +  newLenI,optLen[i] == newLenI);
      }

      // test the internal optAlign data structure:

      int[][][] optAln1 = afpChain.getOptAln();
      int[][][] optAln2 = newChain.getOptAln();

      for(int i = 0; i < blockNum; i ++)  {        	
         for(int j = 0; j < optLen[i]; j ++) {
            int p1 = optAln1[i][0][j];
            int p2 = optAln1[i][1][j];

            int n1 = optAln2[i][0][j];
            int n2 = optAln2[i][1][j];

            assertTrue(p1 == n1);
            assertTrue(p2 == n2);

         }
      }

      String fatcat2 = newChain.toFatcat(ca1, ca3);
      //System.out.println("*** RESULT2 "+result2);

      assertEquals(fatcat,fatcat2);

      String xml2 = AFPChainXMLConverter.toXML(newChain, ca1, ca3);
      assertEquals(xml, xml2);

   }


   public AFPChain doAlign(String name1, String name2, Atom[] ca1, Atom[] ca2 , boolean doRigid) throws StructureException,IOException{
      FatCatParameters params = new FatCatParameters();

      FatCat fatCat = new FatCat();
      AFPChain afpChain;
      if ( doRigid)
         afpChain = fatCat.alignRigid(ca1,ca2,params);
      else
         afpChain = fatCat.alignFlexible(ca1,ca2,params);

      afpChain.setName1(name1);
      afpChain.setName2(name2);
      return afpChain;

   }

   public String[] align (String name1, String name2, Atom[] ca1, Atom[] ca2 , boolean doRigid) throws StructureException,IOException{

      AFPChain afpChain = doAlign(name1, name2, ca1, ca2, doRigid);
      // flexible original results:
      String fatcat = afpChain.toFatcat(ca1,ca2);
      //System.out.println(result1);


      String xml = AFPChainXMLConverter.toXML(afpChain,ca1,ca2);

      String[] result = new String[2];
      result[0] = fatcat;
      result[1] = xml;
      return result;
   }


   public void testMulti() {
      Atom[] ca1 = null;
      Atom[] ca2 = null;
      Atom[] ca3 = null;
      Atom[] ca4 = null;
      Atom[] ca5 = null;
      Atom[] ca6 = null;
      String[] result1 = null;
      String[] result2 = null;

      String name1 = "5pti.A";
      String name2 = "1znf.A";

      String name3 ="1hiv.A";
      String name4 ="1a4w.H";
      try {
         Structure s1 = getStructure("5pti","A");
         Structure s2 = getStructure("1znf","A");
         ca1 = StructureTools.getAtomCAArray(s1);
         ca2 = StructureTools.getAtomCAArray(s2);
         ca3 = StructureTools.cloneCAArray(ca2);

         result1 = align(name1,name2,ca1, ca2,true);

         Structure s3 = getStructure("1hiv","A");
         Structure s4 = getStructure("1a4w","H");
         ca4 = StructureTools.getAtomCAArray(s3);
         ca5 = StructureTools.getAtomCAArray(s4);
         ca6 = StructureTools.cloneCAArray(ca5);

         result2 = align(name3,name4,ca4, ca5,true);


      } catch (Exception e){
         e.printStackTrace();
         fail(e.getMessage());
      }

      String xmlNew = "<multi>"+result1[1]+ result2[1] +"</multi>";
      //System.out.println(xmlNew);

      testSaxParser(xmlNew,name1,name2,result1,ca1,ca3,name3,name4,result2,ca4,ca6);


      //WARNING: THE ORDER CAN CHANGE: order of elements in XML is not necessarily the same!
      AFPChain[] chains = AFPChainXMLParser.parseMultiXML(xmlNew);

      assertTrue(chains.length == 2);

      // recreate the correct chains...
      AFPChain new1 = getAfpFor(name1,chains);
      AFPChain new2 = getAfpFor(name3,chains);

      assertNotNull(new1);
      assertNotNull(new2);

      assertTrue(new1.getName1().equals(name1));
      //System.out.println(new2.getName1() + " " + new2.getName2() + " "+ name3);
      assertTrue(new2.getName1().equals(name3));

      try {
         AFPChainXMLParser.rebuildAFPChain(new1, ca1, ca3);				
         String fatcat1 = new1.toFatcat(ca1, ca3); 
         assertEquals(fatcat1, result1[0]);				
         String xmlnew1 = AFPChainXMLConverter.toXML(new1, ca1, ca3);
         assertTrue(xmlnew1.equals(result1[1]));

         AFPChainXMLParser.rebuildAFPChain(new2, ca4, ca6);
         String fatcat2 = new2.toFatcat(ca4, ca6);
        
         assertEquals(fatcat2,result2[0]);			
         String xmlnew2 = AFPChainXMLConverter.toXML(new2, ca4, ca6);
         assertEquals(xmlnew2,result2[1]);
      } catch (IOException e){
         fail(e.getMessage());
      }


   }

   private void testSaxParser(String xmlNew, String name1, String name2,
         String[] result1, Atom[] ca1, Atom[] ca3, String name3,
         String name4, String[] result2, Atom[] ca4, Atom[] ca6) {




   }

   private AFPChain getAfpFor(String name1, AFPChain[] chains) {

      for ( AFPChain c : chains){
         //System.out.println("comparing with " + c.getName1());
         if (c.getName1().equals(name1))
            return c;
      }
      return null;
   }
}
