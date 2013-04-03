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
 * Created on Nov 5, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.gui;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.util.SortedSet;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.logging.Logger;
import java.util.zip.GZIPOutputStream;

import org.biojava.bio.structure.Atom;

import org.biojava.bio.structure.Structure;

import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.ce.CeSideChainMain;
import org.biojava.bio.structure.align.client.FarmJobParameters;
import org.biojava.bio.structure.align.client.JFatCatClient;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;

public class AlignmentCalcDB implements AlignmentCalculationRunnable {
   public static Logger logger =  Logger.getLogger("org.biojava");

   AtomicBoolean interrupted ;


   String name1;

   Structure structure1;

   AlignmentGui parent;

   UserConfiguration config;
   SortedSet<String> representatives;

   String outFile;
   File resultList;

   public AlignmentCalcDB(AlignmentGui parent, Structure s1,  String name1, UserConfiguration config,String outFile) {

      this.parent= parent;

      structure1 = s1;

      this.name1 = name1;

      this.config = config;
      //this.representatives = representatives;
      interrupted = new AtomicBoolean(false);
      this.outFile = outFile;
   }

   public void run() {
      System.out.println("running AlignmentCalcDB. Results will be in " + outFile);
      AtomCache cache = new AtomCache(config);
      StructureAlignment algorithm = parent.getStructureAlignment();	
      String serverLocation = FarmJobParameters.DEFAULT_SERVER_URL;
      if ( representatives == null){
         SortedSet<String> repre = JFatCatClient.getRepresentatives(serverLocation,40);
         System.out.println("got  " + repre.size() + " representatives for comparison");
         representatives = repre;
      }

      String header = "# algorithm:" + algorithm.getAlgorithmName(); 

      String legend = "# name1\tname2\tscore\tprobability\trmsd\tlen1\tlen2\tsim1\tsim2\t " ;
      if (    algorithm.getAlgorithmName().equalsIgnoreCase(CeMain.algorithmName) || 
            algorithm.getAlgorithmName().equalsIgnoreCase(CeSideChainMain.algorithmName)){
         legend =  "# name1\tname2\tscore\tz-score\trmsd\tlen1\tlen2\tsim1\tsim2\t " ;
      }


      File outFileF = new File(outFile);
      if ( ! outFileF.isDirectory()){
         System.err.println( outFileF.getAbsolutePath() + " is not a directory, can't create result files in there... ");
         interrupt();
         cleanup();
      }

      resultList = new File(outFileF,"results_" + name1 + ".out");
      BufferedWriter out;
      try {
         structure1 = cache.getStructure(name1);
         out = new BufferedWriter(new FileWriter(resultList));
         out.write(header);
         out.write(AFPChain.newline);
         out.write(legend);
         out.write(AFPChain.newline);
      } catch (Exception e){
         System.err.println("Error while loading representative structure " + name1);
         e.printStackTrace();
         interrupt();
         cleanup();
         return;
      }

      for (String repre : representatives){
         if ( interrupted.get()) {
            System.err.println("User interrupted alignments.");
            break;
         }
         try {
            
            Structure structure2 = cache.getStructure(repre);

            Atom[] ca1;
            Atom[] ca2;

            ca1 = StructureTools.getAtomCAArray(structure1);
            ca2 = StructureTools.getAtomCAArray(structure2);

            AFPChain afpChain;

            afpChain = algorithm.align(ca1, ca2);
            afpChain.setName1(name1);
            afpChain.setName2(repre);

            String result = afpChain.toDBSearchResult();
            System.out.print(result);

            out.write(result);

            String xml = AFPChainXMLConverter.toXML(afpChain, ca1, ca2);
            writeXML(outFileF,name1, repre, xml);

         } catch ( Exception e){
            e.printStackTrace();
         }

      }

      try {
         out.flush();
         out.close();
      } catch (Exception e) {
         e.printStackTrace();
      }
      parent.notifyCalcFinished();
      DBResultTable table = new DBResultTable();
      table.show(resultList,config);
   }

   private void writeXML(File outFileF, String name1, String name2, String xml)
   {
      try{
         // Create file 
         File newF = new File(outFileF, "dbsearch_" +name1+"_" + name2+".xml.gz");
                  
         FileOutputStream fstream = new FileOutputStream(newF);
        
         GZIPOutputStream gz = new GZIPOutputStream(fstream);
         OutputStreamWriter writer = new OutputStreamWriter(gz);
         writer.write(xml);
         writer.close();
      }catch (Exception e){//Catch exception if any
         System.err.println("Error: " + e.getMessage());
      }

   }

   /** stops what is currently happening and does not continue
    * 
    *
    */
   public void interrupt() {
      interrupted.set(true);
     
   }

   public void cleanup()
   {
      parent.notifyCalcFinished();
     
      parent=null;
      // cleanup...

      structure1 = null;
      config = null;
      
      
      
   }

}
