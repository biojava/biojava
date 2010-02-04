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
 * Created on Nov 2, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.ce;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;

import java.io.InputStreamReader;
import java.io.PrintStream;

import java.net.URISyntaxException;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.CliTools;
import org.biojava.bio.structure.align.util.ConfigurationException;

import org.biojava.bio.structure.align.xml.AFPChainXMLConverter;
import org.biojava.bio.structure.io.PDBFileReader;



public abstract class AbstractUserArgumentProcessor implements UserArgumentProcessor {

   public static String newline = System.getProperty("line.separator");

   protected StartupParameters params ;

   public static final List<String> mandatoryArgs= new ArrayList<String>();

   /** the system property PDB_DIR can be used to configure the 
    * default location for PDB files.
    */
   public static final String PDB_DIR = "PDB_DIR";

   protected AbstractUserArgumentProcessor(){ 
      params = new StartupParameters();
   }


   public abstract StructureAlignment getAlgorithm();
   public abstract Object getParameters();
   public abstract String getDbSearchLegend();

   public void process(String[] argv){

      List<String> mandatoryArgs = getMandatoryArgs();

      for (int i = 0 ; i < argv.length; i++){
         String arg   = argv[i];

         String value = null;
         if ( i < argv.length -1)
            value = argv[i+1];

         // if value starts with - then the arg does not have a value.
         if (value != null && value.startsWith("-"))
            value = null;
         else 
            i++;


         String[] tmp = {arg,value};

         //System.out.println(arg + " " + value);

         try {

            CliTools.configureBean(params, tmp);  

         } catch (ConfigurationException e){

            e.printStackTrace();

            if ( mandatoryArgs.contains(arg) ) {
               // there must not be a ConfigurationException with mandatory arguments.
               return;
            } else {
               // but there can be with optional ...
            }
         }           
      }

      if ( params.showMenu){

         if ( params.getPdbFilePath() != null){
            System.setProperty(PDB_DIR,params.getPdbFilePath());
         }

         try {
            GuiWrapper.showAlignmentGUI();
         } catch (Exception e){
            System.err.println(e.getMessage());
            e.printStackTrace();
         }
      }

      if ( params.getShowDBresult() != null){
         // user wants to view DB search results:

         System.err.println("not implemented full yet");
         //DBResultTable table = new DBResultTable();
         //UserConfiguration config = UserConfiguration.fromStartupParams(params);
         //table.show(new File(params.getShowDBresult()),config);
      }

      String pdb1  = params.getPdb1();
      String file1 = params.getFile1();
      
     
      if (pdb1 != null || file1 != null){
         runPairwise();
      }

      if ( params.getAlignPairs() != null){
         runDBSearch();
      }
   }


   private void runDBSearch(){



      String pdbFilePath = params.getPdbFilePath();

      if ( pdbFilePath == null || pdbFilePath.equals("")){

         System.err.println("Please specify mandatory argument -pdbFilePath!");
         return;

      }


      AtomCache cache = new AtomCache(pdbFilePath, params.isPdbDirSplit());

      String inputFile = params.getAlignPairs();

      if ( inputFile == null || inputFile.equals("")){
         System.err.println("Please specify the mandatory argument -inputFile!");
         return;
      }

      String outputFile = params.getOutFile();

      if ( outputFile == null || outputFile.equals("")){
         System.err.println("Please specify the mandatory argument -outFile!");
         return;
      }

      System.out.println("running DB search with parameters:" + params);


      try {
         File f = new File(inputFile);

         BufferedReader is = new BufferedReader (new InputStreamReader(new FileInputStream(f)));

         BufferedWriter out = new BufferedWriter(new FileWriter(outputFile, true));
         out.write("#Legend: " + newline );
         String legend = getDbSearchLegend();
         out.write(legend + newline );
         System.out.println(legend);
         String line = null;
         while ( (line = is.readLine()) != null){
            if ( line.startsWith("#"))
               continue;

            String[] spl = line.split(" ");

            if ( spl.length != 2) {
               System.err.println("wrongly formattted line. Expected format: 4hhb.A 4hhb.B but found " + line);
               continue;
            } 

            String pdb1 = spl[0];
            String pdb2 = spl[1];


            Structure structure1 = cache.getStructure(pdb1);
            Structure structure2 = cache.getStructure(pdb2);

            Atom[] ca1;
            Atom[] ca2;


            ca1 = StructureTools.getAtomCAArray(structure1);
            ca2 = StructureTools.getAtomCAArray(structure2);


            StructureAlignment algorithm =  getAlgorithm();
            Object jparams = getParameters();

            
            AFPChain afpChain;

            afpChain = algorithm.align(ca1, ca2, jparams);
            afpChain.setName1(pdb1);
            afpChain.setName2(pdb2);

            String result = getDbSearchResult(afpChain);
            out.write(result);
            System.out.print(result);
         }

         out.close();
         is.close();
      } catch(Exception e){
         e.printStackTrace();
      }
   }


   private void runPairwise(){

      String name1 = params.getPdb1();
      String file1 = params.getFile1();

      if ( name1 == null && file1 == null){
         throw new RuntimeException("You did not specify the -pdb1 or -file1 parameter. Can not find query PDB id for alignment.");
      }

      if ( file1 == null) {
         if ( name1.length() < 4) {
            throw new RuntimeException("-pdb1 does not look like a PDB ID. Please specify PDB code or PDB.chainId.");
         }
      }

      String name2 = params.getPdb2();
      String file2 = params.getFile2();
      if ( name2 == null && file2 == null ){
         throw new RuntimeException("You did not specify the -pdb2 or -file2 parameter. Can not find target PDB id for alignment.");
      }

      if ( file2 == null ){
         if ( name2.length() < 4) {
            throw new RuntimeException("-pdb2 does not look like a PDB ID. Please specify PDB code or PDB.chainId.");
         }
      }

      // first load two example structures

      Structure structure1 = null;
      Structure structure2 = null;

      String path = params.getPdbFilePath();

      if ( file1 == null || file2 == null) {
         if ( path == null){
            throw new RuntimeException("You did not specify the -pdbFilePath parameter. Can not find the PDB files in your file system.");
         }

         AtomCache cache = new AtomCache(path, params.isPdbDirSplit());
         cache.setAutoFetch(params.isAutoFetch());   
         structure1 = getStructure(cache, name1, file1);
         structure2 = getStructure(cache, name2, file2);
      } else {
         
         structure1 = getStructure(null, name1, file1);
         structure2 = getStructure(null, name2, file2);
      }      

      if ( structure1 == null){
         System.err.println("structure 1 is null, can't run alignment.");
         return;
      }

      if ( structure2 == null){
         System.err.println("structure 2 is null, can't run alignment.");
         return;
      }

      //                   default:      new:
      // 1buz - 1ali : time: 8.3s eqr 68 rmsd 3.1 score 161 | time 6.4 eqr 58 rmsd 3.0 scre 168
      // 5pti - 1tap : time: 6.2s eqr 48 rmsd 2.67 score 164 | time 5.2 eqr 49 rmsd 2.9 score 151
      // 1cdg - 8tim
      // 1jbe - 1ord
      // 1nbw.A - 1kid
      // 1t4y - 1rp5


      try {


         Atom[] ca1;
         Atom[] ca2;


         List<Group> hetatms1 = new ArrayList<Group>();
         List<Group> nucs1    = new ArrayList<Group>();
         List<Group> hetatms2 = new ArrayList<Group>();
         List<Group> nucs2    = new ArrayList<Group>();

         ca1 = StructureTools.getAtomCAArray(structure1);
         ca2 = StructureTools.getAtomCAArray(structure2);

         if ( params.isShow3d()){
            hetatms1 = structure1.getChain(0).getAtomGroups("hetatm");
            nucs1    = structure1.getChain(0).getAtomGroups("nucleotide");

         }

         StructureAlignment algorithm =  getAlgorithm();
         Object jparams = getParameters();

         AFPChain afpChain;

         afpChain = algorithm.align(ca1, ca2, jparams);
         afpChain.setName1(name1);
         afpChain.setName2(name2);

         if ( params.isShow3d()){

            if ( (afpChain.getBlockNum() - 1) == 0){
               hetatms2 = structure2.getChain(0).getAtomGroups("hetatm");
               nucs2    = structure2.getChain(0).getAtomGroups("nucleotide");
            }

            if (! GuiWrapper.isGuiModuleInstalled()) {
               System.err.println("The biojava-structure-gui module is not installed. Please install!");
            } else {
               try {
                  Object jmol = GuiWrapper.display(afpChain,ca1,ca2,hetatms1, nucs1, hetatms2, nucs2);
                  GuiWrapper.showAlignmentImage(afpChain, ca1,ca2,jmol);
               } catch (Exception e){
                  System.err.println(e.getMessage());
                  e.printStackTrace();
               }
               //StructureAlignmentJmol jmol = algorithm.display(afpChain,ca1,ca2,hetatms1, nucs1, hetatms2, nucs2);

               //String result = afpChain.toFatcat(ca1, ca2);
               //String rot = afpChain.toRotMat();

               //DisplayAFP.showAlignmentImage(afpChain, ca1,ca2,jmol);
            }
         } 

         if ( params.getOutFile() != null ){

            // write the XML to a file...
            String fatcatXML = AFPChainXMLConverter.toXML(afpChain,ca1,ca2);
            FileOutputStream out; // declare a file output object
            PrintStream p; // declare a print stream object

            try
            {
               // Create a new file output stream
               out = new FileOutputStream(params.getOutFile());

               // Connect print stream to the output stream
               p = new PrintStream( out );

               p.println (fatcatXML);

               p.close();
            }
            catch (Exception e)
            {
               System.err.println ("Error writing to file " + params.getOutFile());
            }
         }

         if ( params.isPrintXML()){
            String fatcatXML = AFPChainXMLConverter.toXML(afpChain,ca1,ca2);
            System.out.println(fatcatXML);
         }
         if ( params.isPrintFatCat()) {
            // default output is to XML on sysout...
            System.out.println(afpChain.toFatcat(ca1, ca2));
         }
         if ( params. isPrintCE()){
            System.out.println(afpChain.toCE(ca1, ca2));
         }


      } catch (Exception e) {
         e.printStackTrace();
         return;
      }
   }

   private Structure getStructure(AtomCache cache, String name1, String file) 
   {

      PDBFileReader reader = new PDBFileReader();
      if ( file != null ){
         try {
            // check if it is a URL:
            try {
               URL url = new URL(file);
               System.out.println(url);
                              
               return reader.getStructure(url);
               
            } catch ( Exception e){
               System.err.println(e.getMessage());
            }
            File f= new File(file);
            System.out.println("file from local " + f.getAbsolutePath());
            return reader.getStructure(f);
         } catch (Exception e){
            System.err.println("general exception:" + e.getMessage());
            System.err.println("unable to load structure from " + file);
            return null;
         }
      }
      try {
         Structure s = cache.getStructure(name1);
         return s;
      } catch ( Exception e){
         System.err.println(e.getMessage());
         System.err.println("unable to load structure from dir: " + cache.getPath() + "/"+ name1);
         return null;
      }

   }


   public List<String> getMandatoryArgs() {
      return mandatoryArgs;
   }

   public String getDbSearchResult(AFPChain afpChain){
      return afpChain.toDBSearchResult();
   }





}
