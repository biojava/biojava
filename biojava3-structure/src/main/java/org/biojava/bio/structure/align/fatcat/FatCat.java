/* This class is based on the original FATCAT implementation by
 * <pre>
 * Yuzhen Ye & Adam Godzik (2003)
 * Flexible structure alignment by chaining aligned fragment pairs allowing twists.
 * Bioinformatics vol.19 suppl. 2. ii246-ii255.
 * http://www.ncbi.nlm.nih.gov/pubmed/14534198
 * </pre>
 * 
 * Thanks to Yuzhen Ye and A. Godzik for granting permission to freely use and redistribute this code.
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
 * Created on Jun 17, 2009
 * Created by Andreas Prlic - RCSB PDB 
 * 
 */

package org.biojava.bio.structure.align.fatcat;

import org.biojava.bio.structure.Atom;


import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.fatcat.calc.FatCatAligner;
import org.biojava.bio.structure.align.fatcat.calc.FatCatParameters;
import org.biojava.bio.structure.align.model.AFPChain;


public class FatCat
{

   public static final float VERSION = 1.0f;

   public static String newline = System.getProperty("line.separator");

   FatCatAligner aligner;




   /** See demo/FatCatDemo.java for an example how to run.
    * 
    * Launch FatCat from command line.
    * 
    * Parameters are:
    * 
    * @param argv
    */
   public static void main(String[] argv){
      FatCat cat = new FatCat();
      if (argv.length  == 0 ) {			
         System.out.println(cat.printHelp());
         return;			
      }

      if ( argv.length == 1){
         if (argv[0].equalsIgnoreCase("-h") || argv[0].equalsIgnoreCase("-help")|| argv[0].equalsIgnoreCase("--help")){
            System.out.println(cat.printHelp());								
            return;
         }
//         if ( argv[0].equalsIgnoreCase("-showMenu")){
//
//            AlignmentGui.getInstance();
//            return;
//         }
      }

      FatCatUserArgumentProcessor processor = new FatCatUserArgumentProcessor();
      processor.process(argv);


   }

   public String toString(){
      return "JFatCat v. " + VERSION;
   }


   /** print the -h help text.
    * 
    */
   public String printHelp(){
      StringBuffer buf = new StringBuffer();

      buf.append("-------------------").append(newline);
      buf.append("jFatCat v." + VERSION + " help: " + newline);
      buf.append("-------------------").append(newline);
      buf.append(newline);
      buf.append("JFatCat accepts the following parameters:").append(newline);
      buf.append(newline);
      buf.append("--- pairwise alignents ---");
      buf.append("-pdbFilePath (mandatory) Path to the directory in your file system that contains the PDB files.").append(newline);
      buf.append("-pdb1 (mandatory) PDB ID of target structure. Chain IDs are optional. In order to specify chain IDs write e.g: 5pti.A").append(newline);
      buf.append("-pdb2 (mandatory) PDB ID of query structure. Chain IDs are optional. In order to specify chain IDs write e.g: 5pti.A").append(newline);
      buf.append("-h / -help / --help : print this help string.").append(newline);
      buf.append("-printXML true/false print the XML representation of the alignment on stdout.").append(newline);
      buf.append("-printFatCat true/false print the original FATCAT output to stdout.").append(newline);
      buf.append("-printCE true/false print the result in CE style").append(newline);
      buf.append("-show3d print a 3D visualisation of the alignment (requires jmolapplet.jar in classpath)").append(newline);
      buf.append("-outFile file to write the output to (writes XML representation).").append(newline);
      buf.append("-autoFetch true/false if set to true PDB files will automatically get downloaded and stored in the right location. (default: false)").append(newline);
      buf.append("-flexible true/false run flexible alignment (default: rigid body alignment, false). ").append(newline);
      buf.append("-pdbDirSplit true/false the directory containing PDB files has all PDBs in one level or is split into multiple subdirs, like the ftp site. (default: true)").append(newline);
      buf.append("-showMenu displays the menu that allows to run alignments through a user interface.");
      buf.append(newline);
      buf.append("--- database searches ---");
      buf.append("-alignPairs (mandatory) path to a file that contains a set of pairs to compair");		
      buf.append("-outFile (mandatory) a file that will contain the summary of all the pairwise alignments");
      buf.append("-pdbFilePath (mandatory) Path to the directory in your file system that contains the PDB files.").append(newline);

      buf.append(newline);
      buf.append("For boolean arguments: if neither the text >true< or >false< is provided it is assumed to mean >true<. Instead of >-argument false< it is also possible to write -noArgument.").append(newline);
      buf.append("--- How to specify what to align ---");
      buf.append(newline);
      buf.append(" If only a PDB code is provided, the whole structure will be used for the alignment.").append(newline);
      buf.append(" To specify a particular chain write as: 4hhb.A (chain IDs are case sensitive, PDB ids are not)").append(newline);
      buf.append(" To specify that the 1st chain in a structure should be used write: 4hhb:0 .").append(newline);
      return buf.toString();
   }



   public AFPChain alignRigid(Atom[] ca1, Atom[] ca2) throws StructureException{
      StructureAlignment fatCat = new FatCatRigid();
      return fatCat.align(ca1,ca2);
   }

   public AFPChain alignRigid(Atom[] ca1, Atom[] ca2, FatCatParameters params) throws StructureException{

      AFPChain afpChain = align(ca1,ca2,params,true);
      afpChain.setAlgorithmName(FatCatRigid.algorithmName);
      afpChain.setVersion(VERSION+"");
      return afpChain;
   }

   public AFPChain alignFlexible(Atom[] ca1, Atom[] ca2, FatCatParameters params) throws StructureException{

      AFPChain afpChain = align(ca1,ca2,params,false);
      afpChain.setAlgorithmName(FatCatFlexible.algorithmName);
      afpChain.setVersion(VERSION+"");
      return afpChain;
   }


   protected AFPChain align(Atom[] ca1, Atom[] ca2, FatCatParameters params, boolean doRigid) throws StructureException{

      aligner = new FatCatAligner();

      aligner.align(ca1, ca2, doRigid, params);

      return aligner.getAfpChain();


   }

   public FatCatAligner getFatCatAligner(){
      if ( aligner == null)
         aligner = new FatCatAligner();
      return aligner;
   }

   /** Display the results of an alignment. All input coordinate are not rotated. rotations are done automatically.
    * 
    * @param afpChain
    * @param ca1
    * @param ca2
    * @param hetatms
    * @param nucs
    * @param hetatms2
    * @param nucs2
    * @throws StructureException
    */
//   public StructureAlignmentJmol display(AFPChain afpChain, Atom[] ca1, Atom[] ca2,List<Group> hetatms, List<Group>nucs,List<Group> hetatms2, List<Group>nucs2 )
//   throws StructureException {
//
//
//      FatCatAligner ali = getFatCatAligner();
//      Group[] twistedGroups = ali.getTwistedGroups();
//
//      if ( twistedGroups == null) {	
//         System.out.println("twisting");
//         Atom[] ca3 = StructureTools.cloneCAArray(ca2);
//
//         // this can happen if the alignment got loaded from a flat file
//         twistedGroups = AFPTwister.twistOptimized(afpChain, ca1, ca3);
//         ali.setTwistedGroups(twistedGroups);
//      }
//
//      return DisplayAFP.display(afpChain,twistedGroups, ca1, ca2,hetatms, nucs, hetatms2, nucs2);
//
//   }


}
