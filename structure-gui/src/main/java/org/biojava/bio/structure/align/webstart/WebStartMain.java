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
 * Created on Jun 13, 2009
 * Created by Andreas Prlic
 *
 */

package org.biojava.bio.structure.align.webstart;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JTextField;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.AFPTwister;
import org.biojava.bio.structure.align.StrucAligParameters;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.StructurePairAligner;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.client.FarmJobParameters;
import org.biojava.bio.structure.align.client.JFatCatClient;
import org.biojava.bio.structure.align.client.PdbPair;
import org.biojava.bio.structure.align.fatcat.FatCatRigid;
import org.biojava.bio.structure.align.gui.AlignmentGui;
import org.biojava.bio.structure.align.gui.ChooseDirAction;
import org.biojava.bio.structure.align.gui.DisplayAFP;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.gui.BiojavaJmol;
import org.biojava.bio.structure.gui.util.PDBDirPanel;
import org.biojava.bio.structure.io.PDBFileReader;



public class WebStartMain
{


   static UserConfiguration userConfig;

   /** 
    *  If no arguments, shows AlignmentGui for pairwise alignments.
    *  if 1 argument: dbmenu shows the DBSearchGUI, otherwise the AlignentGUI is shown.
    *  
    * if more than 3 arguments takes the following arguments
    * arg0 : fatcat or biojava . 
    * arg1 : pdb1.X
    * arg2 ; pdb2.X
    * 
    * The 4th argument is optional and it could be the serverLocation which the client should talk to.
    *
    * @param args
    */
   public static void main(String[] args){
      
      AligUIManager.setLookAndFeel();

      if ( args.length == 0){

         //JOptionPane.showMessageDialog(null,
         //		"Not enough arguments. Need at least 3, but only got " + args.length);

         // we did not get enough arguments, show the general user interface...

         javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
               
               AlignmentGui.getInstance();
            }
         });

         return;

      }

      else if ( args.length < 3){
         //String arg0 = args[0];

         javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
               AlignmentGui.getInstance();
            }
         });

         return;

      }

      String arg0 = args[0];

      if (! (arg0.equals("fatcat") || 
            arg0.equals("biojava") ||
            arg0.equals("fatcat_flexible") ||
            arg0.equals("ce") ||
            arg0.equals("sw")
      )){
         JOptionPane.showMessageDialog(null,
               "Wrong arguments. First argument has to be \"fatcat\", \"ce\", ,\"sw\", \"fatcat_flexible\", or \"biojava\", but got " + arg0);
         return;
      }

      String serverLocation = FarmJobParameters.DEFAULT_SERVER_URL;

      if ( args.length  > 3 ) {
         // we have 4 arguments.
         // in this case the 4th has to be the server URL
         serverLocation = args[3];
      }

      try {

         String name1 = args[1];
         String name2 = args[2];

         // check arguments
         getPdbChain(name1);
         getPdbChain(name2);

         PdbPair pair = new PdbPair(name1, name2);
         System.out.println("### user provided: " + pair);

         UserConfiguration config = getWebStartConfig();

         System.setProperty("PDB_DIR",config.getPdbFilePath());
         System.out.println("using PDB file path: " + config.getPdbFilePath());

         AtomCache cache = new AtomCache(config);
         
         JFrame frame = new JFrame();
         frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

         showProgressBar(frame,"Loading PDB files...", "Loading files from local directory or via FTP.");

         Structure structure1 = cache.getStructure(pair.getPDBCode1());
         Structure structure2 = cache.getStructure(pair.getPDBCode2());

         frame.dispose();

         System.out.println("done reading structures");

         String chainId1 = pair.getChainId1();
         String chainId2 = pair.getChainId2();

         if ( ! structure1.hasChain(chainId1)){
            System.out.println("no chain id " + chainId1 + " found. trying upper case...");
            chainId1 = chainId1.toUpperCase();
         }

         if ( ! structure2.hasChain(chainId2)){
            System.out.println("no chain id " + chainId2 + " found. trying upper case...");
            chainId2 = chainId1.toUpperCase();
         }

         Chain c1 = structure1.getChainByPDB(chainId1);
         Chain c2 = structure2.getChainByPDB(chainId2);

         if (arg0.equalsIgnoreCase("ce") || 
               arg0.equalsIgnoreCase("sw") ||
               arg0.equalsIgnoreCase("fatcat_flexible")){
            try {

               StructureAlignment algorithm ;
               if ( arg0.equalsIgnoreCase("ce"))
                  algorithm = new CeMain();
               else if ( arg0.equalsIgnoreCase("fatcat_flexible"))
                  algorithm = StructureAlignmentFactory.getAlgorithm("jFatCat_flexible");
               else
                  algorithm = new SmithWaterman3Daligner();
               showStructureAlignment(serverLocation,algorithm ,c1,c2, pair.getName1(),pair.getName2());
            } catch (Exception e){
               e.printStackTrace();
               JOptionPane.showMessageDialog(null,
                     "Something went wrong! : " + e.getMessage());
            }
         } else if ( arg0.equalsIgnoreCase("fatcat")){
            try {
               showFatcat(serverLocation,c1,c2, pair.getName1(),pair.getName2());
            } catch (Exception e){
               e.printStackTrace();
               JOptionPane.showMessageDialog(null,
                     "Something went wrong! : " + e.getMessage());
            }
         } else if ( arg0.equalsIgnoreCase("biojava")){
            try {
               showBiojava(structure1,c1,structure2,c2);
            } catch (Exception e){
               e.printStackTrace();
               JOptionPane.showMessageDialog(null,
                     "Something went wrong! : " + e.getMessage());
            }
         }


      } catch (Exception e) {
         e.printStackTrace();
         JOptionPane.showMessageDialog(null,
               "Error: " + e.getMessage());
         return;
      }
   }

   private static JProgressBar showProgressBar(JFrame frame, String title, String displayTxt){

      frame.setTitle(title);
      JProgressBar progressBar;

      JPanel content = new JPanel();
      content.setOpaque(true);

      //Where the GUI is constructed:
      progressBar = new JProgressBar();
      progressBar.setToolTipText(title);
      progressBar.setValue(0);
      progressBar.setStringPainted(true);
      progressBar.setIndeterminate(true);

      JTextField txt = new JTextField(displayTxt);
      txt.setEditable(false);
      content.add(txt);
      content.add(progressBar);

      frame.getContentPane().add(content);
      frame.pack();
      frame.setVisible(true);

      return progressBar;
   }


   public static UserConfiguration getWebStartConfig(){

      if ( userConfig == null) {
         try {
            PersistentConfig webstartConfig = new PersistentConfig();

            userConfig = webstartConfig.load();

         } catch (Exception e){
            System.err.println(e.getMessage());
         }
      }

      // check if we could load it (i.e. we are running in web start mode)
      if ( userConfig == null ) {
         userConfig = WebStartMain.getDefaultConfig();

         try {
            PersistentConfig webstartConfig = new PersistentConfig();

            webstartConfig.save(userConfig);

         } catch (Exception e){
            System.err.println(e.getMessage());
         }
      }

      return userConfig;
   }


   public static UserConfiguration getDefaultConfig(){
      userConfig = new UserConfiguration();

      String pdbDir = System.getProperty(PDBDirPanel.PDB_DIR);
      if ( pdbDir != null) {
         userConfig.setPdbFilePath(pdbDir);

      }

      return userConfig;
   }

   public static void persistConfig(UserConfiguration config){

      try {
          PersistentConfig webstartConfig = new PersistentConfig();

          webstartConfig.save(config);

      } catch (Exception e){
          e.printStackTrace();
      }

  }
   
   public static UserConfiguration requestUserConfig(){

      if ( userConfig == null) {

         //UserConfiguration config = new UserConfiguration();
         userConfig = new UserConfiguration();

         String pdbDir = System.getProperty(PDBDirPanel.PDB_DIR);
         if ( pdbDir != null) {
            userConfig.setPdbFilePath(pdbDir);
            return userConfig;
         }

      }
      JTextField textField = new JTextField();
      ChooseDirAction action = new ChooseDirAction(textField, userConfig);
      action.actionPerformed(null);
      if ( textField.getText() == null) {

         // accessing temp. OS directory:         
         String property = "java.io.tmpdir";

         String tempdir = System.getProperty(property);

         if ( !(tempdir.endsWith(PDBFileReader.lineSplit ) ) )
            tempdir = tempdir + PDBFileReader.lineSplit;

         userConfig.setPdbFilePath(tempdir);
         return userConfig;
      }

      File file = new File(textField.getText());
      if ( ! file.isDirectory() ){
         // should not happen
         System.err.println("did not provide directory, going on level higher! " + file.getAbsolutePath());
         file = file.getParentFile();
      }
      System.setProperty(PDBDirPanel.PDB_DIR, file.getAbsolutePath());
      userConfig.setPdbFilePath(file.getAbsolutePath());

      return userConfig;
   }

   private static String[] getPdbChain(String arg) throws StructureException{
      String[] spl = arg.split("\\.");
      if ( spl.length!=2)
         throw new StructureException("Argument does not look like a pdb code and chain ID. expected format: 1abc.A (chains are case sensitive)");
      if ( spl[0].length() != 4 ){
         throw new StructureException("Argument does not look like a pdb code and chain ID. expected format: 1abc.A (chains are case sensitive)");
      }
      if ( spl[1].length() != 1 ){
         throw new StructureException("Argument does not look like a pdb code and chain ID. expected format: 1abc.A (chains are case sensitive)");
      }
      return spl;
   }

   private static void showStructureAlignment(String serverLocation, StructureAlignment algorithm, Chain c1, Chain c2, String name1, String name2) throws StructureException{
      JFrame tmpFrame = new JFrame();
      tmpFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

      String title = "Calculating " + algorithm.getAlgorithmName() + " V." + algorithm.getVersion()+" alignment... ";
     

      showProgressBar(tmpFrame,title, "Calculating the structure alignment.");

      Atom[] ca1 = StructureTools.getAtomCAArray(c1);
      Atom[] ca2 = StructureTools.getAtomCAArray(c2);

      
      //do the actual alignment 
      AFPChain afpChain = null;
      
      try {
         // using 10 sec as timeout on server now, since we expect the server to be able to complete the calculation within that time...
         afpChain =  JFatCatClient.getAFPChainFromServer(serverLocation,algorithm.getAlgorithmName(), name1, name2, ca1, ca2, 10000);
      } catch (Exception e){
         e.printStackTrace();
      }
      
      if ( afpChain == null )  {
         System.out.println(title);
         afpChain = algorithm.align(ca1, ca2);
      }

      afpChain.setName1(name1);
      afpChain.setName2(name2);

      tmpFrame.dispose();

      List<Group> hetatms = c1.getAtomGroups("hetatm");
      List<Group> nucs    = c1.getAtomGroups("nucleotide");

      List<Group> hetatms2 = new ArrayList<Group>();
      List<Group> nucs2    = new ArrayList<Group>();

     // System.out.println("got afpChain blocknum: " + afpChain.getBlockNum());
      if ( (afpChain.getBlockNum() - 1) == 0){
         hetatms2 = c2.getAtomGroups("hetatm");
         nucs2    = c2.getAtomGroups("nucleotide");
      }

      // show results
      StructureAlignmentJmol jmol = StructureAlignmentDisplay.display(afpChain,ca1,ca2,hetatms,nucs, hetatms2, nucs2);

      //jmol.setTitle(title);
      // make sure the AlignmentGUI is always displayed, because it is the only way to close the application.
      //AlignmentGui.getInstance();


      //String result = afpChain.toFatcat(ca1, ca2);
      //String rot = afpChain.toRotMat();

      System.out.println(afpChain.toCE(ca1, ca2));


      DisplayAFP.showAlignmentImage(afpChain,ca1,ca2,jmol);

   }
   private static void showFatcat(String serverLocation, Chain c1, Chain c2, String name1, String name2) throws StructureException{

      JFrame tmpFrame = new JFrame();
      tmpFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

      showProgressBar(tmpFrame,"Calculating JFatCat (rigid) alignment... ", "Calculating the structure alignment.");

      StructureAlignment fatCatRigid = StructureAlignmentFactory.getAlgorithm(FatCatRigid.algorithmName);

      Atom[] ca1 = StructureTools.getAtomCAArray(c1);
      Atom[] ca2 = StructureTools.getAtomCAArray(c2);

      // check if we already have a precalculated alignment...
      AFPChain afpChain = null;

      try {
         //System.out.println("requesting alignment from server...");

         afpChain = JFatCatClient.getAFPChainFromServer(serverLocation, name1, name2, ca1, ca2);

         if ( afpChain != null){
            // twist groups for visualisation
            // we only have data for the optimized alignment so far...

            //  System.out.println("sending afpChain to fatcatRigid proxy");

            //JFatCatProxy proxy = new JFatCatProxy();
            //proxy.setStructureAlignment(fatCatRigid);
            //proxy.twistGroups(afpChain,ca1,ca2);

         } else {
            System.out.println("calculating");
            afpChain = fatCatRigid.align(ca1,ca2);

            afpChain.setName1(name1);
            afpChain.setName2(name2);

            //System.out.println("submitting...");
            //JFatCatClient.sendAFPChainToServer(serverLocation,afpChain, ca1, ca2);
         }
      } catch (Exception e){
         e.printStackTrace();
      }

      if ( afpChain == null) {
         //do the actual alignment 
         afpChain = fatCatRigid.align(ca1, ca2);
         afpChain.setName1(name1);
         afpChain.setName2(name2);

      }


      tmpFrame.dispose();


      Group[] twistedGroups = AFPTwister.twistOptimized(afpChain,ca1,ca2);

      List<Group> hetatms = ca1[0].getParent().getParent().getAtomGroups("hetatm");
      List<Group> nucs    = ca1[0].getParent().getParent().getAtomGroups("nucleotide");

      List<Group> hetatms2 = new ArrayList<Group>();
      List<Group> nucs2    = new ArrayList<Group>();


      if ( (afpChain.getBlockNum() - 1) == 0){           
         hetatms2 = ca2[0].getParent().getParent().getAtomGroups("hetatm");
         nucs2    = ca2[0].getParent().getParent().getAtomGroups("nucleotide");

         // they are not rotated at this point.. the display will do it for us...

      }
      StructureAlignmentJmol jmol = DisplayAFP.display(afpChain, twistedGroups, ca1, ca2, hetatms, nucs, hetatms2, nucs2);
      DisplayAFP.showAlignmentImage(afpChain,ca1,ca2,jmol);

   }



   private static void showBiojava(Structure structure1, Chain c1,Structure structure2, Chain c2) throws StructureException{

      StructurePairAligner aligner = new StructurePairAligner();

      StrucAligParameters params = StrucAligParameters.getDefaultParameters();
      //params.setJoinFast(true);
      //params.setMaxIter(1);
      try {
         // align the full 2 structures with default parameters.
         // see StructurePairAligner for more options and how to align
         // any set of Atoms
         long start = System.currentTimeMillis();

         Structure s3 = new StructureImpl();
         s3.addChain(c1);
         Structure s4 = new StructureImpl();
         s4.addChain(c2);

         aligner.align(s3,s4,params);
         long end = System.currentTimeMillis();
         System.out.println("calculation time:" + (end-start));

         AlternativeAlignment[] aligs = aligner.getAlignments();
         AlternativeAlignment a = aligs[0];
         System.out.println(a);

         if (! BiojavaJmol.jmolInClassPath()){
            System.err.println("Could not find Jmol in classpath, please install first!");
            return;
         }
         // display the alignment in Jmol

         // first get an artificial structure for the alignment
         Structure artificial = a.getAlignedStructure(structure1, structure2);

         // and then send it to Jmol (only will work if Jmol is in the Classpath)
         BiojavaJmol jmol = new BiojavaJmol();
         jmol.setTitle(artificial.getName());
         jmol.setStructure(artificial);

         // color the two structures


         jmol.evalString("select *; backbone 0.4; wireframe off; spacefill off; " +
         "select not protein and not solvent; spacefill on;");
         jmol.evalString("select */1 ; color red; model 1; ");


         // now color the equivalent residues ...

         String[] pdb1 = a.getPDBresnum1();
         for (String res : pdb1 ){
            jmol.evalString("select " + res + "/1 ; backbone 0.6; color orange;");
         }

         jmol.evalString("select */2; color blue; model 2;");
         String[] pdb2 = a.getPDBresnum2();
         for (String res :pdb2 ){
            jmol.evalString("select " + res + "/2 ; backbone 0.6; color cyan;");
         }


         // now show both models again.
         jmol.evalString("model 0;");

      } catch (StructureException e){
         e.printStackTrace();
      }

   }
}
