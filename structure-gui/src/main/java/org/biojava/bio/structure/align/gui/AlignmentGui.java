/*
 *                  BioJava development code
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
 * Created on Jul 16, 2006
 *
 */
package org.biojava.bio.structure.align.gui;

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.util.logging.Logger;
import javax.swing.AbstractAction;
import javax.swing.Action;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuBar;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JTextField;

import javax.swing.JTabbedPane;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.gui.AlignmentGui;
import org.biojava.bio.structure.align.gui.ProgressThreadDrawer;
import org.biojava.bio.structure.align.util.ResourceManager;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.gui.util.PDBUploadPanel;
import org.biojava.bio.structure.gui.util.StructurePairSelector;


/** A JFrame that allows to trigger a pairwise structure alignment,
 * either from files in a directory,
 * or after manual upload.
 *
 * @author Andreas Prlic
 *
 * @since 1.7
 *
 *
 *
 */
public class AlignmentGui extends JFrame{

   private final static long serialVersionUID =0l;

   public static Logger logger =  Logger.getLogger("org.biojava.spice");

   StructureAlignment algorithm;

   JButton abortB;

   SelectPDBPanel tab1 ;
   PDBUploadPanel tab2;

   Thread thread;
   AlignmentCalculationRunnable alicalc;
   JTabbedPane masterPane;
   JTabbedPane tabPane;
   JProgressBar progress;

   private JTextField pdbDir;

   private DBSearchGUI dbsearch;


   public static void main(String[] args){
      AlignmentGui.getInstance();

   }

   static final ResourceManager resourceManager = ResourceManager.getResourceManager("ce");

   private static final String MAIN_TITLE = "Pairwise Structure Alignment - Main - V." + resourceManager.getString("ce.version");;

   private static final AlignmentGui me = new AlignmentGui();

   public static AlignmentGui getInstance(){

      if (!  me.isVisible())
         me.setVisible(true);

      if ( ! me.isActive())
         me.requestFocus();


      return me;
   }

   public static AlignmentGui getInstanceNoVisibilityChange(){

      return me;
   }


   private AlignmentGui() {
      super();

      thread = null;

      JMenuBar menu = MenuCreator.initAlignmentGUIMenu(this);

      this.setJMenuBar(menu);

      this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

      this.setTitle(MAIN_TITLE);

      /*this.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent evt) {
				JFrame frame = (JFrame) evt.getSource();
				frame.setVisible(false);
				frame.dispose();
			}
		});*/

      String[] algorithms = StructureAlignmentFactory.getAllAlgorithmNames();
      try {
         algorithm = StructureAlignmentFactory.getAlgorithm(algorithms[0]);
      } catch (StructureException e){
         e.printStackTrace();
      }

      JLabel algoLabel = new JLabel("Select alignment algorithm: ");

      JComboBox algorithmList = new JComboBox(algorithms);
      algorithmList.setSelectedIndex(0);

      Action actionAlgorithm = new AbstractAction("Algorithm") {
         public static final long serialVersionUID = 0l;
         // This method is called when the button is pressed
         public void actionPerformed(ActionEvent evt) {
            JComboBox cb = (JComboBox)evt.getSource();
            String algorithmName = (String)cb.getSelectedItem();
            // Perform action...
            //System.out.println("calc structure alignment");
            updateAlgorithm(algorithmName);

         }
      };


      algorithmList.addActionListener(actionAlgorithm);


      Action paramAction = new AbstractAction("Parameters") {
         public static final long serialVersionUID = 0l;
         // This method is called when the button is pressed
         public void actionPerformed(ActionEvent evt) {
            // Perform action...
            //System.out.println("calc structure alignment");
            configureParameters();

         }

      };

      JButton parameterButton = new JButton(paramAction);


      tab1 = new SelectPDBPanel();
      //tab1 = new PDBDirPanel();

      tab2 = new PDBUploadPanel();

      //tab3 = new PDBServerPanel();

      tabPane = new JTabbedPane();

      tabPane.addTab("Select PDB ID", null, tab1,
      "select PDB ID to align");

      tabPane.addTab("Custom files",null, tab2,"Align your own files.");

      //tabPane.addTab("Local PDB config",null,tab3,"Configure local PDB installation.");



      // BUILD UP THE UI out of the components defined above.

      Box vBox = Box.createVerticalBox();

      Box hBoxAlgo = Box.createHorizontalBox();
      hBoxAlgo.add(Box.createGlue());
      hBoxAlgo.add(algoLabel);		
      hBoxAlgo.add(algorithmList);
      hBoxAlgo.add(Box.createGlue());
      hBoxAlgo.add(parameterButton);
      hBoxAlgo.add(Box.createGlue());
      //vBox.add(hBoxAlgo);


      vBox.add(tabPane);

      //		Box hBox42 = Box.createHorizontalBox();
      progress =new JProgressBar();
      progress.setIndeterminate(false);
      progress.setMaximumSize(new Dimension(10,100));


      //		hBox42.add(Box.createGlue());
      //		hBox42.add(progress);		
      //		hBox42.add(Box.createGlue());
      //		vBox.add(hBox42);
      Action action1 = new AbstractAction("Align") {
         public static final long serialVersionUID = 0l;
         // This method is called when the button is pressed
         public void actionPerformed(ActionEvent evt) {
            // Perform action...
            //System.out.println("calc structure alignment");
            int selectedIndex = masterPane.getSelectedIndex();
            if (selectedIndex == 0)
               calcAlignment();
            else if ( selectedIndex == 1)
               calcDBSearch();
            else {
               System.err.println("Unknown TAB: " + selectedIndex);
            }

         }
      };

      JButton submitB = new JButton(action1);

      Action action3 = new AbstractAction("Abort") {
         public static final long serialVersionUID = 0l;
         // This method is called when the button is pressed
         public void actionPerformed(ActionEvent evt) {
            // Perform action...
            abortCalc();
         }
      };
      
      abortB = new JButton(action3);

      abortB.setEnabled(false);

      Action action2 = new AbstractAction("Exit") {
         public static final long serialVersionUID = 0l;
         // This method is called when the button is pressed
         public void actionPerformed(ActionEvent evt) {
            // Perform action...
            abortCalc();
            dispose();
            System.exit(0);
         }
      };
      JButton closeB = new JButton(action2);
      Box hBox = Box.createHorizontalBox();
      hBox.add(closeB);
      hBox.add(Box.createGlue());
      hBox.add(progress);
      hBox.add(abortB);
      //hBox.add(Box.createGlue());
      hBox.add(submitB);


      //vBox.setBorder(BorderFactory.createEmptyBorder(20, 20, 20, 20));

      masterPane = new JTabbedPane();
      masterPane.addTab("Pairwise Comparison", vBox);

      pdbDir = new JTextField(20);
      tab1.setPDBDirField(pdbDir);
      
      dbsearch = new DBSearchGUI();
      dbsearch.getSelectPDBPanel().setPDBDirField(pdbDir);
     
      masterPane.addTab("Database Search",dbsearch);



      

      String conf = System.getProperty(SelectPDBPanel.PDB_DIR);
      if ( conf != null){
         pdbDir.setText(conf);
      }
      JPanel dir = tab1.getPDBDirPanel(pdbDir);

      //JTabbedPane configPane = new JTabbedPane();

      //configPane.addTab("Local PDB install", null, dir,
      //"Configure your local PDB setup.");

      //vBox.add(configPane);


      Box vBoxMain = Box.createVerticalBox();

      // pairwise or db search
      vBoxMain.add(masterPane);
      // algorithm selection
      vBoxMain.add(hBoxAlgo);
      // PDB install config
      vBoxMain.add(dir);
      // buttons
      vBoxMain.add(hBox);

      this.getContentPane().add(vBoxMain);
      this.pack();
      this.setVisible(true);
   }




   protected void configureParameters() {
      StructureAlignment algorithm = getStructureAlignment();
      System.out.println("configure parameters for " + algorithm.getAlgorithmName());

      // show a new config GUI
      new ParameterGUI(algorithm);

   }



   public void cleanUp() {

      if ( alicalc != null) {
         alicalc.cleanup();
      }
   }



   private void calcAlignment() {

      int pos = tabPane.getSelectedIndex();
      StructurePairSelector tab = null;

      if (pos == 0 ){
         tab = tab1;
         tab1.persistCurrentConfig();
      } else if (pos == 1){
         tab = tab2;
      }
      try {
         Structure s1 = tab.getStructure1();
         Structure s2 = tab.getStructure2();

         if ( s1 == null) {
            System.err.println("please select structure 1");
            return ;
         }

         if ( s2 == null) {
            System.err.println("please select structure 2");
            return;
         }

         String name1 = "custom1"; 
         String name2 = "custom2";

         if  ( pos == 0){
            name1 = tab1.getName1();
            name2 = tab1.getName2();
         } else {
            name1 = s1.getPDBCode();
            name2 = s2.getPDBCode();
         }


         alicalc = new AlignmentCalc(this,s1,s2, name1, name2);


         thread = new Thread(alicalc);
         thread.start();
         abortB.setEnabled(true);
         progress.setIndeterminate(true);
         ProgressThreadDrawer drawer = new ProgressThreadDrawer(progress);
         drawer.start();
      } catch (StructureException e){
         JOptionPane.showMessageDialog(null,"Could not align structures. Exception: " + e.getMessage());
      }

   }

   private void calcDBSearch() {

      JTabbedPane tabPane = dbsearch.getTabPane();
      System.out.println("run DB search " + tabPane.getSelectedIndex());

      Structure s = null;
      
      dbsearch.getSelectPDBPanel().setFromFtp(tab1.getFromFtp());
      dbsearch.getSelectPDBPanel().setPdbSplit(tab1.getPdbSplit());
      
      StructurePairSelector tab = null;
      int pos = tabPane.getSelectedIndex();
      
      if (pos == 0 ){
      
         tab = dbsearch.getSelectPDBPanel();
         tab1.persistCurrentConfig();
      
      } else if (pos == 1){
      
         tab = dbsearch.getPDBUploadPanel();
      
      }
      
      try {
      
         s = tab.getStructure1();

         if ( s == null) {
            System.err.println("please select structure 1");
            return ;
         }
      
      } catch (Exception e){
         e.printStackTrace();
      }
   
      String name1 = "custom1"; 
      

      if  ( pos == 0) {
         
         name1 = dbsearch.getSelectPDBPanel().getName1();
         
      } else {
         
         name1 = s.getPDBCode();
       
      }
      
      System.out.println("name1 in alig gui:" + name1);
      String file = dbsearch.getOutFileLocation();
      if ( file == null || file.equals("")){
         JOptionPane.showMessageDialog(null,"Plrease select a directory to contain the DB search results.");
         return;
      }

      UserConfiguration config = tab1.getConfiguration();
   
      alicalc = new AlignmentCalcDB(this, s,  name1,config,file);

      abortB.setEnabled(true);
      progress.setIndeterminate(true);
      ProgressThreadDrawer drawer = new ProgressThreadDrawer(progress);
      drawer.start();

      Thread t = new Thread(alicalc);
      t.start();
   }


   public void notifyCalcFinished(){
      abortB.setEnabled(false);
      thread = null;
      progress.setIndeterminate(false);
      this.repaint();
   }

   private void abortCalc(){
      System.err.println("Interrupting alignment ...");
      if ( alicalc != null )
         alicalc.interrupt();
      notifyCalcFinished();


   }


   public StructureAlignment getStructureAlignment() {

      return algorithm;
   }

   private void updateAlgorithm(String algorithmName) {

      //String algorithmName = (String)algorithmList.getSelectedItem();
      try {
         algorithm = StructureAlignmentFactory.getAlgorithm(algorithmName);
      } catch (StructureException ex){
         ex.printStackTrace();
      }

   }

}

class ProgressThreadDrawer extends Thread {

   JProgressBar progress;
   static int interval = 300;

   public ProgressThreadDrawer(JProgressBar progress) {
      this.progress = progress;
   }


   public void run() {
      boolean finished = false;
      while ( ! finished) {
         try {
            progress.repaint();
            if ( ! progress.isIndeterminate() ){
               finished =false;
               break;
            }

            sleep(interval);
         } catch (InterruptedException e){
         }
         progress.repaint();
      }
      progress = null;
   }

}
