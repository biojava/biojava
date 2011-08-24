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
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JTextField;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.StructureAlignmentFactory;
import org.biojava.bio.structure.align.ce.CeCPMain;
import org.biojava.bio.structure.align.ce.CeMain;
import org.biojava.bio.structure.align.client.FarmJobParameters;
import org.biojava.bio.structure.align.client.JFatCatClient;
import org.biojava.bio.structure.align.client.PdbPair;
import org.biojava.bio.structure.align.fatcat.FatCatFlexible;
import org.biojava.bio.structure.align.fatcat.FatCatRigid;
import org.biojava.bio.structure.align.gui.AlignmentGui;
import org.biojava.bio.structure.align.gui.ChooseDirAction;
import org.biojava.bio.structure.align.gui.DisplayAFP;
import org.biojava.bio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.seq.SmithWaterman3Daligner;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.UserConfiguration;
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
				arg0.equals("ce_cp") ||
				arg0.equals("sw")
		)){
			JOptionPane.showMessageDialog(null,
					"Wrong arguments. First argument has to be \"fatcat\", \"ce\", \"ce_cp\", \"sw\", \"fatcat_flexible\", or \"biojava\", but got " + arg0);
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

			PdbPair pair = new PdbPair(name1, name2);
			System.out.println("### user provided: " + pair);

			UserConfiguration config = getWebStartConfig();

			System.setProperty("PDB_DIR",config.getPdbFilePath());
			System.out.println("using PDB file path: " + config.getPdbFilePath());

			AtomCache cache = new AtomCache(config);

			JFrame frame = new JFrame();
			frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

			showProgressBar(frame,"Loading PDB files...", "Loading files from local directory or via FTP.");

			Atom[] ca1 = cache.getAtoms(pair.getName1());
			Atom[] ca2 = cache.getAtoms(pair.getName2());

			frame.dispose();

			System.out.println("done reading structures");


			if (arg0.equalsIgnoreCase("ce") || 
					arg0.equalsIgnoreCase("ce_cp") ||
					arg0.equalsIgnoreCase("sw") ||
					arg0.equalsIgnoreCase("fatcat") ||
					arg0.equalsIgnoreCase("fatcat_flexible")){
				try {

					StructureAlignment algorithm ;
					if ( arg0.equalsIgnoreCase("ce"))
						algorithm = StructureAlignmentFactory.getAlgorithm(CeMain.algorithmName);
					else if ( arg0.equalsIgnoreCase("ce_cp"))
						algorithm = StructureAlignmentFactory.getAlgorithm(CeCPMain.algorithmName);
					else if ( arg0.equalsIgnoreCase("fatcat"))
						algorithm = StructureAlignmentFactory.getAlgorithm(FatCatRigid.algorithmName);
					else if ( arg0.equalsIgnoreCase("fatcat_flexible"))
						algorithm = StructureAlignmentFactory.getAlgorithm(FatCatFlexible.algorithmName);
					else
						algorithm = new SmithWaterman3Daligner();

					showStructureAlignment(serverLocation,algorithm ,ca1,ca2,pair.getName1(),pair.getName2());

				} catch (Exception e){
					e.printStackTrace();
					JOptionPane.showMessageDialog(null,
							"Something went wrong! : " + e.getMessage());
				}

			
			} else if ( arg0.equalsIgnoreCase("biojava")){
				try {
					//showBiojava(ca1,ca2);
				} catch (Exception e){
					e.printStackTrace();
					JOptionPane.showMessageDialog(null,
							"Something went wrong! : " + e.getMessage());
					System.exit(0);
				}
			}


		} catch (Exception e) {
			e.printStackTrace();
			JOptionPane.showMessageDialog(null,
					"Error: " + e.getMessage());
			System.exit(0);
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



	private static void showStructureAlignment(String serverLocation, StructureAlignment algorithm, Atom[] ca1, Atom[] ca2, String name1, String name2) throws StructureException{
		JFrame tmpFrame = new JFrame();
		tmpFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		String title = "Calculating " + algorithm.getAlgorithmName() + " V." + algorithm.getVersion()+" alignment... ";


		showProgressBar(tmpFrame,title, "Calculating the structure alignment.");


		//do the actual alignment 
		AFPChain afpChain = null;

		try {
			// using 10 sec as timeout on server now, since we expect the server to be able to complete the calculation within that time...
			afpChain =  JFatCatClient.getAFPChainFromServer(serverLocation,algorithm.getAlgorithmName(), name1, name2, ca1, ca2, 10000);
		} catch (Exception e){
			e.printStackTrace();
		}

		if ( afpChain == null )  {         
			afpChain = algorithm.align(ca1, ca2);
		}

		afpChain.setName1(name1);
		afpChain.setName2(name2);

		tmpFrame.dispose();


		// show results
		StructureAlignmentJmol jmol =  StructureAlignmentDisplay.display(afpChain,ca1,ca2);

		System.out.println(afpChain.toCE(ca1, ca2));

		DisplayAFP.showAlignmentImage(afpChain,ca1,ca2,jmol);

	}


//	private static void showBiojava(Atom[] ca1, Atom[] ca2) throws StructureException{
//
//		StructurePairAligner aligner = new StructurePairAligner();
//
//		StrucAligParameters params = StrucAligParameters.getDefaultParameters();
//		//params.setJoinFast(true);
//		//params.setMaxIter(1);
//		try {
//			// align the full 2 structures with default parameters.
//			// see StructurePairAligner for more options and how to align
//			// any set of Atoms
//			long start = System.currentTimeMillis();
//
//
//			aligner.align(ca1,ca2,params);
//			long end = System.currentTimeMillis();
//			System.out.println("calculation time:" + (end-start));
//
//			AlternativeAlignment[] aligs = aligner.getAlignments();
//			AlternativeAlignment a = aligs[0];
//			System.out.println(a);
//
//			if (! BiojavaJmol.jmolInClassPath()){
//				System.err.println("Could not find Jmol in classpath, please install first!");
//				return;
//			}
//			// display the alignment in Jmol
//
//
//
//			// first get an artificial structure for the alignment
//			Structure artificial = a.getAlignedStructure(structure1, structure2);
//
//			// and then send it to Jmol (only will work if Jmol is in the Classpath)
//			BiojavaJmol jmol = new BiojavaJmol();
//			jmol.setTitle(artificial.getName());
//			jmol.setStructure(artificial);
//
//			// color the two structures
//
//
//			jmol.evalString("select *; backbone 0.4; wireframe off; spacefill off; " +
//			"select not protein and not solvent; spacefill on;");
//			jmol.evalString("select */1 ; color red; model 1; ");
//
//
//			// now color the equivalent residues ...
//
//			String[] pdb1 = a.getPDBresnum1();
//			for (String res : pdb1 ){
//				jmol.evalString("select " + res + "/1 ; backbone 0.6; color orange;");
//			}
//
//			jmol.evalString("select */2; color blue; model 2;");
//			String[] pdb2 = a.getPDBresnum2();
//			for (String res :pdb2 ){
//				jmol.evalString("select " + res + "/2 ; backbone 0.6; color cyan;");
//			}
//
//
//			// now show both models again.
//			jmol.evalString("model 0;");
//
//		} catch (StructureException e){
//			e.printStackTrace();
//		}
//
//	}
}
