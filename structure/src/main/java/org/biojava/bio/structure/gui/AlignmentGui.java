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
package org.biojava.bio.structure.gui;

import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.logging.Logger;
import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JProgressBar;
import javax.swing.JTabbedPane;


import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.gui.util.AlignmentCalc;
import org.biojava.bio.structure.gui.util.PDBDirPanel;
import org.biojava.bio.structure.gui.util.PDBServerPanel;
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


	JButton abortB;

	PDBDirPanel    tab1 ;
	PDBUploadPanel tab2;
	PDBServerPanel tab3;

	Thread thread;
	AlignmentCalc alicalc;
	JTabbedPane tabPane;
	JProgressBar progress;


	public static void main(String[] args){
		new AlignmentGui();

	}

	public AlignmentGui() {
		super();

		thread = null;


		this.setTitle("Pairwise Structure Alignment");
		this.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent evt) {
				JFrame frame = (JFrame) evt.getSource();
				frame.setVisible(false);
				frame.dispose();
			}
		});


		tab3 = new PDBServerPanel();


		tab1 = new PDBDirPanel();

		tab2 = new PDBUploadPanel();

		tabPane = new JTabbedPane();

		tabPane.addTab("Fetch files from FTP server",null,tab3,"fetch files from remote");

		tabPane.addTab("From local directory", null, tab1,
				"find files in a local directory");

		tabPane.addTab("Upload files",null, tab2,"Upload files manually");

		Action action1 = new AbstractAction("Submit") {
			public static final long serialVersionUID = 0l;
			// This method is called when the button is pressed
			public void actionPerformed(ActionEvent evt) {
				// Perform action...
				//System.out.println("calc structure alignment");
				calcAlignment();

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

		Action action2 = new AbstractAction("Close") {
			public static final long serialVersionUID = 0l;
			// This method is called when the button is pressed
			public void actionPerformed(ActionEvent evt) {
				// Perform action...
				abortCalc();
				dispose();
			}
		};
		JButton closeB = new JButton(action2);



		Box vBox = Box.createVerticalBox();

		Box hBox42 = Box.createHorizontalBox();
		progress =new JProgressBar();
		progress.setIndeterminate(false);

		hBox42.add(Box.createGlue());
		hBox42.add(progress);


		vBox.add(hBox42);
		vBox.add(tabPane);

		Box hBox = Box.createHorizontalBox();
		hBox.add(submitB);
		hBox.add(abortB);
		hBox.add(closeB);
		vBox.add(hBox);
		this.getContentPane().add(vBox);
		this.pack();
		this.setVisible(true);
	}




	public void cleanUp() {

		if ( alicalc != null) {
			alicalc.cleanup();
		}
	}



	private void calcAlignment() {

		//TODO: see which panel is active and get the structures ...
		int pos = tabPane.getSelectedIndex();
		StructurePairSelector tab = null;
		if (pos == 0) {
			tab = tab3;
		} else if (pos == 1 ){
			tab = tab1;
		} else if (pos == 2){
			tab = tab2;
		}
		try {
			Structure s1 = tab.getStructure1();
			Structure s2 = tab.getStructure2();

			if ( s1 == null) {
				System.err.println("please define structure 1");
				return ;
			}

			if ( s2 == null) {
				System.err.println("please define structure 2");
				return;
			}

			alicalc = new AlignmentCalc(this,s1,s2);


			thread = new Thread(alicalc);
			thread.start();
			abortB.setEnabled(true);
			progress.setIndeterminate(true);
			ProgressThreadDrawer drawer = new ProgressThreadDrawer(progress);
			drawer.start();
		} catch (Exception e){
			e.printStackTrace();
		}

	}

	public void notifyCalcFinished(){
		abortB.setEnabled(false);
		thread = null;
		progress.setIndeterminate(false);
		this.repaint();
	}

	private void abortCalc(){
		if ( alicalc != null )
			alicalc.interrupt();

	}



}

class ProgressThreadDrawer extends Thread {

	JProgressBar progress;
	static int interval = 100;

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
