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
 */
package org.biojava.nbio.structure.align.gui;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.quaternary.BioAssemblyTools;

import javax.swing.*;
import java.awt.*;

public class StructureLoaderThread extends SwingWorker<String, Object> {

	String name;
	boolean showBiolAssembly;

	UserConfiguration config;

	StructureLoaderThread(UserConfiguration config, String name, boolean showBiolAssembly){
		this.name = name;
		this.showBiolAssembly = showBiolAssembly;

		this.config = config;
	}

	@Override
	protected String doInBackground() {

		System.out.println("loading " + name );

		AtomCache cache = new AtomCache(config.getPdbFilePath(),config.getCacheFilePath());
		Structure s = null;
		try {
			if ( showBiolAssembly) {
				s= StructureIO.getBiologicalAssembly(name);

				int atomCount = StructureTools.getNrAtoms(s);

				if ( atomCount > 200000){
					// uh oh, we are probably going to exceed 512 MB usage...
					// scale down to something smaller
					System.err.println("Structure very large. Reducing display to C alpha atoms only");
					s = BioAssemblyTools.getReducedStructure(s);
				}

			} else {
				s = cache.getStructure(name);
			}

			System.out.println("done loading structure...");


			StructureAlignmentJmol jmol = new StructureAlignmentJmol();
			jmol.setStructure(s);

			jmol.evalString("set antialiasDisplay on; select all;spacefill off; wireframe off; backbone off; cartoon;color cartoon chain; select ligand;wireframe 0.16;spacefill 0.5; select all; color cartoon structure;");
			jmol.evalString("save STATE state_1");



		} catch (Exception e){
			e.printStackTrace();

			JOptionPane.showMessageDialog(null, "Error while loading " + name + ":" + e.getMessage());

			s = new StructureImpl();
		}
		hideProgressBar();
		return "Done.";

	}

	public static void showProgressBar() {

		if ( progressFrame == null){

			SwingUtilities.invokeLater(new Runnable() {

				@Override
				public void run() {
					// TODO Auto-generated method stub


					final JFrame frame = new JFrame("Loading ...");
					final JProgressBar progressBar = new JProgressBar();

					progressBar.setIndeterminate(true);

					final JPanel contentPane = new JPanel();
					contentPane.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
					contentPane.setLayout(new BorderLayout());
					contentPane.add(new JLabel("Loading ..."), BorderLayout.NORTH);
					contentPane.add(progressBar, BorderLayout.CENTER);
					frame.setContentPane(contentPane);
					frame.pack();
					frame.setLocationRelativeTo(null);
					progressFrame = frame;
					frame.setVisible(true);
				}
			});

		}


	}

	static JFrame progressFrame = null;
	private void hideProgressBar() {

		SwingUtilities.invokeLater(new Runnable() {

			@Override
			public void run() {
				if ( progressFrame != null){
					progressFrame.dispose();
					progressFrame = null;
				}
			}
		});

	}
}
