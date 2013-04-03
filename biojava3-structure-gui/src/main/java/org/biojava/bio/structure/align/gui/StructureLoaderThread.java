package org.biojava.bio.structure.align.gui;

import java.awt.BorderLayout;

import javax.swing.BorderFactory;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;

import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.align.util.UserConfiguration;
import org.biojava.bio.structure.quaternary.BioAssemblyTools;
import org.biojava3.structure.StructureIO;

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

		AtomCache cache = new AtomCache(config.getPdbFilePath(),config.isSplit());
		Structure s = null;
		try {
			if ( showBiolAssembly) {
				s= StructureIO.getBiologicalAssembly(name);

				int atomCount = StructureTools.getNrAtoms(s);

				if ( atomCount > 200000){
					// uh oh, we are probably going to exceed 512 MB usage...
					// scale down to something smaller
					System.err.println("Structure very large. Reducing display to C alpha atoms only");
					s = BioAssemblyTools.getReducedCAStructure(s);
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
