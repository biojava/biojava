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

import org.biojava.nbio.structure.align.FarmJob;
import org.biojava.nbio.structure.align.client.FarmJobParameters;

import javax.swing.*;
import java.awt.*;

public class GUIFarmJobRunnable implements Runnable{
	FarmJobParameters params;
	GUIAlignmentProgressListener progressListener ;
	public GUIFarmJobRunnable(FarmJobParameters params){
		this.params = params;


	}

	/**
	 * Create the GUI and show it. As with all GUI code, this must run
	 * on the event-dispatching thread.
	 */
	private static void createAndShowGUI(GUIAlignmentProgressListener progressListener) {
		//Create and set up the window.
		JFrame frame = new JFrame("Monitor alignment process");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		//Create and set up the content pane.
		JComponent newContentPane = progressListener;
		newContentPane.setOpaque(true); //content panes must be opaque
		newContentPane.setSize(new Dimension(400,400));
		frame.setContentPane(newContentPane);

		//Display the window.
		frame.pack();
		frame.setVisible(true);
	}

	@Override
	public void run() {

		progressListener = new GUIAlignmentProgressListener();
		progressListener.logStatus(params.toString());

		//createAndShowGUI(progressListener);

		FarmJob job = new FarmJob();

		progressListener.setFarmJob(job);

		job.addAlignmentProgressListener(progressListener);
		job.setParams(params);

		Thread t = new Thread(job);
		t.start();


		javax.swing.SwingUtilities.invokeLater(new Runnable() {
	        @Override
			public void run() {
	            createAndShowGUI(progressListener);
	        }
		});

	}

}
