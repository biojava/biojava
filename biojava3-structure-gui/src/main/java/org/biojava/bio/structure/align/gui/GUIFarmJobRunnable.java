package org.biojava.bio.structure.align.gui;

import java.awt.Dimension;

import javax.swing.JComponent;
import javax.swing.JFrame;

import org.biojava.bio.structure.align.FarmJob;
import org.biojava.bio.structure.align.client.FarmJobParameters;

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
	        public void run() {
	            createAndShowGUI(progressListener);
	        }
		});
		
	}	

}
