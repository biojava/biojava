package org.biojava.bio.structure.align.gui;

import java.awt.BorderLayout;
import java.awt.Cursor;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JButton;

import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

import org.biojava.bio.structure.align.FarmJob;
import org.biojava.bio.structure.align.events.AlignmentProgressListener;

/** a GUI that allows to watch progress as multiple alignments are being processed.
 * 
 * @author Andreas Prlic
 *
 */
public class GUIAlignmentProgressListener extends JPanel implements AlignmentProgressListener, ActionListener {




	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	int alignmentsProcessed;

	JProgressBar progressBar;
	private JTextArea taskOutput;
	private JButton stopButton;

	FarmJob farmJob;

	public GUIAlignmentProgressListener(){

		super(new BorderLayout());
		stopButton = new JButton("Stop");
		stopButton.setActionCommand("Stop");
		stopButton.addActionListener(this);

		progressBar = new JProgressBar(0, 100);
		progressBar.setValue(0);
		progressBar.setStringPainted(true);

		taskOutput = new JTextArea(5, 20);
		taskOutput.setMargin(new Insets(5,5,5,5));
		taskOutput.setEditable(false);

		JPanel panel = new JPanel();
		panel.add(stopButton);
		panel.add(progressBar);

		add(panel, BorderLayout.PAGE_START);
		add(new JScrollPane(taskOutput), BorderLayout.CENTER);
		setBorder(BorderFactory.createEmptyBorder(20, 20, 20, 20));

	}
	
	

	 public FarmJob getFarmJob() {
		return farmJob;
	}



	public void setFarmJob(FarmJob parent) {
		this.farmJob = parent;
	}



	/**
     * Invoked when the user presses the stop button.
     */
    public void actionPerformed(ActionEvent evt) {
    	
    	//System.out.println("stopping!");
    	logStatus("terminating");
    	logStatus(" Total alignments processed: " + alignmentsProcessed);
        stopButton.setEnabled(false);
        setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
        progressBar.setIndeterminate(true);
        progressBar.setStringPainted(false);
        System.out.println("terminating jobs");
        
        farmJob.terminate();
       
        setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
        progressBar.setIndeterminate(false);
        
    }

   
	public void alignmentEnded() {
		
		alignmentsProcessed++;
		
		//System.out.println("aligned " + alignmentsProcessed );
		int v = progressBar.getValue();
		
		progressBar.setValue(v+1);
		progressBar.setString(v+"");
        synchronized(this){notifyAll();}

	}

	public void alignmentStarted(String name1, String name2) {	
		logStatus("#" + progressBar.getValue() + " starting alignment of " + name1 + " " + name2);
	}

	public void downloadingStructures(String name) {		
		logStatus("Downloading " + name );
	}

	public void logStatus(String message) {		
		taskOutput.append(message+"\n");
	}

	public void requestingAlignmentsFromServer(int nrAlignments) {
		logStatus("Requesting " + nrAlignments + " alignments to be calculated");
		progressBar.setMaximum(nrAlignments);
		progressBar.setValue(0);	
        synchronized(this){notifyAll();}

	}

	public void sentResultsToServer(int nrAlignments, String serverMessage) {
		logStatus("sent alignment results back to server. Server responded: >"+serverMessage+"<");
		progressBar.setValue(0);
		synchronized(this){notifyAll();}
	}


}
