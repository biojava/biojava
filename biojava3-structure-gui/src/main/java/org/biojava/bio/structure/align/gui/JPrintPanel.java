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
 * Created on Oct 9, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.gui;

import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;

import javax.swing.JPanel;

public class JPrintPanel extends JPanel implements Printable,ActionListener{
	/**
	 * 
	 */
	private static final long serialVersionUID = -3337337068138131455L;

	public int print(Graphics g, PageFormat pf, int pi) throws PrinterException {

		if (pi >= 1) {
			return Printable.NO_SUCH_PAGE;
		}
		Graphics2D g2D = (Graphics2D) g;
		g.translate(20, 20);
		Font  f = new Font("Monospaced",Font.PLAIN,10);
		g.setFont (f);
		
		double scale = pf.getImageableWidth()/this.getSize().getWidth();

	    g2D.scale(scale,scale);
		
		paint (g);
		
		
		
	    
		
		return Printable.PAGE_EXISTS;
	}

	public void actionPerformed(ActionEvent e) {
		
		PrinterJob printJob = PrinterJob.getPrinterJob();
		printJob.setPrintable(this);

		try { 
			if(printJob.printDialog()){
				printJob.print();
			}
		} catch (Exception printException) { 
			System.err.println("Error during printing: " +printException.getMessage());
			printException.printStackTrace();
		}
	}

}