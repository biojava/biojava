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


import java.awt.Color;

import javax.swing.JEditorPane;

public class AlignmentTextPanel extends JPrintPanel {
	/**
	 * 
	 */
	private static final long serialVersionUID = 5092386365924879073L;
	JEditorPane tp;
	
	public AlignmentTextPanel(){
		super();
		String html = "<html><body><pre></pre></body></html>";

		tp = new JEditorPane("text/html", html);
		tp.setEditable(false);
			
		//this.setBorder(null);
		this.setBackground(Color.white);
		this.add(tp);

	}
	
	public void setText(String result){
		String html = "<html><body><pre>"+result+"</pre></body></html>";

		tp.setText(html);
	}
}
