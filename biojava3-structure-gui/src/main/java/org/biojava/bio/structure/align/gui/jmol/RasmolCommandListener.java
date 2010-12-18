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
 * created at Oct 27, 2007
 */
package org.biojava.bio.structure.align.gui.jmol;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JTextField;

import org.biojava.bio.structure.align.gui.jmol.JmolPanel;
import org.biojava.bio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.bio.structure.gui.BiojavaJmol;




/** a utility class that listens to Ramsol script commands in the @link {@link BiojavaJmol} class
 * 
 * @author Andreas Prlic
 *
 */
public class RasmolCommandListener 
extends KeyAdapter
implements ActionListener,
MouseListener {

	JTextField textfield;
	JmolPanel jmolPanel;

	List<String> history;
	int historyPosition;

	public RasmolCommandListener(JmolPanel panel, JTextField field){
		textfield = field;
		jmolPanel = panel;
		history = new ArrayList<String>();
		historyPosition = -2; // -2 = history = empty;
	}	
	
	public void actionPerformed(ActionEvent event) {
		/*
	        if ( spice.isLoading() ) {
	            logger.finest("loading data, please be patient");
	            return ;
	        }
		 */
		String cmd = textfield.getText();
		jmolPanel.executeCmd(cmd);
		textfield.setText("");

		// now comes history part:

		// no need for history:
		if ( cmd.equals("")) return;

		// check last command in history
		// if equivalent, don't add,
		// otherwise add               
		if (history.size()>0){
			String txt=(String)history.get(history.size()-1);
			if (! txt.equals(cmd)) {
				history.add(cmd);  
			}
		} else {             
			// the first time always add
			history.add(cmd);
		}
		historyPosition=history.size();


	}

	public void  mouseClicked(MouseEvent e){
		String cmd = textfield.getText();
		if ( cmd.equals(StructureAlignmentJmol.COMMAND_LINE_HELP)){
			textfield.setText("");
			textfield.repaint();
		}
	};


	public void  mouseExited(MouseEvent e){};
	public void  mouseReleased(MouseEvent e){};
	public void  mousePressed(MouseEvent e){};

	public void  mouseEntered(MouseEvent e){};

	/** takes care of the cursor up/down keys. triggers copying of stored 
	 * commands into the current textfield
	 * 
	 */ 
	 

	public void keyReleased(KeyEvent e){

		int code = e.getKeyCode();
		//String s = e.getKeyText(code);
		//System.out.println(s);
		if (( code == KeyEvent.VK_UP ) || 
				( code == KeyEvent.VK_KP_UP)) {
			// go one back in history;
			if ( historyPosition > 0){
				historyPosition= historyPosition-1;              
			} 
		} else if (( code == KeyEvent.VK_DOWN ) || 
				( code == KeyEvent.VK_KP_DOWN)) {            
			if ( historyPosition < (history.size()-1) ){
				historyPosition++;                
			} else {
				// clear command if at beginning of history
				textfield.setText("");
				historyPosition=history.size();
				return;
			}
		} else if ( code == KeyEvent.VK_PAGE_UP) {
			if ( historyPosition > 0) {
				historyPosition = 0;
			}
		} else if ( code == KeyEvent.VK_PAGE_DOWN) {
			if ( historyPosition >= 0) {
				historyPosition = history.size()-1;
			}
		} else {
			// some other key has been pressed, do nothing
			return;
		}

		if ( historyPosition >= 0) {
			String txt = (String)history.get(historyPosition);
			textfield.setText(txt);
		}


	}


}
