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
 * Created on Jul 17, 2006
 *
 */
package org.biojava.nbio.structure.gui.util;

import javax.swing.*;
import javax.swing.table.TableColumnModel;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

class JTableMouseButtonListener implements MouseListener {

	private JTable __table;

	private void __forwardEventToButton(MouseEvent e) {
		TableColumnModel columnModel = __table.getColumnModel();
		int column = columnModel.getColumnIndexAtX(e.getX());
		int row    = e.getY() / __table.getRowHeight();

		Component component;


		//System.out.println("row " + row + " col " + column);
		if(row >= __table.getRowCount() || row < 0 ||
				column >= __table.getColumnCount() || column < 0)
			return;

		Object value = __table.getValueAt(row, column);

		if(!(value instanceof Component))
			return;

		//System.out.println("converting event");
		component = (Component)value;

		MouseEvent mevent = SwingUtilities.convertMouseEvent(__table, e, component);

		//System.out.println(mevent);

		component.dispatchEvent(mevent);
		// This is necessary so that when a button is pressed and released
		// it gets rendered properly.  Otherwise, the button may still appear
		// pressed down when it has been released.
		__table.repaint();
	}

	public JTableMouseButtonListener(JTable table) {
		__table = table;
	}

	@Override
	public void mouseClicked(MouseEvent e) {
		//System.out.println("clicked!");
		__forwardEventToButton(e);
	}

	@Override
	public void mouseEntered(MouseEvent e) {
		__forwardEventToButton(e);
	}

	@Override
	public void mouseExited(MouseEvent e) {
		__forwardEventToButton(e);
	}

	@Override
	public void mousePressed(MouseEvent e) {
		__forwardEventToButton(e);
	}

	@Override
	public void mouseReleased(MouseEvent e) {
		__forwardEventToButton(e);
	}
}
