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

import javax.swing.table.AbstractTableModel;

class JTableDataButtonModel extends AbstractTableModel {

		public static final long serialVersionUID = 0l;

		Object[][] __rows;
		String[]   __columns;

		public JTableDataButtonModel(Object[][] rows, String[] columns){
				__rows = rows;
				__columns = columns;
		}


		@Override
	public String getColumnName(int column) {
			return __columns[column];
		}

		@Override
	public int getRowCount() {
			return __rows.length;
		}

		@Override
	public int getColumnCount() {
			return __columns.length;
		}

		@Override
	public Object getValueAt(int row, int column) {
				return __rows[row][column];
		}

		@Override
	public boolean isCellEditable(int row, int column) {
			return false;
		}

		@Override
	@SuppressWarnings({ "unchecked" })
	public Class getColumnClass(int column) {
			return getValueAt(0, column).getClass();
		}
	}
