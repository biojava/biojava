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
package org.biojava.bio.structure.gui.util;

import java.awt.Component;

import javax.swing.JTable;
import javax.swing.table.TableCellRenderer;

class JButtonTableCellRenderer implements TableCellRenderer {
    private TableCellRenderer __defaultRenderer;

    public JButtonTableCellRenderer(TableCellRenderer renderer) {
      __defaultRenderer = renderer;
    }

    public Component getTableCellRendererComponent(JTable table, Object value,
                         boolean isSelected,
                         boolean hasFocus,
                         int row, int column)
    {
      if(value instanceof Component)
        return (Component)value;
      return __defaultRenderer.getTableCellRendererComponent(
       table, value, isSelected, hasFocus, row, column);
    }
  }