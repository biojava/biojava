// $Id: ColorSchemeChooser.java,v 1.5 2009/08/03 16:02:17 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
// Copyright (C) 2003-2007 Ethalinda K.S. Cannon
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.archaeopteryx;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;
import javax.swing.event.ListDataEvent;
import javax.swing.event.ListDataListener;

final class ColorSchemeChooser extends JDialog implements ActionListener {

    private static final long  serialVersionUID = 6150960100859081126L;
    private final TreeColorSet _colorset;
    private final JComboBox    _selector;
    private final JPanel       _color_panel;
    private final JPanel       _color_labels[];
    private final JButton      _ok_btn;
    private final JButton      _cancel_btn;
    private final MainPanel    _main_panel;
    private final int          _prev_selected_scheme;
    private int                _selected_scheme;

    ColorSchemeChooser( final MainPanel parent, final TreeColorSet colorset ) {
        setName( "Color Scheme Chooser" );
        setModal( true );
        _colorset = colorset;
        _prev_selected_scheme = _colorset.getCurrentColorScheme();
        _main_panel = parent;
        setSize( 400, 350 );
        final Container contentpane = getContentPane();
        contentpane.setLayout( new BorderLayout( 5, 15 ) );
        // The scheme selection panel
        final JPanel select_panel = new JPanel();
        final JLabel l = new JLabel( "Choose a color scheme:" );
        select_panel.add( l );
        final Vector<String> list = new Vector<String>();
        for( final String element : TreeColorSet.SCHEME_NAMES ) {
            list.add( element );
        }
        _selector = new JComboBox( list );
        _selector.setMaximumRowCount( list.size() );
        _selector.getModel().addListDataListener( new ListDataListener() {

            public void contentsChanged( final ListDataEvent e ) {
                final int selection = _selector.getSelectedIndex();
                changeDialogColors( selection );
            }

            public void intervalAdded( final ListDataEvent e ) {
                // Not needed.
            }

            public void intervalRemoved( final ListDataEvent e ) {
                // Not needed.
            }
        } );
        select_panel.add( _selector );
        contentpane.add( select_panel, "North" );
        // create color panel
        final int num_colors = TreeColorSet.COLOR_FIELDS.length;
        _color_panel = new JPanel( new GridLayout( num_colors, 2, 8, 0 ) );
        final JLabel headings[] = new JLabel[ num_colors ];
        _color_labels = new JPanel[ num_colors ];
        for( int i = 0; i < num_colors; i++ ) {
            headings[ i ] = new JLabel( TreeColorSet.COLOR_FIELDS[ i ] );
            headings[ i ].setFont( new Font( Configuration.getDefaultFontFamilyName(), Font.PLAIN, 9 ) );
            headings[ i ].setHorizontalAlignment( SwingConstants.RIGHT );
            _color_panel.add( headings[ i ] );
            _color_labels[ i ] = new JPanel();
            _color_labels[ i ].setPreferredSize( new Dimension( 15, 40 ) );
            _color_panel.add( _color_labels[ i ] );
        }
        contentpane.add( _color_panel, "Center" );
        setColors( _colorset.getColorSchemes()[ 0 ] );
        // create button panel
        final JPanel btn_panel = new JPanel();
        _ok_btn = new JButton( "OK" );
        _ok_btn.addActionListener( new ActionListener() {

            public void actionPerformed( final ActionEvent e ) {
                ok();
            }
        } );
        btn_panel.add( _ok_btn );
        _cancel_btn = new JButton( "Cancel" );
        _cancel_btn.addActionListener( new ActionListener() {

            public void actionPerformed( final ActionEvent e ) {
                cancel();
            }
        } );
        btn_panel.add( _cancel_btn );
        btn_panel.setPreferredSize( new Dimension( 400, 30 ) );
        getContentPane().add( btn_panel, "South" );
        setCurrentColor( colorset.getCurrentColorScheme() );
    }

    public void actionPerformed( final ActionEvent e ) {
        // Not needed.
    }

    private void cancel() {
        _colorset.setColorSchema( _prev_selected_scheme );
        for( final TreePanel tree_panel : getMainPanel().getTreePanels() ) {
            tree_panel.setBackground( _colorset.getBackgroundColor() );
        }
        redrawTreePanel();
        setVisible( false );
        dispose();
    }

    private void changeDialogColors( final int scheme_index ) {
        _selected_scheme = scheme_index;
        setColors( _colorset.getColorSchemes()[ scheme_index ] );
        _colorset.setColorSchema( getSelectedScheme() );
        for( final TreePanel tree_panel : getMainPanel().getTreePanels() ) {
            tree_panel.setBackground( _colorset.getBackgroundColor() );
        }
        redrawTreePanel();
    }

    private MainPanel getMainPanel() {
        return _main_panel;
    }

    private int getSelectedScheme() {
        return _selected_scheme;
    }

    private void ok() {
        // set the new color
        _colorset.setColorSchema( getSelectedScheme() );
        // close the window
        setVisible( false );
        dispose();
    }

    private void redrawTreePanel() {
        if ( getMainPanel().getCurrentTreePanel() != null ) {
            getMainPanel().getCurrentTreePanel().repaint();
        }
    }

    private void setColors( final Color colors[] ) {
        for( int i = 0; i < colors.length; i++ ) {
            _color_labels[ i ].setBackground( colors[ i ] );
        }
    }

    private void setCurrentColor( final int color_index ) {
        setColors( _colorset.getColorSchemes()[ color_index ] );
        _selector.setSelectedIndex( color_index );
    }
}
